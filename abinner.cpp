#include <string>
#include <math.h>

#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Time.h"
#include "Teuchos_CommandLineProcessor.hpp"

#include "H5PlanckDataManager.h"
#include "CreateMatrices.h"
#include "ReadParameterFile.h"

#include "H5Cpp.h"

using namespace std;
using namespace H5;
using boost::format;

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif  

Teuchos::CommandLineProcessor My_CLP;
string parameterFilename = "config.xml";
My_CLP.setOption("p", &parameterFilename, "Parameter filename");
My_CLP.parse( argc, argv );

int MyPID = Comm.MyPID();
int i_M, a;

int SAMPLES_PER_PROC;
int MAX_SAMPLES_PER_PROC = 5e6;
int LOOPS = 1;

Epetra_Time time(Comm);
H5PlanckDataManager* dm;

log(MyPID, parameterFilename);
readParameterFile(parameterFilename, dm);

if (dm->DEBUG) {
    SAMPLES_PER_PROC = 5;
}

log(MyPID, format("Number of elements [mil]: %.10d") % (dm->getDatasetLength()/1.e6));
log(MyPID, format("Elements per channel [mil]: %.10d") % (dm->getLengthPerChannel()/1.e6));
do  {
    SAMPLES_PER_PROC = dm->getLengthPerChannel()/Comm.NumProc()/LOOPS+1;
    log(MyPID, format("%d loops: Samples per proc [mil]: %.10d") % LOOPS % (SAMPLES_PER_PROC/1.e6));
    LOOPS++;
} while (SAMPLES_PER_PROC > MAX_SAMPLES_PER_PROC);

log(MyPID, format("Samples per proc [mil]: %.10d") % (SAMPLES_PER_PROC/1.e6));
//Epetra_BlockMap Map(dm->getDatasetLength(), 1, 0, Comm);
Epetra_Map Map(-1, SAMPLES_PER_PROC, 0, Comm);

int NumMyElements = Map.NumMyElements();
//cout << MyPID << " " << Map.MinMyGID() << " " << NumMyElements << endl;

Epetra_Map PixMap(dm->getNPIX(),0,Comm);
Epetra_BlockMap PixBlockMap(dm->getNPIX(),dm->NSTOKES,0,Comm);

// declaring matrices
log(MyPID,"POINTING MATRIX");
Epetra_MultiVector summap(PixMap, dm->NSTOKES*dm->getChannels().size()/2.);
Epetra_MultiVector M(PixMap, dm->NSTOKES * (dm->NSTOKES + 1) / 2);
Epetra_Vector tempmap(PixMap);
Epetra_Vector * hitmap = new Epetra_Vector(PixMap);
Epetra_CrsMatrix * P;
Epetra_CrsGraph * Graph; 

Epetra_MultiVector yqu = Epetra_MultiVector(Map, 3);

double ** yqu_view = new double *[3];
yqu.ExtractView(&yqu_view);

string LABEL[9] = {"I", "Q", "U", "S1", "S2", "S3", "S4", "S5", "S6"};

Epetra_IntVector pix = Epetra_IntVector(Map);

int * pix_view;
pix.ExtractView(&pix_view);

Epetra_Vector tempvec = Epetra_Vector(Map);

double weight = 0;
int ch_number = 0;

//LOOP
BOOST_FOREACH( string channel, dm->getChannels())
    {
        log(MyPID, format("//////////////////Processing channel %s") % channel);

        weight = dm->getWeight(channel);
        log(MyPID, format("Weight %f") % weight);

        for (long offset=0; offset<dm->getLengthPerChannel(); offset=offset+Map.NumGlobalElements()) {
            log(MyPID,"READ POINTING");
            time.ResetStartTime();
            dm->getPointing(channel, Map.MinMyGID() + offset, NumMyElements, pix_view, yqu_view[1], yqu_view[2]);
            log(MyPID, format("Read data timer %f") % time.ElapsedTime());

            log(MyPID,"READ DATA");
            time.ResetStartTime();
            //dm->getData(channel, Map.MinMyGID() + offset,Map.NumMyElements(),yqu_view[0]);
            a = 2 * (ch_number % 2) - 1; // -1 for M, +1 for S
            for (int i=0; i<NumMyElements; i++) {
                yqu[0][i] = a * (1 + yqu[1][i] + yqu[2][i]); // rescanning a +1/-1 map
            }
            log(MyPID, format("Read data timer %f") % time.ElapsedTime());

            log(MyPID,format("Offset [mil]: %d") % (offset/1.e6));
            time.ResetStartTime();
            createGraph(Map, PixMap, pix, Graph);
            log(MyPID, format("Create Graph timer %f") % time.ElapsedTime());
            time.ResetStartTime();
            P = new Epetra_CrsMatrix(Copy, *Graph);
            P->PutScalar(1.);
            P->FillComplete(PixMap, Map);
            log(MyPID, format("Create P timer %f") % time.ElapsedTime());

            log(MyPID,"SUM MAP");

            time.ResetStartTime();
            ////// I
            log(MyPID,"I");
            P->Multiply1(true,*(yqu(0)),tempmap); //SUMMAP = Pt y
            summap(3*(ch_number/2))->Update(weight, tempmap, 1.);
            log(MyPID, format("%f") % time.ElapsedTime());

            time.ResetStartTime();
            log(MyPID,"HitMap");
            tempvec.PutScalar(1.);
            P->Multiply1(true,tempvec,tempmap);
            hitmap->Update(1., tempmap, 1.);
            log(MyPID, format("%f") % time.ElapsedTime());
            time.ResetStartTime();


            //// Q U
            log(MyPID,"QU");
            for (int i=1; i<3; ++i) { // Q=1 U=2
                tempvec.Multiply(1., *(yqu(0)), *(yqu(i)), 0.);
                P->Multiply1(true,tempvec,tempmap); //SUMMAP = Pt y
                summap(3*(ch_number/2)+i)->Update(weight, tempmap, 1.);
            }
            log(MyPID, format("%f") % time.ElapsedTime());

            log(MyPID, "Create M");
            log(MyPID, "IQU loop");
            time.ResetStartTime();
            for (int j=0; j<3; ++j) {
                for (int k=j; k<3; ++k) {
                    log(MyPID,format("M %d %d") % j % k);
                    if (k == 0) { //also j=0
                        tempvec.PutScalar(1.);
                    } else if (j == 0 ) {
                        tempvec.Update(1., *(yqu(k)), 0.);
                    } else {
                        tempvec.Multiply(1., *(yqu(j)), *(yqu(k)), 0.);
                    }
                    P->Multiply1(true,tempvec,tempmap); //SUMMAP = Pt y
                    i_M = j * (2*dm->NSTOKES-1 - j)/2 + k;
                    log(MyPID, format("Setting M %d") % i_M );
                    M(i_M)->Update(weight, tempmap, 1.);
                }
            } // M loop IQU
            log(MyPID, format("IQU loop: %f") % time.ElapsedTime());

            delete P;
            delete Graph;
        } // chunck loop
    ch_number++;
    } // channel loop
//end LOOP

log(MyPID,"Writing SUMMAP");

time.ResetStartTime();
for (int j=0; j<dm->NSTOKES; ++j) {
    WriteH5Vec(summap(j), "summap_" + LABEL[j], dm->outputFolder);
}
log(MyPID, format("%f") % time.ElapsedTime());

WriteH5Vec(hitmap, "hitmap", dm->outputFolder);

log(MyPID,"Computing RCOND and Inverting");
time.ResetStartTime();
Epetra_Vector rcond(PixMap);

Epetra_Vector binmap(PixMap);
invertM(PixMap, dm, M, rcond, summap);
log(MyPID, format("%f") % time.ElapsedTime());

WriteH5Vec(&rcond, "rcondmap", dm->outputFolder);

log(MyPID,"Writing BINMAP");
time.ResetStartTime();

for (int j=0; j<dm->NSTOKES; ++j) {
    for (int ch_number=0; ch_number<dm->getChannels().size(); ch_number+=2) {
        WriteH5Vec(summap(3*(ch_number/2) + j), "binmap_" + dm->getChannels()[ch_number] + "_" + LABEL[j], dm->outputFolder);
    }
}
log(MyPID, format("%f") % time.ElapsedTime());

log(MyPID,"Writing M");
time.ResetStartTime();
for (int j=0; j<dm->NSTOKES; ++j) {
    for (int k=j; k<dm->NSTOKES; ++k) {
        i_M = dm->getIndexM(j, k);
        WriteH5Vec(M(i_M), "M_" + LABEL[j] + "_" + LABEL[k], dm->outputFolder);
    }
}
log(MyPID, format("%f") % time.ElapsedTime());

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
return(0);

};
