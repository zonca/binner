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
#include "Epetra_FEVbrMatrix.h"

#include "H5PlanckDataManager.h"
#include "MapWriter.h"
#include "CreateMatrices.h"
#include "ReadParameterFile.h"

#include "H5Cpp.h"

using namespace std;
using namespace H5;
using boost::format;

int WriteH5Vec(Epetra_Vector * vec, string filename) {
    int MyPID = vec->Comm().MyPID();
    H5std_string  FILE_NAME( str( format("%s_%03d.h5") % filename % MyPID ) );
    H5File file(FILE_NAME, H5F_ACC_TRUNC );
    hsize_t dimsf[1];
    dimsf[0] = vec->Map().NumMyElements();
    DataSpace dataspace( 1, dimsf );
    DataSet dataset = file.createDataSet( "Vector", PredType::NATIVE_DOUBLE, dataspace );
    double * data;
    vec->ExtractView(&data);
    dataset.write( data, PredType::NATIVE_DOUBLE );
    return 0;
}


int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif  

int MyPID = Comm.MyPID();
int i_M;

int SAMPLES_PER_PROC = 1.2 * 1e6;

Epetra_Time time(Comm);
H5PlanckDataManager* dm;

string parameterFilename = "notimplemented.dat";
readParameterFile(parameterFilename, dm);

log(MyPID, format("Number of elements [mil]: %.10d") % (dm->getDatasetLength()/1.e6));
log(MyPID, format("Elements per channel [mil]: %.10d") % (dm->getLengthPerChannel()/1.e6));
log(MyPID, format("Samples per proc [mil]: %.10d") % (SAMPLES_PER_PROC/1.e6));
//Epetra_BlockMap Map(dm->getDatasetLength(), 1, 0, Comm);
Epetra_Map Map(-1, SAMPLES_PER_PROC, 0, Comm);

int NumMyElements = Map.NumMyElements();
//cout << MyPID << " " << Map.MinMyGID() << " " << NumMyElements << endl;

Epetra_Map PixMap(dm->getNPIX(),0,Comm);
Epetra_BlockMap PixBlockMap(dm->getNPIX(),dm->NSTOKES,0,Comm);

// declaring matrices
log(MyPID,"POINTING MATRIX");
Epetra_MultiVector summap(PixMap, dm->NSTOKES);
Epetra_MultiVector M(PixMap, dm->NSTOKES * (dm->NSTOKES + 1) / 2);
Epetra_Vector tempmap(PixMap);
Epetra_Vector * hitmap = new Epetra_Vector(PixMap);
Epetra_CrsMatrix * P;
Epetra_CrsGraph * Graph; 

Epetra_MultiVector yqu = Epetra_MultiVector(Map, 3);

double ** yqu_view = new double *[3];
yqu.ExtractView(&yqu_view);

string LABEL[5] = {"I", "Q", "U", "S1", "S2"};

Epetra_IntVector pix = Epetra_IntVector(Map);

int * pix_view;
pix.ExtractView(&pix_view);

Epetra_Vector tempvec = Epetra_Vector(Map);

double weight = 0;

//LOOP
BOOST_FOREACH( string channel, dm->getChannels())
    {
        log(MyPID, format("Processing channel %s") % channel);

        weight = dm->getWeight(channel);
        log(MyPID, format("Weight %f") % weight);

        for (long offset=0; offset<dm->getLengthPerChannel(); offset=offset+Map.NumGlobalElements()) {
            log(MyPID,"READ POINTING");
            time.ResetStartTime();
            dm->getPointing(channel, Map.MinMyGID() + offset, NumMyElements, pix_view, yqu_view[1], yqu_view[2]);
            log(MyPID, format("Read data timer %f") % time.ElapsedTime());

            log(MyPID,"READ DATA");
            time.ResetStartTime();
            dm->getData(channel, Map.MinMyGID() + offset,Map.NumMyElements(),yqu_view[0]);
            log(MyPID, format("Read data timer %f") % time.ElapsedTime());


            log(MyPID,format("Offset [mil]: %d") % (offset/1.e6));
            time.ResetStartTime();
            createGraph(Map, PixMap, pix, Graph);
            log(MyPID, format("Create Graph timer %f") % time.ElapsedTime());
            time.ResetStartTime();
            P = new Epetra_CrsMatrix(Copy, *Graph);
            P->PutScalar(1.);
            log(MyPID, format("Create P timer %f") % time.ElapsedTime());


            log(MyPID,"SUM MAP");

            time.ResetStartTime();
            ////// I
            log(MyPID,"I");
            P->Multiply1(true,*(yqu(0)),tempmap); //SUMMAP = Pt y
            summap(0)->Update(1., tempmap, weight);
            log(MyPID, format("%f") % time.ElapsedTime());

            time.ResetStartTime();
            log(MyPID,"HitMap");
            tempvec.PutScalar(1.);
            P->Multiply1(true,tempvec,tempmap);
            hitmap->Update(1., tempmap, weight);
            log(MyPID, format("%f") % time.ElapsedTime());
            time.ResetStartTime();

            //// Q U
            log(MyPID,"QU");
            for (int i=1; i<3; ++i) { // Q=1 U=2
                tempvec.Multiply(1., *(yqu(0)), *(yqu(i)), 0.);
                P->Multiply1(true,tempvec,tempmap); //SUMMAP = Pt y
                summap(i)->Update(1., tempmap, weight);
            }
            log(MyPID, format("%f") % time.ElapsedTime());

            log(MyPID, "Create M");
            i_M = 0;
            for (int j=0; j<dm->NSTOKES; ++j) {
                for (int k=j; k<dm->NSTOKES; ++k) {
                    log(MyPID,format("M %d %d") % j % k);
                    time.ResetStartTime();
                    if (k == 0) { //also j=0
                        tempvec.PutScalar(1.);
                    } else if (j == 0 ) {
                        tempvec.Update(1., *(yqu(k)), 0.);
                    } else {
                        tempvec.Multiply(1., *(yqu(j)), *(yqu(k)), 0.);
                    }
                    P->Multiply1(true,tempvec,tempmap); //SUMMAP = Pt y
                    log(MyPID, format("Setting M %d") % i_M );
                    M(i_M)->Update(1., tempmap, weight);
                    log(MyPID, format("M %d %d: %f") % j % k % time.ElapsedTime());
                    i_M++;
                }
            } // M loop
        } // chunck loop
    } // channel loop
//end LOOP
//
for (int j=0; j<dm->NSTOKES; ++j) {
    WriteH5Vec(summap(j), "summap_" + LABEL[j]);
}

WriteH5Vec(hitmap, "hitmap");

log(MyPID,"Computing RCOND and Inverting");
Epetra_Vector rcond(PixMap);
Epetra_Vector binmap(PixMap);
invertM(PixMap, M, rcond, summap);

log(MyPID,"Writing MAPS");

for (int j=0; j<dm->NSTOKES; ++j) {
    WriteH5Vec(summap(j), "binmap_" + LABEL[j]);
}

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
return(0);

};
