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

int WriteH5Vec(Epetra_Vector& vec, string filename) {
    int MyPID = vec.Comm().MyPID();
    H5std_string  FILE_NAME( str( format("%s_%03d.h5") % filename % MyPID ) );
    H5File file(FILE_NAME, H5F_ACC_TRUNC );
    hsize_t dimsf[1];
    dimsf[0] = vec.Map().NumMyElements();
    DataSpace dataspace( 1, dimsf );
    DataSet dataset = file.createDataSet( "Vector", PredType::NATIVE_DOUBLE, dataspace );
    double * data;
    vec.ExtractView(&data);
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

int SAMPLES_PER_PROC = 2. * 1e6;

Epetra_Time time(Comm);
H5PlanckDataManager* dm;

string parameterFilename = "notimplemented.dat";
readParameterFile(parameterFilename, dm);

log(Comm.MyPID(), format("Number of elements [mil]: %.10d") % (dm->getDatasetLength()/1.e6));
log(Comm.MyPID(), format("Elements per channel [mil]: %.10d") % (dm->getLengthPerChannel()/1.e6));
//Epetra_BlockMap Map(dm->getDatasetLength(), 1, 0, Comm);
Epetra_Map Map(-1, SAMPLES_PER_PROC, 0, Comm);

int NumMyElements = Map.NumMyElements();
cout << Comm.MyPID() << " " << Map.MinMyGID() << " " << NumMyElements << endl;

Epetra_Map PixMap(dm->getNPIX(),0,Comm);

// declaring matrices
log(Comm.MyPID(),"POINTING MATRIX");
Epetra_MultiVector summap(PixMap, dm->NSTOKES);
Epetra_Vector temp_summap(PixMap);
Epetra_CrsMatrix * Ptemp;
Epetra_CrsMatrix * PQ;
Epetra_CrsMatrix * PU;
Epetra_CrsGraph * Graph; 
Epetra_Vector * y;
double* data; 

Epetra_FEVbrMatrix invM(Copy, PixMap, 1);

// initialize M
log(Comm.MyPID(),"Initializing M");
initM(PixMap, dm->NSTOKES, invM);
log(Comm.MyPID(),"GlobalAssemble");
invM.GlobalAssemble();
log(Comm.MyPID(),"GlobalAssemble DONE");
cout << time.ElapsedTime() << endl;

//data = new double[Map.NumMyElements()];
y = new Epetra_Vector(Map);
y->ExtractView(&data);

double weight = 0;

//LOOP
BOOST_FOREACH( string channel, dm->getChannels())
    {
        log(Comm.MyPID(), format("Processing channel %s") % channel);

        weight = dm->getWeight(channel);
        log(Comm.MyPID(), format("Weight %f") % weight);

        for (long offset=0; offset<dm->getLengthPerChannel(); offset=offset+Map.NumGlobalElements()) {

            log(Comm.MyPID(),format("Offset [mil]: %d") % (offset/1.e6));
            createP(channel, Map, PixMap, dm, offset, PQ, PU, Graph);

            log(Comm.MyPID(),"READ DATA");
            dm->getData(channel, Map.MinMyGID() + offset,Map.NumMyElements(),data);

            cout << time.ElapsedTime() << endl;

            log(Comm.MyPID(),"SUM MAP");
            Ptemp = new Epetra_CrsMatrix(Copy, *Graph);

            ////// I
            log(Comm.MyPID(),"I");
            Ptemp->PutScalar(1.);
            Ptemp->Multiply1(true,*y,temp_summap); //SUMMAP = Pt y
            for (int i=0; i<PixMap.NumMyElements(); i++) {
                (*(summap(0)))[i] += weight * temp_summap[i];
            }
            delete Ptemp;

            //// Q
            log(Comm.MyPID(),"Q");
            PQ->Multiply1(true,*y,temp_summap); //SUMMAP = Pt y
            for (int i=0; i<PixMap.NumMyElements(); i++) {
                (*(summap(1)))[i] += weight * temp_summap[i];
            }

            //// U
            log(Comm.MyPID(),"U");
            PQ->Multiply1(true,*y,temp_summap); //SUMMAP = Pt y
            for (int i=0; i<PixMap.NumMyElements(); i++) {
                (*(summap(2)))[i] += weight * temp_summap[i];
            }

            //log(Comm.MyPID(),"Creating M");
            createM(PixMap, Map, PQ, PU, weight, dm->NSTOKES, invM);
            delete PQ;
            delete PU;
            delete Graph;
            log(Comm.MyPID(),"GlobalAssemble");
            invM.GlobalAssemble();
            log(Comm.MyPID(),"GlobalAssemble DONE");
        }
    }
//end LOOP
delete y;
//
//log(Comm.MyPID(),"HITMAP");
//Epetra_Vector * hitmap;
//hitmap = new Epetra_Vector(PixMap);
//createHitmap(PixMap, *hitmap, invM);
//WriteH5Vec(*hitmap, "hitmap");
//delete hitmap;
WriteH5Vec(*summap(0), "summap");

//Epetra_Vector binmap(PixMap);
//log(Comm.MyPID(),"Computing RCOND and Inverting");
//Epetra_Vector rcond(PixMap);
//invertM(PixMap, invM, rcond);
//
//cout << time.ElapsedTime() << endl;
//
//log(Comm.MyPID(),"BINMAP");
//invM.Apply(summap, binmap);
//
//log(Comm.MyPID(),"Writing MAPS");
//
//WriteH5Vec(binmap, "binmap");
//WriteH5Vec(rcond, "rcondmap");
//cout << time.ElapsedTime() << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
return(0);

};
