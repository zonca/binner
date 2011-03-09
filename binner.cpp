#include <string>
#include <math.h>

#include <boost/scoped_ptr.hpp>

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

int WriteH5Vec(Epetra_Vector& vec, string filename, int NPIX) {
    int MyPID = vec.Comm().MyPID();
    H5std_string  FILE_NAME( str( format("%s_%03d.h5") % filename % MyPID ) );
    H5File file(FILE_NAME, H5F_ACC_TRUNC );
    hsize_t dimsf[1];
    dimsf[0] = vec.Map().NumMyPoints();
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

int SAMPLES_PER_PROC = 4000000;

Epetra_Time time(Comm);
H5PlanckDataManager* dm;

string parameterFilename = "notimplemented.dat";
readParameterFile(parameterFilename, dm);

log(Comm.MyPID(), format("Number of elements: %.10d") % dm->getDatasetLength());
//Epetra_BlockMap Map(dm->getDatasetLength(), 1, 0, Comm);
Epetra_Map Map(-1, SAMPLES_PER_PROC, 0, Comm);

int NumMyElements = Map.NumMyElements();
cout << Comm.MyPID() << " " << Map.MinMyGID() << " " << NumMyElements << endl;

Epetra_BlockMap PixMap(dm->getNPIX(),dm->NSTOKES,0,Comm);

// declaring matrices
log(Comm.MyPID(),"POINTING MATRIX");
Epetra_Vector summap(PixMap);
Epetra_Vector temp_summap(PixMap);
Epetra_VbrMatrix * P;
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

//LOOP
for (int offset=0; offset<dm->getDatasetLength(); offset=offset+Map.NumGlobalElements()) {

    log(Comm.MyPID(),format("Offset: %d") % offset);
    P = new Epetra_VbrMatrix(Copy, Map, 1);
    createP(Map, PixMap, dm, offset, P);

    log(Comm.MyPID(),"READ DATA");
    data = new double[Map.NumMyElements()];
    dm->getData(Map.MinMyGID() + offset,Map.NumMyElements(),data);
    y = new Epetra_Vector(Copy,Map,data);
    delete[] data;

    cout << time.ElapsedTime() << endl;

    log(Comm.MyPID(),"SUM MAP");
    P->Multiply1(true,*y,temp_summap); //SUMMAP = Pt y
    for (int i=0; i<PixMap.NumMyPoints(); i++) {
        summap[i] += temp_summap[i];
    }
    delete y;

    log(Comm.MyPID(),"Creating M");
    createM(PixMap, Map, P, dm->NSTOKES, invM);
    log(Comm.MyPID(),"GlobalAssemble");
    invM.GlobalAssemble();
    log(Comm.MyPID(),"GlobalAssemble DONE");
    delete P;
}
//end LOOP
//
log(Comm.MyPID(),"HITMAP");
Epetra_Vector * hitmap;
hitmap = new Epetra_Vector(PixMap);
createHitmap(PixMap, *hitmap, invM);
WriteH5Vec(*hitmap, "hitmap", dm->getNPIX());
delete hitmap;
WriteH5Vec(summap, "summap", dm->getNPIX());

Epetra_Vector binmap(PixMap);
log(Comm.MyPID(),"Computing RCOND and Inverting");
Epetra_Vector rcond(PixMap);
invertM(PixMap, invM, rcond);

cout << time.ElapsedTime() << endl;

log(Comm.MyPID(),"BINMAP");
invM.Apply(summap, binmap);

log(Comm.MyPID(),"Writing MAPS");

WriteH5Vec(binmap, "binmap", dm->getNPIX());
WriteH5Vec(rcond, "rcondmap", dm->getNPIX());
cout << time.ElapsedTime() << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
return(0);

};
