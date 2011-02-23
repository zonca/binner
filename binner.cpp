#include <string>
#include <math.h>
#include "PlanckDataManager.h"
#include "MapWriter.h"
#include "CreateMatrices.h"
#include "ReadParameterFile.h"
#include "Utils.h"

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
#include "Epetra_DataAccess.h"
#include "Epetra_Time.h"
#include "Epetra_FEVbrMatrix.h"

using namespace std;

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif  

Epetra_Time time(Comm);
PlanckDataManager* dm;

string parameterFilename = "notimplemented.dat";
readParameterFile(parameterFilename, dm);

log(Comm.MyPID(), format("Number of elements: %d") % dm->getDatasetLength());
Epetra_Map Map(dm->getDatasetLength(), 0, Comm);

Epetra_BlockMap PixMap(dm->getNPIX(),dm->NSTOKES,0,Comm);

log(Comm.MyPID(),"POINTING MATRIX");
Epetra_VbrMatrix P(Copy,Map,1);
createP(Map, PixMap, dm, P);
cout << time.ElapsedTime() << endl;

log(Comm.MyPID(),"READ DATA");
double* data; data = new double[Map.NumMyElements()];
dm->getData("data", Map.MinMyGID() ,Map.NumMyElements(),data);
Epetra_Vector y(Copy,Map,data);
delete data;
cout << time.ElapsedTime() << endl;

log(Comm.MyPID(),"SUM MAP");
Epetra_Vector summap(PixMap);
P.Multiply1(true,y,summap); //SUMMAP = Pt y

Epetra_FEVbrMatrix invM(Copy, PixMap, 1);

log(Comm.MyPID(),"Initializing M");
initM(PixMap, dm->NSTOKES, invM);

log(Comm.MyPID(),"Creating M");
createM(PixMap, Map, P, dm->NSTOKES, invM);

log(Comm.MyPID(),"GlobalAssemble");
invM.GlobalAssemble();
log(Comm.MyPID(),"GlobalAssemble DONE");
cout << time.ElapsedTime() << endl;

Epetra_Vector binmap(PixMap);
log(Comm.MyPID(),"Computing RCOND and Inverting");
Epetra_Vector rcond(PixMap);
invertM(PixMap, invM, rcond);

cout << time.ElapsedTime() << endl;

log(Comm.MyPID(),"BINMAP");
invM.Apply(summap, binmap);

log(Comm.MyPID(),"Writing MAPS");

MapWriter mapwriter(PixMap, Comm, dm->getNPIX());
//mapwriter.write(binmap, "binmap.fits");
mapwriter.write(rcond, "rcondmap.fits");
mapwriter.write(summap, "summap.fits");
cout << time.ElapsedTime() << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
return(0);

};
