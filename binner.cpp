#include <list>
#include <string>
#include <math.h>
#include "PlanckDataManager.h"
#include "MapWriter.h"
#include "CreateMatrices.h"



#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/assign/list_of.hpp> // for 'list_of()'

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
#include "Epetra_CrsMatrix.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Time.h"
#include "Epetra_FEVbrMatrix.h"



//#include "AztecOO.h"

#include <EpetraExt_MatrixMatrix.h>

extern "C" {
    #include "read_fits_lib.h"
}

using boost::format;
using namespace std;
using namespace boost::assign; // bring 'list_of()' into scope

int log(int MyPID, string message) {

    if( MyPID == 0 ) {
        cout << "******** " << message << endl << endl << flush;
    }
    return true;
} 

int log(int MyPID, format message) {

    log(MyPID, message.str());
    return true;

} 

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif  

Epetra_Time time(Comm);

int NSIDE;
int NPIX ;
PlanckDataManager* dm;
bool DEBUG = true;
int NSTOKES = 3;

string dataPath, pointingPath;

dataPath = "/home/zonca/p/testdata/lfi_ops_dx4";
//dataPath = "/global/homes/z/zonca/planck/data/mission/lfi_ops_dx4";

if (DEBUG) {
        pointingPath = "/home/zonca/p/testdata/dx4_1_nest";
//        pointingPath = "/global/homes/z/zonca/p/pointing/dx4_1_nest";

    NSIDE = 1;
} else {
    pointingPath = "/global/homes/z/zonca/p/pointing/dx4_1024_nest";
    NSIDE = 1024;
}
//const list<string> channels = list_of( "LFI28M" );
const list<string> channels = list_of("LFI27M" )( "LFI27S" )( "LFI28M" )( "LFI28S" );
dm = new PlanckDataManager(92, 93, channels, dataPath, pointingPath);

//MAPS
// Construct a Map with NumElements and index base of 0
//
int NumElements;
if (DEBUG) {
    NumElements = 1000;
} else {
    NumElements = dm->getDatasetLength();
}

log(Comm.MyPID(), format("Number of elements: %d") % NumElements);
Epetra_Map Map(NumElements, 0, Comm);

NPIX = 12. * NSIDE * NSIDE +1; //total pixel size, each pixel is an element which contains 3 floats which are IQU
Epetra_BlockMap PixMap(NPIX,NSTOKES,0,Comm);

int * MyGlobalElements = Map.MyGlobalElements();
int * PixMyGlobalElements = PixMap.MyGlobalElements();

log(Comm.MyPID(),"POINTING MATRIX");
Epetra_VbrMatrix P(Copy,Map,1);
createP(Map, PixMap, dm, P);
//cout << P;
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

//cout << P << endl;
//log(Comm.MyPID(),"/////////////// Creating M pix x pix");
Epetra_FEVbrMatrix invM(Copy, PixMap, 1);
//Epetra_FEVbrMatrix invM(Copy, PixMap, PixMap, 1);

int err;

log(Comm.MyPID(),"Initializing M");

initM(PixMap, NSTOKES, invM);

log(Comm.MyPID(),"Creating M");
createM(PixMap, Map, P, NSTOKES, invM);

log(Comm.MyPID(),"GlobalAssemble");
invM.GlobalAssemble();
log(Comm.MyPID(),"GlobalAssemble DONE");
cout << time.ElapsedTime() << endl;

Epetra_Vector binmap(PixMap);
log(Comm.MyPID(),"Computing RCOND and Inverting");
Epetra_Vector rcond(PixMap);
invertM(PixMap, invM, rcond);

cout << time.ElapsedTime() << endl;

Comm.Barrier();

log(Comm.MyPID(),"BINMAP");
invM.Apply(summap, binmap);

log(Comm.MyPID(),"Writing MAPS");

MapWriter mapwriter(PixMap, Comm, NPIX);
////Write map to file on proc 0
//mapwriter.write(binmap, "binmap.fits");
mapwriter.write(rcond, "rcondmap.fits");
mapwriter.write(summap, "summap.fits");
cout << time.ElapsedTime() << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
return(0);

};
