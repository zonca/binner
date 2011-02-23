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
#include <Epetra_SerialDenseSolver.h>



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
bool DEBUG = false;
int NSTOKES = 3;

string dataPath, pointingPath;

//dataPath = "/home/zonca/p/testdata/lfi_ops_dx4";
dataPath = "/global/homes/z/zonca/planck/data/mission/lfi_ops_dx4";


if (DEBUG) {
//        pointingPath = "/home/zonca/p/testdata/dx4_1_nest";
        pointingPath = "/global/homes/z/zonca/p/pointing/dx4_1_nest";

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
    NumElements = 10000;
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

int BlockIndices[1];
Epetra_SerialDenseMatrix *Prow;
int RowDim, NumBlockEntries;
int *BlockIndicesOut;
int err;
Epetra_SerialDenseMatrix * Mpp;
Epetra_SerialDenseMatrix * Zero;

log(Comm.MyPID(),"Initializing M");

for( int i=0 ; i<PixMap.NumMyElements(); ++i ) { //loop on local pixel
    BlockIndices[0] = PixMyGlobalElements[i];
    Zero = new Epetra_SerialDenseMatrix(NSTOKES, NSTOKES);
    invM.BeginInsertGlobalValues(BlockIndices[0], 1, BlockIndices);
    err = invM.SubmitBlockEntry(Zero->A(), Zero->LDA(), NSTOKES, NSTOKES);
            if (err != 0) {
                cout << "PID:" << Comm.MyPID() << "Error in inserting init zero values in M, error code:" << err << endl;
                }
    err = invM.EndSubmitEntries();
    delete Zero;
    }

for( int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local pointing

    P.BeginExtractMyBlockRowView(i, RowDim, NumBlockEntries, BlockIndicesOut);
    P.ExtractEntryView(Prow);

    BlockIndices[0] = PixMyGlobalElements[BlockIndicesOut[0]];

    Mpp = new Epetra_SerialDenseMatrix(NSTOKES, NSTOKES);

    err = Mpp->Multiply('T','N', 1., *Prow, *Prow, 0.);
            if (err != 0) {
                cout << "Error in computing Mpp, error code:" << err << endl;
                }

    invM.BeginSumIntoGlobalValues(BlockIndices[0], 1, BlockIndices);

    err = invM.SubmitBlockEntry(Mpp->A(), Mpp->LDA(), NSTOKES, NSTOKES); //FIXME check order
            if (err != 0) {
                cout << "PID:" << Comm.MyPID() << "Error in inserting values in M, error code:" << err << endl;
                }

    err = invM.EndSubmitEntries();
            if (err != 0) {
                cout << "PID:" << Comm.MyPID() << " LocalRow[i]:" << i << " Error in ending submit entries in M, error code:" << err << endl;
                }
    delete Mpp;

}

log(Comm.MyPID(),"GlobalAssemble");
invM.GlobalAssemble();
log(Comm.MyPID(),"GlobalAssemble DONE");
cout << time.ElapsedTime() << endl;

log(Comm.MyPID(),"Computing RCOND and Inverting");
Epetra_SerialDenseSolver * SSolver;

double rcond_blockM;
Epetra_SerialDenseMatrix * blockM;

Epetra_Vector binmap(PixMap);
Epetra_Vector rcond(PixMap);
int RCondIndices[1];
double RCondValues[1];

for( int i=0 ; i<PixMap.NumMyElements(); ++i ) { //loop on local pointing

    invM.BeginExtractMyBlockRowView(i, RowDim, NumBlockEntries, BlockIndicesOut);
    invM.ExtractEntryView(blockM);

    SSolver = new Epetra_SerialDenseSolver();
    SSolver->SetMatrix(*blockM);
    //cout << "PID:" << Comm.MyPID() << " localPIX:" << i << " globalPIX:" << PixMyGlobalElements[i] << endl;
    if ((*blockM)(0,0) > 0) {
        rcond_blockM = PixMyGlobalElements[i];
        //err = SSolver->ReciprocalConditionEstimate(rcond_blockM);
        //if (err != 0) {
        //    cout << "PID:" << Comm.MyPID() << " LocalRow[i]:" << i << " cannot compute RCOND, error code:" << err << " estimated:"<< rcond_blockM << endl;
        //}
        //if (rcond_blockM > 1e-5) {
        //    err = SSolver->Invert();
        //    if (err != 0) {
        //        cout << "PID:" << Comm.MyPID() << " LocalRow[i]:" << i << " cannot invert matrix, error code:" << err << endl;
        //    }
        //} else {
        //    for (int r=0; r<3; ++r) {
        //        for (int c=0; c<3; ++c) {
        //            (*blockM)(r,c) = 0.;
        //        }
        //    }
        //}
    } else {
        rcond_blockM = -1;
    }
    RCondValues[0] = rcond_blockM;
    RCondIndices[0] = PixMyGlobalElements[i];
    rcond.ReplaceGlobalValues(1, 0, RCondValues, RCondIndices);
}

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
