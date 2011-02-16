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
dataPath = "/global/homes/z/zonca/planck/data/mission/lfi_ops_dx4";

if (DEBUG) {
    pointingPath = "/global/homes/z/zonca/p/pointing/dx4_1_nest";
    NSIDE = 1;
} else {
    pointingPath = "/global/homes/z/zonca/p/pointing/dx4_1024_nest";
    NSIDE = 1024;
}
//const list<string> channels = list_of( "ch1q" )( "ch1u" );
const list<string> channels = list_of( "LFI28M" );
//const list<string> channels = list_of("LFI27M" )( "LFI27S" )( "LFI28M" )( "LFI28S" );
dm = new PlanckDataManager(92, 93, channels, dataPath, pointingPath);

//MAPS
// Construct a Map with NumElements and index base of 0
//
int NumElements;
if (DEBUG) {
    NumElements = 80;
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
int Values[NSTOKES * NSTOKES];
Epetra_SerialDenseMatrix Mpp(NSTOKES, NSTOKES);

for( int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local pointing

    P.BeginExtractMyBlockRowView(i, RowDim, NumBlockEntries, BlockIndicesOut);
    P.ExtractEntryView(Prow);
    //cout << *Prow << endl;
    err = Mpp.Multiply('T','N', 1., *Prow, *Prow, 1.);
            if (err != 0) {
                cout << "Error in computing Mpp, error code:" << err << endl;
                }
    //cout << Mpp << endl;
    BlockIndices[0] = 0;
    //invM.BeginSumIntoGlobalValues(BlockIndicesOut[0], 1, BlockIndicesOut);
    invM.BeginSumIntoGlobalValues(0, 1, BlockIndices);

    err = invM.SubmitBlockEntry(Mpp.A(), Mpp.LDA(), NSTOKES, NSTOKES); //FIXME check order
            if (err != 0) {
                cout << "Error in inserting values in M, error code:" << err << endl;
                }

    err = invM.EndSubmitEntries();
            if (err != 0) {
                cout << "Error in ending submit entries in M, error code:" << err << endl;
                }

}
invM.GlobalAssemble();

cout << invM << endl;

//log(Comm.MyPID(),"M-M");
//int err = EpetraExt::MatrixMatrix::Multiply(P, true, P, false, invM);
//log(Comm.MyPID(),"M-M END");
//cout << time.ElapsedTime() << endl;
//if (err != 0) {
//    cout << "err "<<err<<" from MatrixMatrix::Multiply"<<endl;
//    return(err);
//}
//
//log(Comm.MyPID(),"invM");
//Epetra_Vector diag(PixMap);
//
//invM.ExtractDiagonalCopy(diag);
//
//log(Comm.MyPID(),"Exporting HITMAP");
//
MapWriter mapwriter(PixMap, Comm, NPIX);
//mapwriter.write(diag, "hitmap.fits");
//log(Comm.MyPID(),"inverting");
//
//invertM(PixMap, invM);
//Epetra_Vector binmap(PixMap);
//invM.Multiply(false, summap, binmap); //B = inv(M) * Pt
//
////Write map to file on proc 0
//mapwriter.write(binmap, "binmap.fits");
mapwriter.write(summap, "summap.fits");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
return(0);

};
