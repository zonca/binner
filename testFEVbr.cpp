#include <list>
#include <string>
#include <math.h>

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_FEVbrMatrix.h"

//#include "AztecOO.h"

#include <EpetraExt_MatrixMatrix.h>

using namespace std;

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif  

int NSIDE = 1;
int NPIX ;
int NSTOKES = 3;

int NumElements;
NumElements = 11;

NPIX = 12. * NSIDE * NSIDE +1; //total pixel size, each pixel is an element which contains 3 floats which are IQU
Epetra_BlockMap PixMap(NPIX,NSTOKES,0,Comm);

int * PixMyGlobalElements = PixMap.MyGlobalElements();

cout <<  PixMap << endl;
Epetra_FEVbrMatrix invM(Copy, PixMap, 1);

int BlockIndices[1];
BlockIndices[0] = 2;
Epetra_SerialDenseMatrix *Prow;
int RowDim, NumBlockEntries;
int err;
Epetra_SerialDenseMatrix Mpp(NSTOKES, NSTOKES);
Mpp[0][0] = 1.;
cout << Mpp << endl;

int debugPID = 1;

int NumHits = 2*Comm.MyPID() + 5;
for( int i=0 ; i<NumHits; ++i ) { //loop on local pointing

    invM.BeginSumIntoGlobalValues(BlockIndices[0], 1, BlockIndices);

    err = invM.SubmitBlockEntry(Mpp.A(), Mpp.LDA(), NSTOKES, NSTOKES); //FIXME check order
            if (err != 0) {
                cout << "PID:" << Comm.MyPID() << "Error in inserting values in M, error code:" << err << endl;
                }

    err = invM.EndSubmitEntries();
            if (err != 0) {
                cout << "PID:" << Comm.MyPID() << " LocalRow[i]:" << i << " Error in ending submit entries in M, error code:" << err << endl;
                }

}
invM.GlobalAssemble();

cout << invM << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

return(0);

};
