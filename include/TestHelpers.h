#include <math.h>
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

using namespace std;

int createFakeM(const Epetra_Map Map, Epetra_CrsMatrix& M) {
    double * Values = new double[2];
    int * Indices = new int[2];
    double * ValuesSecRow = new double[2];
    int * IndicesSecRow = new int[2];
    unsigned int secrow;
    for(unsigned int row=0; row<Map.NumGlobalElements(); row=row+2) {

        Values[0] = (row + 1) * 2.5;
        Values[1] = (row + 1) * 1.5;
        Indices[0] = row;
        Indices[1] = row + 1;
        M.InsertGlobalValues(row, 2, Values, Indices);

        secrow = row + 1;
        ValuesSecRow[0] = (row + 1) * 1.5;
        ValuesSecRow[1] = 3.5;
        IndicesSecRow[0] = secrow - 1;
        IndicesSecRow[1] = secrow;
        M.InsertGlobalValues(row + 1, 2, ValuesSecRow, IndicesSecRow);
    }
    M.FillComplete();
}

int createFakeP(const Epetra_Map& Map, const Epetra_BlockMap& PixMap, Epetra_VbrMatrix& P) {

    Epetra_IntSerialDenseVector pointing(Map.NumMyElements());
    pointing[0]=3;
    pointing[1]=1;
    pointing[2]=4;
    pointing[3]=4;
    pointing[4]=4;
    pointing[5]=1;
    pointing[6]=1;
    pointing[7]=1;
    pointing[8]=3;
    pointing[9]=3;

    double ang;

    if (Map.Comm().MyPID() == 0) {
        ang = 10. / 180.* 3.14;
    } else {
        ang = 100. / 180.* 3.14;
    }

    int err;
    int Indices[1];
    double* Values; Values = new double[3];
    int * MyGlobalElements = Map.MyGlobalElements();

    for( int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local rows

            int GlobalNode = MyGlobalElements[i];
            Indices[0] = int(pointing[i]);

            err = P.BeginInsertGlobalValues(GlobalNode, 1, Indices);

            Values[0] = 1.;
            Values[1] = cos(2 * (ang + 4*i));
            Values[2] = sin(2 * (ang + 4*i));

            err = P.SubmitBlockEntry(Values, 1, 1, 3);
            if (err != 0) {
                cout << "Error in submitting entries for P, error code:" << err << endl;
                }
            err = P.EndSubmitEntries();
    }
    P.FillComplete(PixMap, Map);

}
