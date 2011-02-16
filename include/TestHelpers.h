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
