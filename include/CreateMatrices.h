#ifndef CREATEMATRICES_H
#define CREATEMATRICES_H

#include <string>
#include <vector>

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_DataAccess.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"

#include "PlanckDataManager.h"

using namespace std;

int createF(const Epetra_Map& Map, const Epetra_Map& BaselinesMap, int BaselineLength, Epetra_CrsMatrix& F) {

    double one = 1.0;
    int baselinenum;

    int * MyGlobalElements = Map.MyGlobalElements();
    for(unsigned int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local rows
        
        baselinenum = MyGlobalElements[i] / BaselineLength;
        F.InsertGlobalValues(MyGlobalElements[i], 1, &one, &baselinenum);
    }

    F.FillComplete(BaselinesMap, Map);
    return 0;
}

int createP(const Epetra_Map& Map, const Epetra_BlockMap& PixMap, PlanckDataManager* dm, Epetra_VbrMatrix& P) {

    int * MyGlobalElements = Map.MyGlobalElements();
    int NumMyElements = Map.NumMyElements();
    int err;

    double* pointing; pointing = new double[NumMyElements];
    double* qw; qw = new double[NumMyElements];
    double* uw; uw = new double[NumMyElements];
    dm->getData("pointing", Map.MinMyGID(), NumMyElements, pointing);
    dm->getData("qw", Map.MinMyGID(), NumMyElements, qw);
    dm->getData("uw", Map.MinMyGID(), NumMyElements, uw);

    //double Values[3];
    double* Values; Values = new double[3];

    int Indices[1];
    int MyPID = Map.Comm().MyPID();

    //for( int i=0 ; i<1; ++i ) { //loop on local rows
    for( int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local rows
            int GlobalNode = MyGlobalElements[i];
            Indices[0] = int(pointing[i]);
            //cout << "0: " <<Indices[0] << " 1: " << Indices[1] << " 2: " << Indices[2] << endl;
            err = P.BeginInsertGlobalValues(GlobalNode, 1, Indices);

            if (err != 0) {
                cout << "Error in inserting values in P, error code:" << err << endl;
                }

        
            Values[0] = 1.;
            Values[1] = qw[i];
            Values[2] = uw[i];
            err = P.SubmitBlockEntry(Values, 1, 1, 3);
            if (err != 0) {
                cout << "Error in submitting entries for P, error code:" << err << endl;
                }
            err = P.EndSubmitEntries();
            if (err != 0) {
                cout << "Error in ending submit entries for P, error code:" << err << endl;
            cout << "INDICES 0: " <<Indices[0] << endl;
            cout << "GlobalNode:" << GlobalNode<< endl;
            cout << "VALUES 0: " <<Values[0] << " 1: " << Values[1] << " 2: " << Values[2] << endl;
                }
    }

    P.FillComplete(PixMap, Map);
    //P.FillComplete();
    return 0;
}

int invertM(const Epetra_Map& PixMap, Epetra_CrsMatrix& invM) {

    double * Values = new double[2];
    double * ValuesSecRow = new double[2];
    int * Indices = new int[2];
    int NumMyElements = PixMap.NumMyElements();
    int * MyGlobalElements = new int [NumMyElements];
    PixMap.MyGlobalElements( MyGlobalElements );

    int error;
    double det;
    for(unsigned int i=0; i<PixMap.NumMyElements(); i=i+2) {
        if (invM.NumMyEntries(i) > 0 ){
            Indices[0] = MyGlobalElements[i];
            Indices[1] = MyGlobalElements[i+1];

                det = invM[i][0] * invM[i+1][1] - invM[i][1] * invM[i+1][0];
                //cout << det << endl;
                Values[0] = invM[i+1][1] / det;
                Values[1] = -1. * invM[i][1] / det;
                ValuesSecRow[0] = -1 * invM[i+1][0] / det;
                ValuesSecRow[1] = invM[i][0] / det;

                error = invM.ReplaceGlobalValues(MyGlobalElements[i], 2, Values, Indices);

                if (error != 0){
                    cout << "error:" << error << "-i" << i << "-Elem[" << MyGlobalElements[i] << ",0]=" <<  invM[i][0] << "Values:" << Values[0] << endl;
                }

                error = invM.ReplaceGlobalValues(MyGlobalElements[i+1], 2, ValuesSecRow, Indices);
        }
        
    }
    return 0;
}
#endif
