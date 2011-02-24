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
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include <Epetra_SerialDenseSolver.h>

#include "PlanckDataManager.h"

using namespace std;

void createP(const Epetra_Map& Map, const Epetra_BlockMap& PixMap, PlanckDataManager* dm, Epetra_VbrMatrix& P) {

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
}

void initM(const Epetra_BlockMap& PixMap, int NSTOKES, Epetra_FEVbrMatrix& invM) {

    int BlockIndices[1];
    int * PixMyGlobalElements = PixMap.MyGlobalElements();
    int err;
    Epetra_SerialDenseMatrix * Zero;

    for( int i=0 ; i<PixMap.NumMyElements(); ++i ) { //loop on local pixel
        BlockIndices[0] = PixMyGlobalElements[i];
        Zero = new Epetra_SerialDenseMatrix(NSTOKES, NSTOKES);
        invM.BeginInsertGlobalValues(BlockIndices[0], 1, BlockIndices);
        err = invM.SubmitBlockEntry(Zero->A(), Zero->LDA(), NSTOKES, NSTOKES);
                if (err != 0) {
                    cout << "PID:" << PixMap.Comm().MyPID() << "Error in inserting init zero values in M, error code:" << err << endl;
                    }
        err = invM.EndSubmitEntries();
        delete Zero;
    }
}

void createM(const Epetra_BlockMap& PixMap, const Epetra_BlockMap& Map, const Epetra_VbrMatrix& P, int NSTOKES, Epetra_FEVbrMatrix& invM) {

    Epetra_SerialDenseMatrix * Mpp;
    int err;

    int * PixMyGlobalElements = PixMap.MyGlobalElements();
    Epetra_SerialDenseMatrix *Prow;
    int RowDim, NumBlockEntries;
    int *BlockIndicesOut;
    for( int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local pointing

        P.BeginExtractMyBlockRowView(i, RowDim, NumBlockEntries, BlockIndicesOut);
        P.ExtractEntryView(Prow);

        Mpp = new Epetra_SerialDenseMatrix(NSTOKES, NSTOKES);

        err = Mpp->Multiply('T','N', 1., *Prow, *Prow, 0.);
            if (err != 0) {
                cout << "Error in computing Mpp, error code:" << err << endl;
                }

        invM.BeginSumIntoGlobalValues(BlockIndicesOut[0], 1, BlockIndicesOut);

        err = invM.SubmitBlockEntry(Mpp->A(), Mpp->LDA(), NSTOKES, NSTOKES); //FIXME check order
                if (err != 0) {
                    cout << "PID:" << PixMap.Comm().MyPID() << "Error in inserting values in M, error code:" << err << endl;
                    }

        err = invM.EndSubmitEntries();
                if (err != 0) {
                    cout << "PID:" << PixMap.Comm().MyPID() << " LocalRow[i]:" << i << " Error in ending submit entries in M, error code:" << err << endl;
                    }
        delete Mpp;

    }

}

void createHitmap(const Epetra_BlockMap& PixMap, Epetra_Vector& hitmap, Epetra_FEVbrMatrix& invM) {
    int * PixMyGlobalElements = PixMap.MyGlobalElements();
    int RCondIndices[1], err;
    double RCondValues[1];
    int RowDim, NumBlockEntries;
    int *BlockIndicesOut;
    Epetra_SerialDenseMatrix * blockM;

    for( int i=0 ; i<PixMap.NumMyElements(); ++i ) { //loop on local pointing
        invM.BeginExtractMyBlockRowView(i, RowDim, NumBlockEntries, BlockIndicesOut);
        invM.ExtractEntryView(blockM);

        RCondValues[0] = (*blockM)(0,0);
        RCondIndices[0] = PixMyGlobalElements[i];
        hitmap.ReplaceGlobalValues(1, 0, RCondValues, RCondIndices);
    }
}

void invertM(const Epetra_BlockMap& PixMap, Epetra_FEVbrMatrix& invM, Epetra_Vector& rcond) {

    Epetra_SerialDenseSolver * SSolver;
    int * PixMyGlobalElements = PixMap.MyGlobalElements();

    double rcond_blockM;
    Epetra_SerialDenseMatrix * blockM;

    int RCondIndices[1], err;
    double RCondValues[1];

    int RowDim, NumBlockEntries;
    int *BlockIndicesOut;
    for( int i=0 ; i<PixMap.NumMyElements(); ++i ) { //loop on local pointing

        invM.BeginExtractMyBlockRowView(i, RowDim, NumBlockEntries, BlockIndicesOut);
        invM.ExtractEntryView(blockM);

        SSolver = new Epetra_SerialDenseSolver();
        SSolver->SetMatrix(*blockM);
        //cout << "PID:" << Comm.MyPID() << " localPIX:" << i << " globalPIX:" << PixMyGlobalElements[i] << endl;
        if ((*blockM)(0,0) > 0) {
            //rcond_blockM = PixMyGlobalElements[i];
            err = SSolver->ReciprocalConditionEstimate(rcond_blockM);
            if (err != 0) {
                cout << "PID:" << PixMap.Comm().MyPID() << " LocalRow[i]:" << i << " cannot compute RCOND, error code:" << err << " estimated:"<< rcond_blockM << endl;
            }
            if (rcond_blockM > 1e-5) {
                err = SSolver->Invert();
                if (err != 0) {
                    cout << "PID:" << PixMap.Comm().MyPID() << " LocalRow[i]:" << i << " cannot invert matrix, error code:" << err << endl;
                }
            } else {
                for (int r=0; r<blockM->M(); ++r) {
                    for (int c=0; c<blockM->N(); ++c) {
                        (*blockM)(r,c) = 0.;
                    }
                }
            }
        } else {
            rcond_blockM = -1;
        }
        RCondValues[0] = rcond_blockM;
        RCondIndices[0] = PixMyGlobalElements[i];
        rcond.ReplaceGlobalValues(1, 0, RCondValues, RCondIndices);
    }
}
#endif
