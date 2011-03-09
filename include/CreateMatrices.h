#ifndef CREATEMATRICES_H
#define CREATEMATRICES_H

#include <string>

#include <boost/scoped_array.hpp>
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
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include <Epetra_SerialDenseSolver.h>

#include "H5PlanckDataManager.h"
#include "Utils.h"

using namespace std;

void createP(const Epetra_BlockMap& Map, const Epetra_BlockMap& PixMap, H5PlanckDataManager* dm, int offset,Epetra_VbrMatrix* P) {

    int * MyGlobalElements = Map.MyGlobalElements();

    int MyPID = Map.Comm().MyPID();
    int NumMyElements = Map.NumMyElements();
    int err;

    boost::scoped_array<pointing_t> pointing(new pointing_t[NumMyElements]);
    log(MyPID, "Reading pointing");

    cout << MyPID << " " << Map.MinMyGID() + offset << " " << NumMyElements << endl;
    dm->getPointing(Map.MinMyGID() + offset, NumMyElements, pointing.get());

    boost::scoped_array<double> Values(new double[dm->NSTOKES]);

    int Indices[1];

    log(MyPID, "Assembling P");
    for( int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local rows
            int GlobalNode = MyGlobalElements[i];
            Indices[0] = pointing[i].pix;
            //cout << "0: " <<Indices[0] << " 1: " << Indices[1] << " 2: " << Indices[2] << endl;
            err = P->BeginInsertGlobalValues(GlobalNode, 1, Indices);

            if (err != 0) {
                cout << "Error in inserting values in P, error code:" << err << endl;
            }

            Values[0] = 1.;
            Values[1] = pointing[i].qw;
            Values[2] = pointing[i].uw;

            err = P->SubmitBlockEntry(Values.get(), 1, 1, 3);
            if (err != 0) {
                cout << "Error in submitting entries for P, error code:" << err << endl;
            }
            err = P->EndSubmitEntries();
            if (err != 0) {
                cout << "Error in ending submit entries for P, error code:" << err << endl;
                cout << "INDICES 0: " <<Indices[0] << endl;
                cout << "GlobalNode:" << GlobalNode<< endl;
                cout << "VALUES 0: " <<Values[0] << " 1: " << Values[1] << " 2: " << Values[2] << endl;
            }
    }

    log(MyPID, "FillComplete P");
    P->FillComplete(PixMap, Map);
    log(MyPID, "FillComplete P done");
}

void initM(const Epetra_BlockMap& PixMap, int NSTOKES, Epetra_FEVbrMatrix& invM) {

    int BlockIndices[1];
    int * PixMyGlobalElements = PixMap.MyGlobalElements();
    int err;
    boost::scoped_ptr<Epetra_SerialDenseMatrix> Zero (new Epetra_SerialDenseMatrix(NSTOKES, NSTOKES));

    for( int i=0 ; i<PixMap.NumMyElements(); ++i ) { //loop on local pixel
        BlockIndices[0] = PixMyGlobalElements[i];
        invM.BeginInsertGlobalValues(BlockIndices[0], 1, BlockIndices);
        err = invM.SubmitBlockEntry(Zero->A(), Zero->LDA(), NSTOKES, NSTOKES);
                if (err != 0) {
                    cout << "PID:" << PixMap.Comm().MyPID() << "Error in inserting init zero values in M, error code:" << err << endl;
                    }
        err = invM.EndSubmitEntries();
    }
}

void createM(const Epetra_BlockMap& PixMap, const Epetra_BlockMap& Map, const Epetra_VbrMatrix* P, int NSTOKES, Epetra_FEVbrMatrix& invM) {

    int err;

    int * PixMyGlobalElements = PixMap.MyGlobalElements();
    Epetra_SerialDenseMatrix *Prow;
    int RowDim, NumBlockEntries;
    int * LocalPix;

    boost::scoped_array<int> GlobalPix (new int[1]);
    boost::scoped_ptr<Epetra_SerialDenseMatrix> Mpp (new Epetra_SerialDenseMatrix(NSTOKES, NSTOKES));

    for( int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local pointing

        P->BeginExtractMyBlockRowView(i, RowDim, NumBlockEntries, LocalPix);
        P->ExtractEntryView(Prow);
        GlobalPix[0] = P->GCID(LocalPix[0]);
        //cout << Map.Comm().MyPID() << ": i:" << i << " Loc:" << LocalPix[0] << " Glob:" << GlobalPix[0] << " " << endl;

        err = Mpp->Multiply('T','N', 1., *Prow, *Prow, 0.);
            if (err != 0) {
                cout << "Error in computing Mpp, error code:" << err << endl;
                }

        invM.BeginSumIntoGlobalValues(GlobalPix[0], 1, GlobalPix.get());
        //invM.BeginSumIntoLocalValues(GlobalPix[0], 1, GlobalPix);

        err = invM.SubmitBlockEntry(Mpp->A(), Mpp->LDA(), NSTOKES, NSTOKES); //FIXME check order
                if (err != 0) {
                    cout << "PID:" << PixMap.Comm().MyPID() << "Error in inserting values in M, error code:" << err << endl;
                    }

        err = invM.EndSubmitEntries();
                if (err != 0) {
                    cout << "PID:" << PixMap.Comm().MyPID() << " LocalRow[i]:" << i << " Error in ending submit entries in M, error code:" << err << endl;
                    }

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
        err = hitmap.ReplaceGlobalValues(1, 0, RCondValues, RCondIndices);
                if (err != 0) {
                    cout << "PID:" << PixMap.Comm().MyPID() << " LocalRow[i]:" << i << " Error in creating HITMAP:" << err << endl;
                    }
    }
}

void invertM(const Epetra_BlockMap& PixMap, Epetra_FEVbrMatrix& invM, Epetra_Vector& rcond) {

    int * PixMyGlobalElements = PixMap.MyGlobalElements();
    boost::scoped_ptr<Epetra_SerialDenseSolver> SSolver (new Epetra_SerialDenseSolver());

    double rcond_blockM;
    Epetra_SerialDenseMatrix * blockM;

    int RCondIndices[1], err;
    double RCondValues[1];

    int RowDim, NumBlockEntries;
    int *BlockIndicesOut;
    for( int i=0 ; i<PixMap.NumMyElements(); ++i ) { //loop on local pointing

        invM.BeginExtractMyBlockRowView(i, RowDim, NumBlockEntries, BlockIndicesOut);
        invM.ExtractEntryView(blockM);

        SSolver->SetMatrix(*blockM);
        //cout << "PID:" << Comm.MyPID() << " localPIX:" << i << " globalPIX:" << PixMyGlobalElements[i] << endl;
        if ((*blockM)(0,0) > 0) {
            //rcond_blockM = PixMyGlobalElements[i];
            err = SSolver->ReciprocalConditionEstimate(rcond_blockM);
            if (err != 0) {
                cout << "PID:" << PixMap.Comm().MyPID() << " LocalRow[i]:" << i << " cannot compute RCOND, error code:" << err << " estimated:"<< rcond_blockM << endl;
                rcond_blockM = -2;
            }
            if (rcond_blockM > 1e-5) {
                err = SSolver->Invert();
                if (err != 0) {
                    cout << "PID:" << PixMap.Comm().MyPID() << " LocalRow[i]:" << i << " cannot invert matrix, error code:" << err << endl;
                    rcond_blockM = -3;
                }
            }
        } else {
            rcond_blockM = -1;
        }
        if (rcond_blockM < 1e-5) {
                for (int r=0; r<blockM->M(); ++r) {
                    for (int c=0; c<blockM->N(); ++c) {
                        (*blockM)(r,c) = 0.;
                    }
                }
        }
        RCondValues[0] = rcond_blockM;
        RCondIndices[0] = PixMyGlobalElements[i];
        rcond.ReplaceGlobalValues(1, 0, RCondValues, RCondIndices);
    }
}
#endif
