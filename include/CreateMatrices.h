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
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include <Epetra_SerialDenseSolver.h>
#include "Epetra_Time.h"

#include "H5PlanckDataManager.h"
#include "Utils.h"

using namespace std;

int createGraph(const Epetra_Map& Map, const Epetra_Map& PixMap, const Epetra_IntVector & pix, Epetra_CrsGraph* &Graph) {

    int * MyGlobalElements = Map.MyGlobalElements();

    int MyPID = Map.Comm().MyPID();
    int NumMyElements = Map.NumMyElements();
    int err;

    boost::scoped_array<double> Values(new double[1]);
    int Indices[1];

    Graph = new Epetra_CrsGraph(Copy, Map, 1, true);

    log(MyPID, "Assembling Graph");
    for( int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local rows
            Indices[0] = pix[i];
            //cout << "0: " <<Indices[0] << " 1: " << Indices[1] << " 2: " << Indices[2] << endl;
            Graph->InsertGlobalIndices(MyGlobalElements[i], 1, Indices);
    }

    log(MyPID, "FillComplete Graph");
    Graph->FillComplete(PixMap, Map);
    log(MyPID, "FillComplete Graph done");
    Graph->OptimizeStorage();
    log(MyPID, "OptimizeStorage Graph done");
    return 0;
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

void createM(const Epetra_BlockMap& PixMap, const Epetra_BlockMap& Map, const Epetra_CrsMatrix* PQ, const Epetra_CrsMatrix* PU, double weight, int NSTOKES, int NPIX, Epetra_FEVbrMatrix& invM) {

    int err;

    
    int MyPID = Map.Comm().MyPID();

    int * PixMyGlobalElements = PixMap.MyGlobalElements();
    int NumEntries;
    int LocalPix[1];
    double Values[1];

    boost::scoped_array<int> GlobalPix (new int[1]);
    boost::scoped_ptr<Epetra_SerialDenseMatrix> Mpp (new Epetra_SerialDenseMatrix(NSTOKES, NSTOKES));
    boost::scoped_ptr<Epetra_SerialDenseMatrix> Prow (new Epetra_SerialDenseMatrix(1, NSTOKES));
    (*Prow)(0,0) = 1.;

    for( int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local pointing

        PQ->ExtractMyRowCopy(i, 1, NumEntries, Values, LocalPix);
        GlobalPix[0] = PQ->GCID(LocalPix[0]);

        if (GlobalPix[0] != NPIX) { 

            (*Prow)(0,1) = Values[0];

            //cout << Map.Comm().MyPID() << ": i:" << i << " Loc:" << LocalPix[0] << " Glob:" << GlobalPix[0] << " " << endl;

            PU->ExtractMyRowCopy(i, 1, NumEntries, Values);
            (*Prow)(0,2) = Values[0];


            err = Mpp->Multiply('T','N', weight, *Prow, *Prow, 0.);
                if (err != 0) {
                    cout << "Error in computing Mpp, error code:" << err << endl;
                    }

            invM.BeginSumIntoGlobalValues(GlobalPix[0], 1, GlobalPix.get());

            err = invM.SubmitBlockEntry(Mpp->A(), Mpp->LDA(), NSTOKES, NSTOKES);
                    if (err != 0) {
                        cout << "PID:" << PixMap.Comm().MyPID() << "Error in inserting values in M, error code:" << err << endl;
                        }

            err = invM.EndSubmitEntries();
                    if (err != 0) {
                        cout << "PID:" << PixMap.Comm().MyPID() << " LocalRow[i]:" << i << " Error in ending submit entries in M, error code:" << err << endl;
                        }

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

void invertM(const Epetra_Map& PixMap, int NSTOKES, Epetra_MultiVector& M, Epetra_Vector& rcond, Epetra_MultiVector& summap) {

    int * PixMyGlobalElements = PixMap.MyGlobalElements();
    boost::scoped_ptr<Epetra_SerialDenseSolver> SSolver (new Epetra_SerialDenseSolver());
    boost::scoped_ptr<Epetra_SerialSymDenseMatrix> blockM (new Epetra_SerialSymDenseMatrix());
    boost::scoped_ptr<Epetra_SerialDenseMatrix> PixelArray (new Epetra_SerialDenseMatrix(NSTOKES, 1));

    double rcond_blockM;
    blockM->Shape(NSTOKES);
    blockM->SetUpper();

    int RCondIndices[1], err, i_M;
    double RCondValues[1];

    int RowDim, NumBlockEntries;
    int *BlockIndicesOut;
    for( int i=0 ; i<PixMap.NumMyElements(); ++i ) { //loop on local pointing

        //build blockM
        i_M = 0;
        for (int j=0; j<NSTOKES; ++j) {
            for (int k=j; k<NSTOKES; ++k) {
                (*blockM)(j, k) = M[i_M][i];
            }
        }

        SSolver->SetMatrix(*blockM);
        //cout << "PID:" << MyPID << " localPIX:" << i << " globalPIX:" << PixMyGlobalElements[i] << endl;
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

        //apply to summap
        for (int j=0; j<NSTOKES; ++j) {
            PixelArray(j, 0) = summap[j][i];
        }

        blockM->Apply(*PixelArray, *PixelArray);

        //apply to summap
        for (int j=0; j<NSTOKES; ++j) {
            summap[j][i] = PixelArray(j, 0);
        }
        
    }
}
#endif
