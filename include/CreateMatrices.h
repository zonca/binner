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
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialSymDenseMatrix.h"
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

void invertM(const Epetra_Map& PixMap, int NSTOKES, Epetra_MultiVector& M, Epetra_Vector& rcond, Epetra_MultiVector& summap) {

    int MyPID = PixMap.Comm().MyPID();
    int * PixMyGlobalElements = PixMap.MyGlobalElements();
    boost::scoped_ptr<Epetra_SerialDenseSolver> SSolver (new Epetra_SerialDenseSolver());
    boost::scoped_ptr<Epetra_SerialSymDenseMatrix> blockM (new Epetra_SerialSymDenseMatrix());
    boost::scoped_ptr<Epetra_SerialDenseMatrix> PixelArray (new Epetra_SerialDenseMatrix(NSTOKES, 1));
    boost::scoped_ptr<Epetra_SerialDenseMatrix> WeightedPixelArray (new Epetra_SerialDenseMatrix(NSTOKES, 1));

    double rcond_blockM;
    blockM->Shape(NSTOKES);
    blockM->SetUpper();

    int RCondIndices[1], err, i_M;
    double RCondValues[1];

    //double * A;
    //int LDA; // stride
    //summap.ExtractView(&A, &LDA);

    for( int i=0 ; i<PixMap.NumMyElements(); ++i ) { //loop on local pointing

        //build blockM
        for (int j=0; j<NSTOKES; ++j) {
            for (int k=j; k<NSTOKES; ++k) {
                i_M = j * (2*NSTOKES-1 - j)/2 + k;
                (*blockM)(j, k) = M[i_M][i];
            }
        }


        SSolver->SetMatrix(*blockM);
        if ((*blockM)(0,0) > 0) {
            cout << "PID:" << MyPID << " localPIX:" << i << " globalPIX:" << PixMyGlobalElements[i] << *blockM << endl;
            //rcond_blockM = PixMyGlobalElements[i];
            err = SSolver->ReciprocalConditionEstimate(rcond_blockM);
            if (err != 0) {
                cout << "PID:" << PixMap.Comm().MyPID() << " LocalRow[i]:" << i << " cannot compute RCOND, error code:" << err << " estimated:"<< rcond_blockM << endl;
                rcond_blockM = -2;
            }
            if (rcond_blockM > 1e-7) {
                err = SSolver->Invert();
                if (err != 0) {
                    cout << "PID:" << PixMap.Comm().MyPID() << " LocalRow[i]:" << i << " cannot invert matrix, error code:" << err << endl;
                    rcond_blockM = -3;
                }

            }
        } else {
            rcond_blockM = -1;
        }
        if (rcond_blockM < 1e-7) {
            blockM->Scale(0.);
        }
        RCondValues[0] = rcond_blockM;
        RCondIndices[0] = PixMyGlobalElements[i];
        rcond.ReplaceGlobalValues(1, 0, RCondValues, RCondIndices);

        //apply to summap
        for (int j=0; j<NSTOKES; ++j) {
            (*PixelArray)(j, 0) = summap[j][i];
        }
        
        blockM->Apply(*PixelArray, *WeightedPixelArray);

        //apply to summap
        for (int j=0; j<NSTOKES; ++j) {
            summap[j][i] = (*WeightedPixelArray)(j, 0);
        }
        
    }
}
#endif
