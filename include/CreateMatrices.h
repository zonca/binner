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
#include <Epetra_SerialDenseSolver.h>
#include "Epetra_Time.h"
#include "Epetra_Operator.h"

#include "DataManager.h"
#include "Utils.h"

using namespace std;


int SumMap(const Epetra_CrsMatrix * P, const DataManager * dm, const Epetra_MultiVector * TOD, const Epetra_MultiVector * yqu, Epetra_MultiVector * summap, Epetra_Vector * tempvec, Epetra_Vector * tempmap) {
    summap->PutScalar(0.);
    int qu_i=0;
    for (int i=1; i<3; ++i) { // Q=1 U=2
        tempvec->Multiply(1., *TOD, *((*yqu)(i)), 0.);
        P->Multiply1(true,*tempvec,*tempmap); //SUMMAP = Pt y
        qu_i = i-1;
        (*summap)(qu_i)->Update(1., *tempmap, 1.);
    }
}

int SumMap(const Epetra_CrsMatrix * P, const DataManager * dm, const Epetra_MultiVector * TOD, const Epetra_MultiVector * yqu, Epetra_MultiVector * summap) {
    boost::scoped_ptr<Epetra_Vector> tempvec (new Epetra_Vector(yqu->Map()));
    boost::scoped_ptr<Epetra_Vector> tempmap (new Epetra_Vector(summap->Map()));
    SumMap(P, dm, TOD, yqu, summap, tempvec.get(), tempmap.get());
}

int WeightMap(const Epetra_MultiVector * M, DataManager * dm, Epetra_MultiVector * summap, Epetra_SerialDenseMatrix * blockM, Epetra_SerialDenseMatrix * PixelArray, Epetra_SerialDenseMatrix * WeightedPixelArray) {
    int i_M;
    for( int i=0 ; i<summap->Map().NumMyElements(); ++i ) { //loop on local pointing
        //build blockM
        for (int j=0; j<dm->NSTOKES; ++j) {
            for (int k=j; k<dm->NSTOKES; ++k) {
                i_M = dm->getIndexM(j, k);
                (*blockM)(j, k) = (*M)[i_M][i];
                (*blockM)(k, j) = (*M)[i_M][i];
            }
        }
        for (int j=0; j<dm->NSTOKES; ++j) {
            (*PixelArray)(j, 0) = (*summap)[j][i];
        }
        blockM->Apply(*PixelArray, *WeightedPixelArray);
        for (int j=0; j<dm->NSTOKES; ++j) {
            (*summap)[j][i] = (*WeightedPixelArray)(j, 0);
        }
    }
}

int WeightMap(const Epetra_MultiVector * M, DataManager * dm, Epetra_MultiVector * summap) {
    boost::scoped_ptr<Epetra_SerialDenseMatrix> blockM (new Epetra_SerialDenseMatrix(dm->NSTOKES, dm->NSTOKES));
    boost::scoped_ptr<Epetra_SerialDenseMatrix> PixelArray (new Epetra_SerialDenseMatrix(dm->NSTOKES, 1));
    boost::scoped_ptr<Epetra_SerialDenseMatrix> WeightedPixelArray (new Epetra_SerialDenseMatrix(dm->NSTOKES, 1));
    WeightMap(M, dm, summap, blockM.get(), PixelArray.get(), WeightedPixelArray.get());
}

class DestripingOperator : public Epetra_Operator
{
public:
    DestripingOperator(const Epetra_CrsMatrix * P, const Epetra_MultiVector * yqu, const Epetra_MultiVector * M, const Epetra_CrsMatrix * F, DataManager * dm, Epetra_Map & BaselinesMap, Epetra_Vector * tempvec, Epetra_Vector * tempvec2, Epetra_Vector * tempmap, Epetra_MultiVector * summap):
        P(P), yqu(yqu), M(M), F(F), dm(dm), Map_(BaselinesMap), tempvec(tempvec), tempvec2(tempvec2), tempmap(tempmap), summap(summap)
    { 
    PixelArray = new Epetra_SerialDenseMatrix(dm->NSTOKES, 1);
    WeightedPixelArray = new Epetra_SerialDenseMatrix(dm->NSTOKES, 1);
    blockM = new Epetra_SerialDenseMatrix(dm->NSTOKES, dm->NSTOKES);
}

    int EstimateNoise( Epetra_MultiVector * TOD, Epetra_MultiVector & Y ) const;

    int Apply( const Epetra_MultiVector & X,
	     Epetra_MultiVector & Y ) const
        {
        //Fa
        F->Multiply(false,X,*tempvec2); //baselines TOD
        EstimateNoise(tempvec2, Y);
        }

  // other function
  int SetUseTranspose( bool UseTranspose) 
  {
    return(-1); // not implemented
  }

  int ApplyInverse( const Epetra_MultiVector & X,
		    Epetra_MultiVector & Y ) const
  {
     cout << "CALLING APPLY INVERSE" << endl;
    return(-1); // not implemented
  }

  double NormInf() const
  {
    return(-1);
  }

  const char * Label () const
  {
    return("DestripingOperator");
  }

  bool UseTranspose() const
  {
    return(false);
  }

  bool HasNormInf () const
  {
    return(false);
  }
  
  
  const Epetra_Comm & Comm() const
  {
    return(Map_.Comm());
  }

  const Epetra_Map & OperatorDomainMap() const
  {
    return(Map_);
  }
  
  const Epetra_Map & OperatorRangeMap() const
  {
    return(Map_);
  }

private:
    const Epetra_CrsMatrix * P; 
    const Epetra_MultiVector * yqu;
    const Epetra_MultiVector * M;
    const Epetra_CrsMatrix * F;
    Epetra_Map Map_;
    DataManager * dm;
    Epetra_Vector * tempvec;
    Epetra_Vector * tempvec2;
    Epetra_Vector * tempmap;
    Epetra_MultiVector * summap;
    Epetra_SerialDenseMatrix * PixelArray;
    Epetra_SerialDenseMatrix * blockM;
    Epetra_SerialDenseMatrix * WeightedPixelArray;
};

int DestripingOperator::EstimateNoise( Epetra_MultiVector * TOD,
	     Epetra_MultiVector & Y ) const
  {
    //// PFa
    SumMap(P, dm, TOD, yqu, summap, tempvec, tempmap);
    //M-1P-1Fa   BFa
    cout << dm->NSTOKES << endl;
    WeightMap(M, dm, summap, blockM, PixelArray, WeightedPixelArray);
    // Fa - PBFa
    // Q
    P->Multiply(false, *((*summap)(0)), *tempvec);
    TOD->Multiply(-1., *((*yqu)(1)), *tempvec, 1.);
    // U
    P->Multiply(false, *((*summap)(1)), *tempvec);
    TOD->Multiply(-1., *((*yqu)(2)), *tempvec, 1.);
    // F-1(Fa - PBFa)
    F->Multiply(true,*TOD,Y);
  }


int createGraph(const Epetra_Map& Map, const Epetra_Map& PixMap, const Epetra_Vector & pix, Epetra_CrsGraph* &Graph) {

    int * MyGlobalElements = Map.MyGlobalElements();

    int MyPID = Map.Comm().MyPID();
    int NumMyElements = Map.NumMyElements();
    int err;

    boost::scoped_array<double> Values(new double[1]);
    int Indices[1];

    Graph = new Epetra_CrsGraph(Copy, Map, 1, true);

    log(MyPID, "Assembling Graph");
    for( int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local rows
            Indices[0] = int(pix[i]);
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

void invertM(const Epetra_Map& PixMap, DataManager* dm, Epetra_MultiVector& M, Epetra_Vector& rcond, Epetra_MultiVector& summap) {

    int MyPID = PixMap.Comm().MyPID();
    int * PixMyGlobalElements = PixMap.MyGlobalElements();
    boost::scoped_ptr<Epetra_SerialDenseSolver> SSolver (new Epetra_SerialDenseSolver());
    boost::scoped_ptr<Epetra_SerialDenseMatrix> blockM (new Epetra_SerialDenseMatrix(dm->NSTOKES, dm->NSTOKES));
    boost::scoped_ptr<Epetra_SerialDenseMatrix> PixelArray (new Epetra_SerialDenseMatrix(dm->NSTOKES, 1));
    boost::scoped_ptr<Epetra_SerialDenseMatrix> WeightedPixelArray (new Epetra_SerialDenseMatrix(dm->NSTOKES, 1));

    double rcond_blockM;

    int RCondIndices[1], err, i_M;
    double RCondValues[1];

    //double * A;
    //int LDA; // stride
    //summap.ExtractView(&A, &LDA);

    for( int i=0 ; i<PixMap.NumMyElements(); ++i ) { //loop on local pointing

        //build blockM
        for (int j=0; j<dm->NSTOKES; ++j) {
            for (int k=j; k<dm->NSTOKES; ++k) {
                i_M = dm->getIndexM(j, k);
                (*blockM)(j, k) = M[i_M][i];
                (*blockM)(k, j) = M[i_M][i];
            }
        }


        SSolver->SetMatrix(*blockM);
        if ((*blockM)(0,0) > 0) {
            if ((PixMyGlobalElements[i] % 100000) == 0) {
                cout << "PID:" << MyPID << " localPIX:" << i << " globalPIX:" << PixMyGlobalElements[i] << *blockM << endl;
            }
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
        for (int j=0; j<dm->NSTOKES; ++j) {
            (*PixelArray)(j, 0) = summap[j][i];
        }
        
        blockM->Apply(*PixelArray, *WeightedPixelArray);
        if ((PixMyGlobalElements[i] % 100000) == 0) {
            cout << "PID:" << MyPID << " localPIX:" << i << " globalPIX:" << PixMyGlobalElements[i] << *PixelArray << endl;
            cout << "PID:" << MyPID << " localPIX:" << i << " globalPIX:" << PixMyGlobalElements[i] << *blockM << endl;
            cout << "PID:" << MyPID << " localPIX:" << i << " globalPIX:" << PixMyGlobalElements[i] << *WeightedPixelArray << endl;
        }

        for (int j=0; j<dm->NSTOKES; ++j) {
            summap[j][i] = (*WeightedPixelArray)(j, 0);
        }
        
        // set inverted values back to M
        for (int j=0; j<dm->NSTOKES; ++j) {
            for (int k=j; k<dm->NSTOKES; ++k) {
                i_M = dm->getIndexM(j, k);
                M[i_M][i] = (*blockM)(j, k);
            }
        }
    }
}
#endif
