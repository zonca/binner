#define BOOST_TEST_MODULE datamanager
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <string>

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

#include "CreateMatrices.h"
#include "PlanckDataManager.h"
#include "TestHelpers.h"
#include "Epetra_FEVbrMatrix.h"

//____________________________________________________________________________//

using namespace std;

struct MPIsetUp {

    MPIsetUp() {
      #ifdef HAVE_MPI
      int dummy;
      char ** dummy2;
      MPI_Init(&dummy, &dummy2);
      #endif  
             }
    ~MPIsetUp()
             {
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
             }
};

struct setUp {
    int NumElements, NPIX, NSTOKES, MyPID;
    Epetra_Map* Map;
    Epetra_BlockMap* PixMap;
    setUp() {
                NumElements = 21; 
                NPIX = 5;
                NSTOKES = 3;
                Epetra_MpiComm Comm(MPI_COMM_WORLD);
                MyPID = Comm.MyPID();
                Map = new Epetra_Map(NumElements, 0, Comm);
                PixMap = new Epetra_BlockMap(NPIX, 3, 0, Comm);
             }
    ~setUp()
             { }
};

BOOST_GLOBAL_FIXTURE(MPIsetUp)

BOOST_FIXTURE_TEST_SUITE(test_creatematrices_suite, setUp)

BOOST_AUTO_TEST_CASE( test_numelements )
{
    BOOST_CHECK_EQUAL( Map->NumGlobalElements(), NumElements);
    BOOST_CHECK_EQUAL( PixMap->NumGlobalElements(), NPIX);
}

BOOST_AUTO_TEST_CASE( test_initm )
{
    int RowDim, NumBlockEntries;
    int * BlockIndicesOut;
    Epetra_SerialDenseMatrix * blockM;
    Epetra_FEVbrMatrix invM(Copy, *PixMap, 1);
    initM(*PixMap, NSTOKES, invM);
    invM.GlobalAssemble();
    Epetra_SerialDenseMatrix Zero(NSTOKES, NSTOKES);

    BOOST_CHECK_EQUAL( invM.NumGlobalBlockCols(), NPIX);
    BOOST_CHECK_EQUAL( invM.NumGlobalBlockRows(), NPIX);
    for( int i=0 ; i<PixMap->NumMyElements(); ++i ) { //loop on local pixel
        invM.BeginExtractMyBlockRowView(i, RowDim, NumBlockEntries, BlockIndicesOut);
        invM.ExtractEntryView(blockM);
        BOOST_CHECK_EQUAL( *blockM, Zero);
    }
}

BOOST_AUTO_TEST_CASE( test_createm )
{
    Epetra_VbrMatrix P(Copy,*Map,1);
    createFakeP(*Map, *PixMap, P);
    cout << P << endl;
    cout << "Ind Global " << P.IndicesAreGlobal() << endl;
    cout << "Ind Loc " << P.IndicesAreLocal() << endl;

    Epetra_FEVbrMatrix invM(Copy, *PixMap, 1);
    initM(*PixMap, NSTOKES, invM);
    createM(*PixMap, *Map, P, NSTOKES, invM);
    invM.GlobalAssemble();

    int RowDim, NumBlockEntries;
    int * BlockIndicesOut;
    Epetra_SerialDenseMatrix * blockM;

    Epetra_IntSerialDenseVector hits(PixMap->NumMyElements());

    if (MyPID == 0) {
        hits[0]=1;
        hits[1]=8;
        hits[2]=0;
    } else {
        hits[0]=6;
        hits[1]=6;
    }

    for( int i=0 ; i<PixMap->NumMyElements(); ++i ) { //loop on local pixel
        invM.BeginExtractMyBlockRowView(i, RowDim, NumBlockEntries, BlockIndicesOut);
        invM.ExtractEntryView(blockM);
        BOOST_CHECK_EQUAL( (*blockM)(0,0), hits[i]);
    }
}

//BOOST_AUTO_TEST_CASE( test_invertm )
//{
//    Epetra_VbrMatrix P(Copy,*Map,1);
//    createFakeP(*Map, *PixMap, P);
//    Epetra_FEVbrMatrix invM(Copy, *PixMap, 1);
//    initM(*PixMap, NSTOKES, invM);
//    createM(*PixMap, *Map, P, NSTOKES, invM);
//    invM.GlobalAssemble();
//
//    Epetra_Vector rcond(*PixMap);
//    cout << invM << endl;
//    invertM(*PixMap, invM, rcond);
//    cout << invM << endl;
//    PixMap->Comm().Barrier();
//
//    int RowDim, NumBlockEntries;
//    int * BlockIndicesOut;
//    Epetra_SerialDenseMatrix * blockM;
//
//    invM.BeginExtractMyBlockRowView(0, RowDim, NumBlockEntries, BlockIndicesOut);
//    invM.ExtractEntryView(blockM);
//
//    if (MyPID == 0) {
//        BOOST_CHECK_CLOSE( (*blockM)(0,0), 0.67217016, 0.003);
//        BOOST_CHECK_CLOSE( (*blockM)(0,1), -0.17267239, 0.003);
//        BOOST_CHECK_CLOSE( (*blockM)(1,0), -0.17267239, 0.003);
//        BOOST_CHECK_CLOSE( (*blockM)(1,1), 0.56853723, 0.003);
//    }
//
//}
BOOST_AUTO_TEST_SUITE_END()
