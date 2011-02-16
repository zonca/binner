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
#include "DataManager.h"
#include "TestHelpers.h"

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
    int NumElements;
    Epetra_Map* Map;
    setUp() {
                NumElements = 8; 
                Epetra_MpiComm Comm(MPI_COMM_WORLD);
                Map = new Epetra_Map(NumElements, 0, Comm);
             }
    ~setUp()
             { }
};

BOOST_GLOBAL_FIXTURE(MPIsetUp)

BOOST_FIXTURE_TEST_SUITE(test_creatematrices_suite, setUp)

BOOST_AUTO_TEST_CASE( test_numelements )
{
    BOOST_CHECK_EQUAL( Map->NumGlobalElements(), NumElements);
}

BOOST_AUTO_TEST_CASE( test_inverm )
{
    Epetra_CrsMatrix M(Copy, *Map,2);
    createFakeM(*Map, M);
    //cout << M;
    invertM(*Map, M);
    cout << M;
}

BOOST_AUTO_TEST_SUITE_END()
