#define BOOST_TEST_MODULE datamanager
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <boost/assign/list_of.hpp> // for 'list_of()'
#include <string>
#include "PlanckDataManager.h"

using namespace boost::assign; // bring 'list_of()' into scope

//____________________________________________________________________________//

struct setUp {
    string dataPath, pointingPath;
    list<string> onech, twoch;
    setUp() {
                //dataPath = "/global/homes/z/zonca/planck/data/mission/lfi_ops_dx4";
                dataPath = "/home/zonca/p/testdata/lfi_ops_dx4";
                //pointingPath = "/global/homes/z/zonca/p/pointing/dx4_1024_nest";
                pointingPath = "/home/zonca/p/testdata/dx4_1024_nest";
                onech = list_of( "LFI27S" );
                twoch = list_of( "LFI28M" )( "LFI28S" );
             }
    ~setUp()
             { }
};

BOOST_FIXTURE_TEST_SUITE(test_planckdatamanager_suite, setUp)

BOOST_AUTO_TEST_CASE( test_getdatasetlength )
{
    PlanckDataManager * dm;
    dm = new PlanckDataManager(92, 92, onech, dataPath, pointingPath);
    BOOST_CHECK_EQUAL( dm->getDatasetLength(), 2812782);
    dm = new PlanckDataManager(92, 92, twoch, dataPath, pointingPath);
    BOOST_CHECK_EQUAL( dm->getDatasetLength(), 2812782*2 );
    dm = new PlanckDataManager(92, 94, onech, dataPath, pointingPath);
    BOOST_CHECK_EQUAL( dm->getDatasetLength(), 2812782 + 2812879 + 2812945);
}

BOOST_AUTO_TEST_CASE( test_getdata )
{
    double* data;
    data = new double[10];
    PlanckDataManager dm1(92, 92, onech, dataPath, pointingPath);
    dm1.getData("data", 4, 10, data);
    BOOST_CHECK_CLOSE(data[0], 0.001275042424, 1e-6);
    BOOST_CHECK_CLOSE(data[9], 0.001514075488, 1e-6);
}

BOOST_AUTO_TEST_CASE( test_getdata_twoch )
{
    double* data;
    data = new double[10];
    PlanckDataManager dm1(92, 92, twoch, dataPath, pointingPath);
    dm1.getData("data", 2812782, 10, data);
    BOOST_CHECK_CLOSE(data[0], -0.000420763029471, 1e-6);
    BOOST_CHECK_CLOSE(data[9], -0.001445570012144, 1e-6);
}

BOOST_AUTO_TEST_CASE( test_getqw_twoch )
{
    double* data;
    data = new double[10];
    PlanckDataManager dm1(92, 92, twoch, dataPath, pointingPath);
    dm1.getData("qw", 2812782, 10, data);
    BOOST_CHECK_CLOSE(data[0], 0.40361079, 1e-5);
    BOOST_CHECK_CLOSE(data[9], 0.34403955, 1e-5);
}

BOOST_AUTO_TEST_CASE( test_getpointing )
{
    double* data;
    data = new double[10];
    PlanckDataManager dm1(92, 92, twoch, dataPath, pointingPath);
    dm1.getData("pointing", 4, 10, data);
    BOOST_CHECK_CLOSE(data[0], 688111, 1e-6);
    BOOST_CHECK_CLOSE(data[9], 709496, 1e-6);
}

BOOST_AUTO_TEST_CASE( test_getqw )
{
    double* data;
    data = new double[10];
    PlanckDataManager dm1(92, 92, twoch, dataPath, pointingPath);
    dm1.getData("qw", 4, 10, data);
    BOOST_CHECK_CLOSE(data[0], 0.37720931, 1e-5);
    BOOST_CHECK_CLOSE(data[9], 0.31748405, 1e-5);
}

BOOST_AUTO_TEST_CASE( test_getuw )
{
    double* data;
    data = new double[10];
    PlanckDataManager dm1(92, 92, twoch, dataPath, pointingPath);
    dm1.getData("uw", 4, 10, data);
    BOOST_CHECK_CLOSE(data[0], 0.92612803, 1e-5);
    BOOST_CHECK_CLOSE(data[9], 0.94826361, 1e-5);
}

BOOST_AUTO_TEST_CASE( test_getdata_across )
{
    double* data;
    data = new double[10];
    PlanckDataManager dm1(92, 93, onech, dataPath, pointingPath);
    dm1.getData("data", 2812780, 10, data);
    BOOST_CHECK_CLOSE(data[0], -0.000309952667485, 1e-6); //5 from end of  elem of 92
    BOOST_CHECK_CLOSE(data[9], -0.001144992688151, 1e-6);  //8th elem of 93
}

BOOST_AUTO_TEST_CASE( test_qw_across )
{
    double* data;
    data = new double[10];
    PlanckDataManager dm1(92, 93, onech, dataPath, pointingPath);
    dm1.getData("qw", 2812780, 10, data);
    BOOST_CHECK_CLOSE(data[0], 0.581214755555, 1e-6); //5 from end of  elem of 92
    BOOST_CHECK_CLOSE(data[9], 0.582972722288, 1e-6);  //8th elem of 93
}
BOOST_AUTO_TEST_SUITE_END()
