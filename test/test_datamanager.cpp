#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE datamanager
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <string>
#include "DataManager.h"

//____________________________________________________________________________//

struct setUp {
    Teuchos::Array<string> channels;
    setUp() {
                channels = Teuchos::tuple<string>("data");
             }
    ~setUp()
             { }
};

BOOST_FIXTURE_TEST_SUITE(test_datamanager_suite, setUp)

BOOST_AUTO_TEST_CASE( test_getdatasetlength )
{
    DataManager dm1(1, 1, channels, "data/test_qu","data/test_qu");
    BOOST_CHECK_EQUAL( dm1.getLength(), 30);
    DataManager dm2(2, 2, channels, "data/test_qu","data/test_qu");
    BOOST_CHECK_EQUAL( dm2.getLength(), 18);
    DataManager dm12(1, 2, channels, "data/test_qu","data/test_qu");
    BOOST_CHECK_EQUAL( dm12.getLength(), 48);
}

BOOST_AUTO_TEST_CASE( test_getdata )
{
    double* data;
    data = new double[10];
    DataManager dm1(1, 1, channels, "data/test_qu","data/test_qu");
    dm1.getData("data", 4, 10, data);
    BOOST_CHECK_CLOSE(data[0], -3.8111595753452776, 1e-9);
    BOOST_CHECK_CLOSE(data[9], -1.2044160264027586, 1e-9);
}

BOOST_AUTO_TEST_CASE( test_getdata_across )
{
    double* data;
    data = new double[10];
    DataManager dm1(1, 2, channels, "data/test_qu","data/test_qu");
    dm1.getData("data", 28, 10, data);
    BOOST_CHECK_CLOSE(data[0], 1.8415440693208618, 1e-9); //28th elem of 01
    BOOST_CHECK_CLOSE(data[9], -0.43837158328378062, 1e-9);  //7th elem of 02
}

BOOST_AUTO_TEST_CASE( test_qw_across )
{
    double* data;
    data = new double[10];
    DataManager dm1(1, 2, channels, "data/test_qu","data/test_qu");
    dm1.getData("qw", 28, 10, data);
    BOOST_CHECK_CLOSE(data[0], -0.98480775301220802, 1e-9); //28th elem of 01
    BOOST_CHECK_CLOSE(data[9], -0.64278760968653925, 1e-9);  //7th elem of 02
}
BOOST_AUTO_TEST_CASE( test_adjustdistribution )
{
    DataManager dm(1, 2, channels, "data/test_qu","data/test_qu");
    dm.BaselineLength = 4;
    BOOST_CHECK_EQUAL(dm.adjustDistribution(0,4), 4);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(0,5), 8);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(0,6), 8);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(0,7), 8);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(0,8), 8);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(0,9), 12);

    BOOST_CHECK_EQUAL(dm.adjustDistribution(2,3), 4);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(3,3), 4);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(2,7), 8);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(2,8), 8);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(2,10), 8);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(2,11), 12);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(2,12), 12);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(2,14), 12);

    BOOST_CHECK_EQUAL(dm.adjustDistribution(30,4), 4);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(30,5), 8);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(30,6), 8);

    BOOST_CHECK_EQUAL(dm.adjustDistribution(39,9), 6);
}

BOOST_AUTO_TEST_CASE( test_adjustdistribution_across )
{
    DataManager dm(1, 2, channels, "data/test_qu","data/test_qu");
    dm.BaselineLength = 4;
    BOOST_CHECK_EQUAL(dm.adjustDistribution(26,5), 2+4);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(26,6), 2+4);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(26,7), 2+4);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(26,8), 2+4);
    BOOST_CHECK_EQUAL(dm.adjustDistribution(26,9), 2+8);
}
BOOST_AUTO_TEST_SUITE_END()
