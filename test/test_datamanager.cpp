#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE datamanager
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <string>
#include "DataManager.h"

//____________________________________________________________________________//

BOOST_AUTO_TEST_CASE( test_getdatasetlength )
{
    DataManager dm1(1, 1, "DATA", "data/test_qu","data/test_qu");
    BOOST_CHECK_EQUAL( dm1.getDatasetLength(), 30);
    DataManager dm2(2, 2, "DATA", "data/test_qu","data/test_qu");
    BOOST_CHECK_EQUAL( dm2.getDatasetLength(), 18);
    DataManager dm12(1, 2, "DATA", "data/test_qu","data/test_qu");
    BOOST_CHECK_EQUAL( dm12.getDatasetLength(), 48);
}

BOOST_AUTO_TEST_CASE( test_getdata )
{
    double* data;
    data = new double[10];
    DataManager dm1(1, 1, "DATA", "data/test_qu","data/test_qu");
    dm1.getData("data", 4, 10, data);
    BOOST_CHECK_CLOSE(data[0], -3.8111595753452776, 1e-9);
    BOOST_CHECK_CLOSE(data[9], -1.2044160264027586, 1e-9);
}

BOOST_AUTO_TEST_CASE( test_getdata_across )
{
    double* data;
    data = new double[10];
    DataManager dm1(1, 2, "DATA", "data/test_qu","data/test_qu");
    dm1.getData("data", 28, 10, data);
    BOOST_CHECK_CLOSE(data[0], -2.1584559306791382, 1e-9); //28th elem of 01
    BOOST_CHECK_CLOSE(data[9], 1.5616284167162193, 1e-9);  //7th elem of 02
}

BOOST_AUTO_TEST_CASE( test_qw_across )
{
    double* data;
    data = new double[10];
    DataManager dm1(1, 2, "DATA", "data/test_qu","data/test_qu");
    dm1.getData("qw", 28, 10, data);
    BOOST_CHECK_CLOSE(data[0], -0.98480775301220802, 1e-9); //28th elem of 01
    BOOST_CHECK_CLOSE(data[9], -0.64278760968653925, 1e-9);  //7th elem of 02
}
