#define BOOST_TEST_MODULE test_read_fits_lib
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <cstring>

extern "C" {
#include "read_fits_lib.h"
}

using namespace std;

struct setUp {
    string fileData1, fileData2, filePix1, filePix2;
    setUp() {
                fileData1 = "data/test_qu/01/data.fits";
                fileData2 = "data/test_qu/02/data.fits";
                filePix1 = "data/test_qu/01/pix.fits";
                filePix2 = "data/test_qu/02/pix.fits";
             }
    ~setUp()
             {
             }
};

BOOST_FIXTURE_TEST_SUITE(test_read_fits_lib_suite, setUp)

BOOST_AUTO_TEST_CASE( test_nelem )
{
    long length = 0;
    nelem(fileData1.c_str(), &length);
    BOOST_CHECK_EQUAL( length, 30);
    nelem(fileData2.c_str(), &length);
    BOOST_CHECK_EQUAL( length, 18);
}

BOOST_AUTO_TEST_CASE( test_read_data )
{
    double* d;
    d = new double[10];
    read_data(fileData1.c_str(), "DATA", 2, 2, 10, d);
    BOOST_CHECK_CLOSE( d[0], 1.5464101615137755, 1e-9);
    BOOST_CHECK_CLOSE( d[3], -3.036602540378444, 1e-9);
    BOOST_CHECK_CLOSE( d[9], 2.4188840424654723, 1e-9);
}

BOOST_AUTO_TEST_CASE( test_read_pix )
{
    double* d;
    d = new double[5];
    read_data(filePix2.c_str(), "PIX", 2, 5, 5, d);
    BOOST_CHECK_EQUAL( d[0], 8);
    BOOST_CHECK_EQUAL( d[4], 3);
}

BOOST_AUTO_TEST_CASE( test_read_qw )
{
    double* d;
    d = new double[6];
    read_data(filePix1.c_str(), "PIX", 3, 20, 6, d);
    BOOST_CHECK_CLOSE( d[0], 0.49999999999999978, 1e-9);
    BOOST_CHECK_CLOSE( d[5], -0.64278760968653925, 1e-9);
}

BOOST_AUTO_TEST_CASE( test_read_qu )
{
    double* d;
    d = new double[6];
    read_data(filePix2.c_str(), "PIX", 4, 0, 6, d);
    BOOST_CHECK_CLOSE( d[0], -0.7660444431189779, 1e-9);
    BOOST_CHECK_CLOSE( d[5], -0.1736481776669303, 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()
