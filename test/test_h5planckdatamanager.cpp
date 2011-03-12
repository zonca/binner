#define BOOST_TEST_MODULE datamanager
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <string>

#include "H5PlanckDataManager.h"

typedef std::map<string, double> WeightDict;

//____________________________________________________________________________//

struct setUp {
    string dataPath, pointingPath;
    vector<string> onech, twoch;
    WeightDict Weights;
    setUp() {
                pointingPath = "/home/zonca/p/testdata/dx4_1024_nest_30_9293.h5";
                dataPath = "/home/zonca/p/testdata/lfi_ops_dx4.h5";
                onech.push_back("LFI27S");

                twoch.push_back("LFI28M");
                twoch.push_back("LFI28S");
                Weights["LFI18"] = 5.2814E+04;
                Weights["LFI19"] = 3.9294E+04;
                Weights["LFI23"] = 4.6195E+04;
                Weights["LFI22"] = 4.8167E+04;
                Weights["LFI20"] = 3.4468E+04;
                Weights["LFI21"] = 4.8501E+04;
                Weights["LFI24"] = 1.1531E+05;
                Weights["LFI25"] = 1.3288E+05;
                Weights["LFI26"] = 1.0654E+05;
                Weights["LFI27"] = 3.6656E+05;
                Weights["LFI28"] = 3.4432E+05;
             }
    ~setUp()
             { }
};

BOOST_FIXTURE_TEST_SUITE(test_planckdatamanager_suite, setUp)

BOOST_AUTO_TEST_CASE( test_getdatasetlength )
{
    H5PlanckDataManager * dm;
    dm = new H5PlanckDataManager(92, 92, onech, dataPath, pointingPath, Weights);
    BOOST_CHECK_EQUAL( dm->getDatasetLength(), 2812782);
    dm = new H5PlanckDataManager(92, 92, twoch, dataPath, pointingPath, Weights);
    BOOST_CHECK_EQUAL( dm->getDatasetLength(), 2812782*2 );
    dm = new H5PlanckDataManager(92, 93, onech, dataPath, pointingPath, Weights);
    BOOST_CHECK_EQUAL( dm->getDatasetLength(), 2812782 + 2812879);
}

BOOST_AUTO_TEST_CASE( test_getdata )
{
    double* data;
    data = new double[10];
    H5PlanckDataManager dm1(92, 92, onech, dataPath, pointingPath, Weights);
    dm1.getData("LFI27S", 4, 10, data);
    BOOST_CHECK_CLOSE(data[0], 0.001275042424, 1e-6);
    BOOST_CHECK_CLOSE(data[9], 0.001514075488, 1e-6);

    H5PlanckDataManager dm2(93, 93, onech, dataPath, pointingPath, Weights);
    dm2.getData("LFI28S", 0, 10, data);
    BOOST_CHECK_CLOSE(data[8], 0.0015338654352308992, 1e-6);
}

BOOST_AUTO_TEST_CASE( test_getdata_outofbounds )
{
    double* data;
    data = new double[10];
    H5PlanckDataManager dm1(92, 92, onech, dataPath, pointingPath, Weights);
    dm1.getData("LFI28S", 2812780, 10, data);
    BOOST_CHECK_CLOSE(data[0], -0.00091543088, 1e-6);
    BOOST_CHECK_CLOSE(data[1], 5.8370141748324746e-05, 1e-6);
    BOOST_CHECK_CLOSE(data[2], 0, 1e-6);
    BOOST_CHECK_CLOSE(data[3], 0, 1e-6);
}

BOOST_AUTO_TEST_CASE( test_getpointing )
{
    int pix[10];
    double qw[10], uw[10];
    H5PlanckDataManager dm1(92, 92, twoch, dataPath, pointingPath, Weights);
    dm1.getPointing("LFI28M", 4, 10, pix, qw, uw);
    BOOST_CHECK_EQUAL(pix[0], 688111);
    BOOST_CHECK_EQUAL(pix[9], 709496);
    BOOST_CHECK_CLOSE(qw[0], 0.37720931, 1e-5);
    BOOST_CHECK_CLOSE(qw[9], 0.31748405, 1e-5);
    BOOST_CHECK_CLOSE(uw[0], 0.92612803, 1e-5);
    BOOST_CHECK_CLOSE(uw[9], 0.94826361, 1e-5);
}

//BOOST_AUTO_TEST_CASE( test_getdata_across )
//{
//    double data[10];
//    H5PlanckDataManager dm1(92, 93, onech, dataPath, pointingPath);
//    dm1.getData(2812780, 10, data);
//    BOOST_CHECK_CLOSE(data[0], -0.000309952667485, 1e-6); //5 from end of  elem of 92
//    BOOST_CHECK_CLOSE(data[9], -0.001144992688151, 1e-6);  //8th elem of 93
//}
//
//BOOST_AUTO_TEST_CASE( test_qw_across )
//{
//    pointing_t pointing[10];
//    H5PlanckDataManager dm1(92, 93, onech, dataPath, pointingPath);
//    dm1.getPointing(2812780, 10, pointing);
//    BOOST_CHECK_CLOSE(pointing[0].qw, 0.581214755555, 1e-6); //5 from end of  elem of 92
//    BOOST_CHECK_CLOSE(pointing[9].qw, 0.582972722288, 1e-6);  //8th elem of 93
//}
BOOST_AUTO_TEST_SUITE_END()
