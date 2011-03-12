#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <cstring>
#include <vector>

#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "H5Cpp.h"

#include "H5PlanckDataManager.h"


using boost::format;
using namespace std;
using namespace H5;

H5PlanckDataManager::H5PlanckDataManager(int FirstOd,int  LastOd, vector<string> Channels, string DataPath, string PointingPath, WeightDict Weights) : Channels(Channels),
    DataPath(DataPath), PointingPath(PointingPath), Weights(Weights)
{
    H5File file( PointingPath, H5F_ACC_RDONLY );

    DataSet dataset = file.openDataSet( "OD" );

    //DATASPACE
    DataSpace dataspace = dataset.getSpace();
    hsize_t  offset[1];       // hyperslab offset in memory
    hsize_t  count[1];        // size of the hyperslab in memory
    offset[0] = FirstOd - 91;
    count[0]  = 1;
    dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );

    //MEMSPACE
    hsize_t     dimsm[1];
    dimsm[0] = 1;
    hsize_t      offset_out[1];       // hyperslab offset in memory
    offset_out[0] = 0;
    hsize_t      count_out[1];        // size of the hyperslab in memory
    count_out[0]  = 1;
    DataSpace memspace( 1, dimsm );
    memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out );

    int odstart[1];
    dataset.read( odstart, PredType::NATIVE_INT, memspace, dataspace );

    GlobalOffset = odstart[0]; //computed from first OD

    offset[0] = LastOd + 1 - 91;
    dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );
    dataset.read( odstart, PredType::NATIVE_INT, memspace, dataspace );

    LengthPerChannel = odstart[0];
    //LengthPerChannel = 3000000 * 6 / 4.;
    TotalLength = LengthPerChannel * Channels.size(); //computed from last OD - first OD

};

int H5PlanckDataManager::getPointing(string channel, long iStart, int nElements, int * pix, double * qw, double * uw){
     H5File file( PointingPath, H5F_ACC_RDONLY );

    int NPIX = getNPIX();
    if (iStart >= LengthPerChannel) {
        cout << "istart over " << iStart << endl;
    } else {

        //cout << "istart + nELEm " << iStart + nElements << endl;
        if (iStart + nElements > LengthPerChannel) {
            nElements = LengthPerChannel - iStart;
            cout << "istart " << iStart << " new nElem " << nElements << endl;
        }

        //cout << "Channel " << channel << endl;
        DataSet dataset = file.openDataSet( channel );

        //DATASPACE
        DataSpace dataspace = dataset.getSpace();
        hsize_t  offset[1];       // hyperslab offset in memory
        hsize_t  count[1];        // size of the hyperslab in memory
        offset[0] = GlobalOffset + iStart;
        count[0]  = nElements;
        dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );

        //MEMSPACE
        hsize_t     dimsm[1];
        dimsm[0] = nElements;
        hsize_t      offset_out[1];       // hyperslab offset in memory
        offset_out[0] = 0;
        hsize_t      count_out[1];        // size of the hyperslab in memory
        count_out[0]  = nElements;
        DataSpace memspace( 1, dimsm );
        memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out );

        CompType pointing_h5type( sizeof(int) );
        pointing_h5type.insertMember( MEMBER1, 0, PredType::NATIVE_INT);
        dataset.read(pix, pointing_h5type, memspace, dataspace );

        CompType pointing_h5typeqw( sizeof(double) );
        pointing_h5typeqw.insertMember( MEMBER2, 0, PredType::NATIVE_DOUBLE);
        dataset.read(qw, pointing_h5typeqw, memspace, dataspace );

        CompType pointing_h5typeuw( sizeof(double) );
        pointing_h5typeuw.insertMember( MEMBER3, 0, PredType::NATIVE_DOUBLE);
        dataset.read(uw, pointing_h5typeuw, memspace, dataspace );
    }
    return 0;
}

int H5PlanckDataManager::getData(string channel, long iStart, int nElements, double* data){
    H5File file( DataPath, H5F_ACC_RDONLY );

    if (iStart < LengthPerChannel) {
        if (iStart + nElements >= LengthPerChannel) {
            nElements = LengthPerChannel - iStart;
        }

        DataSet dataset = file.openDataSet( channel );

        //DATASPACE
        DataSpace dataspace = dataset.getSpace();
        hsize_t  offset[1];       // hyperslab offset in memory
        hsize_t  count[1];        // size of the hyperslab in memory
        offset[0] = GlobalOffset + iStart;
        count[0]  = nElements;
        dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );

        //MEMSPACE
        hsize_t     dimsm[1];
        dimsm[0] = nElements;
        hsize_t      offset_out[1];       // hyperslab offset in memory
        offset_out[0] = 0;
        hsize_t      count_out[1];        // size of the hyperslab in memory
        count_out[0]  = nElements;
        DataSpace memspace( 1, dimsm );
        memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out );

        dataset.read( data, PredType::NATIVE_DOUBLE, memspace, dataspace );
    }
    return 0;
}
