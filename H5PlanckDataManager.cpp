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
#include "Utils.h"

using boost::format;
using namespace std;
using namespace H5;

H5PlanckDataManager::H5PlanckDataManager(int FirstOd,int  LastOd, vector<string> Channels, string DataPath, string PointingPath) : Channels(Channels),
    DataPath(DataPath), PointingPath(PointingPath) 
{
    //TODO compute instead of hardcode

    GlobalOffset = 0; //computed from first OD
    TotalLength = 100; //computed from last OD - first OD
    LengthPerChannel = 25;

};

int H5PlanckDataManager::getPointing(long iStart, int nElements, pointing_t* pointing){
     H5File file( PointingPath, H5F_ACC_RDONLY );

    string channel = Channels[iStart/LengthPerChannel];
    cout << "Channel " << channel << endl;
    int Offset = iStart % LengthPerChannel; 

    DataSet dataset = file.openDataSet( channel );

    //DATASPACE
    DataSpace dataspace = dataset.getSpace();
    hsize_t  offset[1];       // hyperslab offset in memory
    hsize_t  count[1];        // size of the hyperslab in memory
    offset[0] = Offset;
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

    CompType pointing_h5type( sizeof(pointing_t) );
    pointing_h5type.insertMember( MEMBER1, HOFFSET(pointing_t, pix), PredType::NATIVE_INT);
    pointing_h5type.insertMember( MEMBER2, HOFFSET(pointing_t, qw), PredType::NATIVE_DOUBLE);
    pointing_h5type.insertMember( MEMBER3, HOFFSET(pointing_t, uw), PredType::NATIVE_DOUBLE);

    dataset.read( pointing, pointing_h5type, memspace, dataspace );
}

int H5PlanckDataManager::getData(string what, long iStart, int nElements, double* data){

    return 0;
}
