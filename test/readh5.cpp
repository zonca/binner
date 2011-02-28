#include "H5Cpp.h"
#include <iostream>

using namespace H5;
using namespace std;

const H5std_string MEMBER1( "PIX" );
const H5std_string MEMBER2( "QW" );
const H5std_string MEMBER3( "UW" );

int main (void)
{
    H5File file( "/home/zonca/p/testdata/dx4_1024_nest_30_9293.h5", H5F_ACC_RDONLY );

    DataSet dataset = file.openDataSet( "OD" );
    H5T_class_t type_class = dataset.getTypeClass();

    DataSpace dataspace = dataset.getSpace();

    hsize_t length[1];
    dataspace.getSimpleExtentDims( length, NULL);
    cout << "Length: " << length[0] << endl;

    dataset = file.openDataSet( "LFI28M" );
    dataspace = dataset.getSpace();

    hsize_t dims_out[2];
    dataspace.getSimpleExtentDims( dims_out, NULL);
    cout << "dimensions " << (unsigned long)(dims_out[0]) <<" x " <<
(unsigned long)(dims_out[1]) << endl;


    
       /* First structure  and dataset*/
   typedef struct pointing_t {
	int    pix;
	double  qw;
	double uw;
   } pointing_t;

   int LENGTH=10;
   pointing_t pointing[LENGTH];

    CompType mtype1( sizeof(pointing_t) );
    mtype1.insertMember( MEMBER1, HOFFSET(pointing_t, pix), PredType::NATIVE_INT);
    mtype1.insertMember( MEMBER2, HOFFSET(pointing_t, qw), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( MEMBER3, HOFFSET(pointing_t, uw), PredType::NATIVE_DOUBLE);

    dataset = DataSet (file.openDataSet( "LFI28M" ));

    dataspace = dataset.getSpace();
    hsize_t      offset[1];       // hyperslab offset in memory
    hsize_t      count[1];        // size of the hyperslab in memory
    offset[0] = 1000;
    count[0]  = 10;
    dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );

    hsize_t     dimsm[1];
    dimsm[0] = 10;
    hsize_t      offset_out[1];       // hyperslab offset in memory
    offset_out[0] = 0;
    hsize_t      count_out[1];        // size of the hyperslab in memory
    count_out[0]  = 10;
    DataSpace memspace( 1, dimsm );
    memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out );


    dataset.read( pointing, mtype1, memspace, dataspace );

    for (int i=0; i<10; ++i) 
    {
        cout << pointing[i].pix << " " << pointing[i].qw << " " << pointing[i].uw << endl;
    }
}
