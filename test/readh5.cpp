#include "H5Cpp.h"
#include <iostream>

using namespace H5;
using namespace std;

int main (void)
{
    H5File file( "dx4_1024_nest_30_9293.h5", H5F_ACC_RDONLY );

    DataSet dataset = file.openDataSet( "OD" );
    H5T_class_t type_class = dataset.getTypeClass();

    DataSpace dataspace = dataset.getSpace();

    hsize_t length[1];
    dataspace.getSimpleExtentDims( length, NULL);
    cout << "Length: " << length[0] << endl;

    dataset = file.openDataSet( "LFI28M" );
    dataspace = dataset.getSpace();

    dataspace.getSimpleExtentDims( dims_out, NULL);
    

}
