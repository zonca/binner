#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <boost/format.hpp>
#include "H5Cpp.h"

using namespace std;
using boost::format;

int log(int MyPID, string message) {

    if( MyPID == 0 ) {
        cout << "******** " << message << endl << endl << flush;
    }
    return true;
} 

int log(int MyPID, format message) {

    log(MyPID, message.str());
    return true;

} 

int WriteH5Vec(Epetra_Vector * vec, string filename, string outputFolder) {
    int MyPID = vec->Comm().MyPID();
    H5std_string  FILE_NAME( str( format(outputFolder + "/%s_%03d.h5") % filename % MyPID ) );
    H5File file(FILE_NAME, H5F_ACC_TRUNC );
    hsize_t dimsf[1];
    dimsf[0] = vec->Map().NumMyElements();
    DataSpace dataspace( 1, dimsf );
    DataSet dataset = file.createDataSet( "Vector", PredType::NATIVE_DOUBLE, dataspace );
    double * data;
    vec->ExtractView(&data);
    dataset.write( data, PredType::NATIVE_DOUBLE );
    return 0;
}

#endif
