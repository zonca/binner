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

#endif
