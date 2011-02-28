#ifndef POINTINGSTRUCT_H
#define POINTINGSTRUCT_H

#include <string>
#include "H5Cpp.h"

using namespace std;
using namespace H5;

const H5std_string MEMBER1( "PIX" );
const H5std_string MEMBER2( "QW" );
const H5std_string MEMBER3( "UW" );

typedef struct pointing_t {
	int    pix;
	double  qw;
	double uw;
} pointing_t;


#endif
