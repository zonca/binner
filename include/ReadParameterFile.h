#ifndef READPARAMETERFILE_H
#define READPARAMETERFILE_H

#include <boost/assign/list_of.hpp> // for 'list_of()'
#include "PlanckDataManager.h"

using namespace boost::assign; // bring 'list_of()' into scope

int readParameterFile(string parameterFilename, PlanckDataManager *& dm) {

    //TODO read parameter file

    bool DEBUG = true;
    int NSIDE;
    string dataPath, pointingPath;

    dataPath = "/home/zonca/p/testdata/lfi_ops_dx4";
    //dataPath = "/global/homes/z/zonca/planck/data/mission/lfi_ops_dx4";

    if (DEBUG) {
        pointingPath = "/home/zonca/p/testdata/dx4_1_nest";
    //        pointingPath = "/global/homes/z/zonca/p/pointing/dx4_1_nest";
        NSIDE = 1;
    } else {
        pointingPath = "/global/homes/z/zonca/p/pointing/dx4_1024_nest";
        NSIDE = 1024;
    }
    //const list<string> channels = list_of( "LFI28M" );
    const list<string> channels = list_of("LFI27M" )( "LFI27S" )( "LFI28M" )( "LFI28S" );
    dm = new PlanckDataManager(92, 93, channels, dataPath, pointingPath);
    if (DEBUG) {
        dm->setDatasetLength(1000);
    }
    dm->NSIDE = NSIDE;
    dm->NSTOKES = 3;


}

#endif

