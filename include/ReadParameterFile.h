#ifndef READPARAMETERFILE_H
#define READPARAMETERFILE_H

#include <boost/assign/list_of.hpp> // for 'list_of()'
#include "H5PlanckDataManager.h"

using namespace boost::assign; // bring 'list_of()' into scope

void readParameterFile(string parameterFilename, H5PlanckDataManager *& dm) {

    //TODO read parameter file

    bool DEBUG = false;
    int NSIDE;
    string dataPath, pointingPath;

    dataPath = "/scratch/scratchdirs/zonca/pointing/lfi_ops_dx4_30.h5";

    if (DEBUG) {
        pointingPath = "/home/zonca/p/testdata/dx4_1_nest";
        pointingPath = "/global/homes/z/zonca/p/pointing/dx4_1_nest";
        NSIDE = 1;
    } else {
        pointingPath = "/scratch/scratchdirs/zonca/pointing/dx4_1024_nest_30.h5";
        NSIDE = 1024;
    }
    //const list<string> channels = list_of( "LFI28M" );
    vector<string> channels;
    channels.push_back("LFI27M");
    channels.push_back("LFI27S");
    channels.push_back("LFI28M");
    channels.push_back("LFI28S");

    dm = new H5PlanckDataManager(91, 94, channels, dataPath, pointingPath);
    if (DEBUG) {
        dm->setDatasetLength(1000);
    }
    dm->NSIDE = NSIDE;
    dm->NSTOKES = 3;
}

#endif

