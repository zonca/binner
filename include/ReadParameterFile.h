#ifndef READPARAMETERFILE_H
#define READPARAMETERFILE_H

#include <boost/assign/list_of.hpp> // for 'list_of()'
#include "H5PlanckDataManager.h"

using namespace boost::assign; // bring 'list_of()' into scope

void readParameterFile(string parameterFilename, H5PlanckDataManager *& dm) {

    //TODO read parameter file

    bool DEBUG = true;
    int NSIDE;
    string dataPath, pointingPath;

    int firstOD, lastOD;

    dataPath = "/scratch/scratchdirs/zonca/pointing/lfi_ops_dx4_30.h5";

    NSIDE = 1024;
    firstOD = 91;
    lastOD = 141;

    if (DEBUG) {
        pointingPath = "/home/zonca/p/testdata/dx4_1024_nest_30_9293.h5";
        dataPath = "/home/zonca/p/testdata/lfi_ops_dx4.h5";
        firstOD = 92;
        lastOD = 92;
    } else {
        pointingPath = "/scratch/scratchdirs/zonca/pointing/dx4_1024_nest_30.h5";
    }
    vector<string> channels;
    channels.push_back("LFI27M");
    channels.push_back("LFI27S");
    channels.push_back("LFI28M");
    channels.push_back("LFI28S");

    dm = new H5PlanckDataManager(firstOD, lastOD, channels, dataPath, pointingPath);
    if (DEBUG) {
        dm->setDatasetLength(1000000);
    }
    dm->NSIDE = NSIDE;
    dm->NSTOKES = 3;
}

#endif

