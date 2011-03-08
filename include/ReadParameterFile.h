#ifndef READPARAMETERFILE_H
#define READPARAMETERFILE_H

#include "H5PlanckDataManager.h"

void readParameterFile(string parameterFilename, H5PlanckDataManager *& dm) {

    //TODO read parameter file

    bool DEBUG = false;
    int NSIDE;
    string dataPath, pointingPath;

    int firstOD, lastOD;

    dataPath = "/scratch/scratchdirs/zonca/pointing/lfi_ops_dx4_30.h5";

    NSIDE = 1024;
    firstOD = 92;
    lastOD = 92;

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

