#ifndef READPARAMETERFILE_H
#define READPARAMETERFILE_H

#include "H5PlanckDataManager.h"

typedef std::map<string, double> WeightDict;

void readParameterFile(string parameterFilename, H5PlanckDataManager *& dm) {

    //TODO read parameter file

    WeightDict Weights;
    Weights["LFI18"] = 5.2814E+04;
    Weights["LFI19"] = 3.9294E+04;
    Weights["LFI23"] = 4.6195E+04;
    Weights["LFI22"] = 4.8167E+04;
    Weights["LFI20"] = 3.4468E+04;
    Weights["LFI21"] = 4.8501E+04;
    Weights["LFI24"] = 1.1531E+05;
    Weights["LFI25"] = 1.3288E+05;
    Weights["LFI26"] = 1.0654E+05;
    Weights["LFI27"] = 3.6656E+05;
    Weights["LFI28"] = 3.4432E+05;

    //TORSTI
    //Weights["LFI27"] = 363767.99;
    //Weights["LFI28"] = 342751.11;

    bool DEBUG = false;
    int NSIDE;
    string dataPath, pointingPath;

    int firstOD, lastOD;
    int freq = 44;

    dataPath = str( format("/scratch/scratchdirs/zonca/pointing/lfi_ops_dx6_%d.h5") % freq );

    NSIDE = 1024;
    firstOD = 91;
    lastOD = 563;

    if (DEBUG) {
        pointingPath = "/home/zonca/p/testdata/dx4_1_nest_30_9293.h5";
        dataPath = "/home/zonca/p/testdata/lfi_ops_dx4.h5";
        NSIDE = 1;
        firstOD = 93;
        lastOD = 93;
    } else {
        pointingPath = str( format("/scratch/scratchdirs/zonca/pointing/dx6_%d_horn_nest_%d.h5") % NSIDE % freq );
    }
    vector<string> channels;
    //channels.push_back("LFI24M");
    //channels.push_back("LFI24S");
    channels.push_back("LFI25M");
    channels.push_back("LFI25S");
    channels.push_back("LFI26M");
    channels.push_back("LFI26S");


    dm = new H5PlanckDataManager(firstOD, lastOD, channels, dataPath, pointingPath, Weights);
    if (DEBUG) {
        dm->setLengthPerChannel(30);
    }
    dm->NSIDE = NSIDE;
    dm->NSTOKES = 5;
    dm->DEBUG = DEBUG;
}

#endif
