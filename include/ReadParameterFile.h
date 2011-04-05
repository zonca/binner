#ifndef READPARAMETERFILE_H
#define READPARAMETERFILE_H

#include "H5PlanckDataManager.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"

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

    Teuchos::RCP<Teuchos::ParameterList> Config=Teuchos::getParametersFromXmlFile(parameterFilename);

    int freq = Config->get("Frequency", 70);

    Teuchos::Array<string> channels = Config->template get<Teuchos::Array<string> >("Channels");

    int NSIDE = Config->get("Nside", 1024);
    string pointingPath = str( format(Config->get("PointingPath", "/scratch/scratchdirs/zonca/pointing/%d/dx6_%d_horn_nest_%d")) % freq % NSIDE % freq );
    string dataPath = str( format(Config->get("DataPath", "/scratch/scratchdirs/zonca/pointing/%d/lfi_ops_dx6_%d")) % freq % freq );

    dm = new H5PlanckDataManager(Config->get("FirstOD", 91), Config->get("LastOD", 563), channels, dataPath, pointingPath, Weights);

    dm->DEBUG = Config->get("Debug", false);
    dm->NSIDE=NSIDE;

    dm->outputFolder = Config->get("OutputFolder", "/scratch/scratchdirs/zonca/out");

    if (dm->DEBUG) {
        dm->setLengthPerChannel(30);
    }
    ////Teuchos::writeParameterListToXmlFile(Config, parameterFilename);
    dm->NSTOKES = 3;
    if (Config->get("Spurious", false)) {
        dm->NSTOKES += channels.size()/2.;
    }

    //cout << *Config << endl;
}

#endif
