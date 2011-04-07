#ifndef READPARAMETERFILE_H
#define READPARAMETERFILE_H

#include "DataManager.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"

typedef std::map<string, double> WeightDict;

void readParameterFile(string parameterFilename, DataManager *& dm) {

    Teuchos::RCP<Teuchos::ParameterList> Config=Teuchos::getParametersFromXmlFile(parameterFilename);

    //Teuchos::Array<string> channels = Teuchos::getParameter< Teuchos::Array<string> >(*Config,"Channels");
    //Teuchos::Array<string> channels = Config->get<Teuchos::Array<string> >("Channels");
    //Teuchos::Array<string> channels =  Teuchos::tuple<string>("ch1q","ch1q");
    //channels = Config->get("Channels",channels);

    Teuchos::Array<string> channels =  Teuchos::tuple<string>("ch1q","ch1u");
    //Teuchos::Array<string> channels =  Teuchos::tuple<string>("data");

    int NSIDE = Config->get("Nside", 1024);
    string pointingPath = str( format(Config->get("PointingPath", "")) );
    string dataPath = str( format(Config->get("DataPath", "")) );

    dm = new DataManager(Config->get("FirstOD", 1), Config->get("LastOD", 1), channels, dataPath, pointingPath);

    dm->DEBUG = Config->get("Debug", false);
    dm->BaselineLength = Config->get("BaselineLength", 2000);
    dm->NSIDE=NSIDE;

    dm->outputFolder = Config->get("OutputFolder", "out");

    ////Teuchos::writeParameterListToXmlFile(Config, parameterFilename);
    dm->NSTOKES = 2;
    //cout << *Config << endl;
}

#endif
