#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <cstring>
#include <list>

#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "H5Cpp.h"

#include "PlanckDataManager.h"
extern "C" {
#include "read_fits_lib.h"
}

using boost::format;
using namespace std;
using namespace H5;

PlanckDataManager::PlanckDataManager(int firstOd,int  lastOd, const list<string>& channels, string dataPath, string pointingPath)
{
    string fileName, pointingFileName;
    long prevLength = 0; 
    long length = 0;
    int frequency = 30; //FIXME get from channel tags

    for (int od=firstOd; od<=lastOd; od++) {

        BOOST_FOREACH ( string channel, channels )
        {
            fileName = getLatestFilename(od, frequency, 'R', dataPath);
            fileNames.push_back(fileName);

            channelNames.push_back(channel);

            pointingFileName =  getLatestFilename(od, frequency, 'P', pointingPath);
            pointingFileNames.push_back(pointingFileName);

            nelem(fileName.c_str(), &length);
            prevLength = prevLength + length;
            fileLengths.push_back(prevLength);
        }
    }

// DEBUG
//    for (unsigned int i=0; i<fileLengths.size();i++) {
//        cout << fileNames[i] << ':' << channelNames[i]  << ':' << fileLengths[i] << endl;
//    }
// /DEBUG
};

int PlanckDataManager::getPointing(long iStart, int nElements, int* pointing, double* qw, double* uw){
     H5File file( "/home/zonca/p/testdata/dx4_1024_nest_30_9293.h5", H5F_ACC_RDONLY );
     DataSet dataset = file.openDataSet( DATASET_NAME );

}

int PlanckDataManager::getData(string what, long iStart, int nElements, double* data){

    //double data[nElements];
    long firstelem, fileStart;
    int firstFilenElements;
    int fileIndex;
    double* d = data;
    bool reloop = true;

    fileIndex = int(lower_bound(fileLengths.begin(), fileLengths.end(), iStart) - fileLengths.begin());
    //cout << "Fileindex:" << fileIndex << endl;

    if (fileIndex == 0) {
        fileStart = 0;
    } else {
        fileStart = fileLengths[fileIndex - 1];
    }

    while (reloop)
    {
        if (iStart + nElements > fileLengths[fileIndex] ) {
            firstFilenElements = fileLengths[fileIndex] - iStart;
        } else {
            firstFilenElements = nElements;
            reloop = false;
        }

        firstelem = iStart - fileStart;
        //cout << "fileStart:" << fileStart << "  fileLengths[fileIndex] " <<  fileLengths[fileIndex] << endl;
        //cout << "firstelem:" << firstelem << " len " << firstFilenElements << endl;
        if (strcmp(what.c_str(), "data") == 0) {
            read_data(fileNames[fileIndex].c_str(), channelNames[fileIndex].c_str(), 1, firstelem, long(firstFilenElements), d);
        } else if  (strcmp(what.c_str(), "pointing") == 0) {
            read_data(pointingFileNames[fileIndex].c_str(), channelNames[fileIndex].c_str(), 1, firstelem, long(firstFilenElements), d);
        } else if  (strcmp(what.c_str(), "qw") == 0) {
            read_data(pointingFileNames[fileIndex].c_str(), channelNames[fileIndex].c_str(), 2, firstelem, long(firstFilenElements), d);
        } else if  (strcmp(what.c_str(), "uw") == 0) {
            read_data(pointingFileNames[fileIndex].c_str(), channelNames[fileIndex].c_str(), 3, firstelem, long(firstFilenElements), d);
        } else {
            cout << "ERROR type of data not recognized" << endl;
        }

        fileIndex++;
        fileStart = fileLengths[fileIndex - 1];
        iStart = fileStart;
        d += firstFilenElements;
        nElements -= firstFilenElements;
    }
    return 0;
}

string PlanckDataManager::getLatestFilename(int od, int frequency, char type, string dataPath){

    char filename[300];
    latest_exchange(od, frequency, type, dataPath.c_str(), filename);
    string latestFilename(filename);
    return latestFilename;
};

