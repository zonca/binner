#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <cstring>
#include <list>

#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "DataManager.h"

using boost::format;
using namespace std;

DataManager::DataManager(int firstOd,int  lastOd,  Teuchos::Array<string> channels, string dataPath, string pointingPath)
{
    string fileName, pointingFileName;
    long prevLength = 0; 
    long length = 0;

    for (int od=firstOd; od<=lastOd; od++) {

        BOOST_FOREACH ( string channel, channels )
        {
            format fName("/%02d/%s.fits");
            fName % od % channel; 
            fileName = dataPath + fName.str(); 
            fileNames.push_back(fileName);

            fName = format("/%02d/%s_pix.fits");
            fName % od % channel; 
            pointingFileName = pointingPath + fName.str();
            pointingFileNames.push_back(pointingFileName);

            nelem(fileName.c_str(), &length);
            prevLength = prevLength + length;
            fileLengths.push_back(prevLength);
        }
    }

// DEBUG
    for (unsigned int i=0; i<fileLengths.size();i++) {
        cout << fileNames[i] << ':' << fileLengths[i] << endl;
    }
// /DEBUG
};

int DataManager::adjustDistribution(int iStart, int nElements) {

    int FirstDatasetIndex = int(upper_bound(fileLengths.begin(), fileLengths.end(), iStart) - fileLengths.begin());
    int LastDatasetIndex = int(upper_bound(fileLengths.begin(), fileLengths.end(), iStart+nElements) - fileLengths.begin());
    int AdjustedNumMyElements = nElements;
    int LostElements = 0;
    int LastFilenElements, fileIndex, NewLen;
    int fileStart=0;

    if (FirstDatasetIndex == 0) {
        fileStart = 0;
    } else {
        fileStart = fileLengths[FirstDatasetIndex - 1];
    }
    //left boundary - every element is LOST to the prev Proc
    LostElements = (iStart-fileStart) % BaselineLength;
    int NewIstart = ceil((float)(iStart-fileStart)/BaselineLength) * BaselineLength;
    //cout << "FirstDatasetIndex:" << FirstDatasetIndex << " iStart:" << iStart << "fileStart:" << fileStart <<" nElements:" << nElements << " NewIstart:" << NewIstart << endl;
    if (NewIstart > 0) {
        AdjustedNumMyElements -= NewIstart - (iStart-fileStart);
    }

    if (LastDatasetIndex == 0) {
        fileStart = 0;
    } else {
        fileStart = fileLengths[LastDatasetIndex - 1];
    }
    //right boundary - even 1 element more gets full baseline
    if (iStart + nElements < fileLengths[LastDatasetIndex] ) {
        LastFilenElements = iStart + nElements - fileStart;
        NewLen = ceil((float)LastFilenElements/BaselineLength) * BaselineLength;
        NewLen = min(NewLen, int(fileLengths[LastDatasetIndex] - fileStart));
        AdjustedNumMyElements += NewLen - LastFilenElements;
    //cout << "iStart:" << iStart << "fileStart:" << fileStart <<" LastFilenElements:" << LastFilenElements << " NewLen:" << NewLen << endl;
    }

    return AdjustedNumMyElements;
}

int DataManager::getData(string what, long iStart, int nElements, double* data){

    //double data[nElements];
    long firstelem, fileStart;
    int firstFilenElements;
    int fileIndex;
    double* d = data;
    bool reloop = true;

    fileIndex = int(lower_bound(fileLengths.begin(), fileLengths.end(), iStart) - fileLengths.begin());
    cout << "Fileindex:" << fileIndex << endl;

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

        cout << "fileName:" << fileNames[fileIndex] << endl;
        cout << "fileStart:" << fileStart << "  fileLengths[fileIndex] " <<  fileLengths[fileIndex] << endl;
        cout << "firstelem:" << firstelem << " len " << firstFilenElements << endl;
        if (strcmp(what.c_str(), "data") == 0) {
            read_data(fileNames[fileIndex].c_str(), "DATA", 2, firstelem, long(firstFilenElements), d);
        } else if  (strcmp(what.c_str(), "pointing") == 0) {
            read_data(pointingFileNames[fileIndex].c_str(), "PIX", 2, firstelem, long(firstFilenElements), d);
        } else if  (strcmp(what.c_str(), "qw") == 0) {
            read_data(pointingFileNames[fileIndex].c_str(), "PIX", 3, firstelem, long(firstFilenElements), d);
        } else if  (strcmp(what.c_str(), "uw") == 0) {
            read_data(pointingFileNames[fileIndex].c_str(), "PIX", 4, firstelem, long(firstFilenElements), d);
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
