#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <string>
#include <list>
#include <vector>
using namespace std;


class DataManager

{
    private:
        long datasetLength_;
        int frequency_;
        string channel_;
        vector<string> fileNames;
        vector<string> pointingFileNames;
        vector<long> fileLengths;


    public:
        DataManager(int firstOd, int lastOd, const list<string> &channels, string dataPath, string pointingPath);
        int getData(string what, long iStart, int nElements, double* data);

        static string getLatestFilename(int od, int frequency, char type, string dataPath);

        double getDatasetLength(void) {
            return fileLengths.back();
        };

};

#endif
