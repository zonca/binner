#ifndef PLANCKDATAMANAGER_H
#define PLANCKDATAMANAGER_H

#include <string>
#include <list>
#include <vector>
using namespace std;


class PlanckDataManager

{
    private:
        long datasetLength_;
        int frequency_;
        string channel_;
        vector<string> fileNames;
        vector<string> channelNames;
        vector<string> pointingFileNames;
        vector<long> fileLengths;


    public:
        PlanckDataManager(int firstOd, int lastOd, const list<string> &channels, string dataPath, string pointingPath);
        int getData(string what, long iStart, int nElements, double* data);

        static string getLatestFilename(int od, int frequency, char type, string dataPath);

        double getDatasetLength(void) {
            return fileLengths.back();
        };
        void setDatasetLength(double length) {
            fileLengths.push_back(length);
        };
        int NSIDE;
        int NSTOKES;
        double getNPIX(void){
            return 12. * NSIDE * NSIDE +1; //total pixel size, each pixel is an element which contains 3 floats which are IQU
        }

};

#endif
