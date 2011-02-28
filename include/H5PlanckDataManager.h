#ifndef PLANCKDATAMANAGER_H
#define PLANCKDATAMANAGER_H

#include <string>
#include <list>
#include <vector>
using namespace std;
#include "H5Cpp.h"
#include "PointingStruct.h"

class H5PlanckDataManager

{
    private:
        long datasetLength_;
        int TotalLength, LengthPerChannel;
        vector<string> Channels;
        string DataPath;
        string PointingPath;
        int GlobalOffset;

    public:
        H5PlanckDataManager(int firstOd, int lastOd, vector<string> channels, string dataPath, string pointingPath);
        int getData(string what, long iStart, int nElements, double* data);
        int getPointing(long iStart, int nElements, pointing_t* pointing);

        double getDatasetLength(void) {
            return TotalLength;
        };
        void setDatasetLength(double length) {
            TotalLength = length;
        };
        int NSIDE;
        int NSTOKES;
        double getNPIX(void){
            return 12. * NSIDE * NSIDE +1; //total pixel size, each pixel is an element which contains 3 floats which are IQU
        }

};

#endif
