#ifndef PLANCKDATAMANAGER_H
#define PLANCKDATAMANAGER_H

#include <string>
#include <list>
#include <vector>
#include <map>
using namespace std;
#include "H5Cpp.h"
#include "PointingStruct.h"

typedef map<string, double> WeightDict;

class H5PlanckDataManager

{
    private:
        long TotalLength, LengthPerChannel, GlobalOffset;
        vector<string> Channels;
        string DataPath;
        string PointingPath;
        WeightDict Weights;

    public:
        H5PlanckDataManager(int firstOd, int lastOd, vector<string> channels, string dataPath, string pointingPath, WeightDict Weights);
        int getData(string channel, long iStart, int nElements, double* data);
        int getPointing(string channel, long iStart, int nElements, int * pix, double * qw, double * uw);
        double getWeight(string channel) {
            return Weights[channel.substr(0, 5)];
        }

        long getDatasetLength(void) {
            return TotalLength;
        };
        long getLengthPerChannel(void) {
            return LengthPerChannel;
        };
        vector<string> getChannels(void) {
            return Channels;
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
