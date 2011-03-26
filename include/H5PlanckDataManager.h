#ifndef PLANCKDATAMANAGER_H
#define PLANCKDATAMANAGER_H

#include <string>
#include <list>
#include <vector>
#include <map>
using namespace std;
#include "H5Cpp.h"
#include "PointingStruct.h"
#include "Teuchos_Array.hpp"

typedef map<string, double> WeightDict;

class H5PlanckDataManager

{
    private:
        long TotalLength, LengthPerChannel, GlobalOffset;
        Teuchos::Array<string> Channels;
        string DataPath;
        string PointingPath;
        WeightDict Weights;

    public:
        H5PlanckDataManager(int firstOd, int lastOd, Teuchos::Array<string> channels, string dataPath, string pointingPath, WeightDict Weights);
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
        Teuchos::ArrayView<string> getChannels(void) {
            return Channels;
        };
        void setDatasetLength(double length) {
            TotalLength = length;
        };
        void setLengthPerChannel(double length) {
            LengthPerChannel = length;
        };
        int NSIDE;
        bool DEBUG;
        int NSTOKES;
        string outputFolder;
        double getNPIX(void){
            return 12. * NSIDE * NSIDE +1; //total pixel size, each pixel is an element which contains 3 floats which are IQU
        }

        int getIndexM(int row, int col) {
                return row * (2*NSTOKES-1 - row)/2 + col;
        }

};

#endif
