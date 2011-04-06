#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <string>
#include <list>
#include <vector>
#include "Teuchos_Array.hpp"
using namespace std;

extern "C" {
#include "read_fits_lib.h"
}

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
        DataManager(int firstOd, int lastOd, Teuchos::Array<string> channels, string dataPath, string pointingPath);
        int getData(string what, long iStart, int nElements, double* data);

        static string getLatestFilename(int od, int frequency, char type, string dataPath);

        double getLength(void) {
            return fileLengths.back();
        };

        int adjustDistribution(int MinMyGID, int NumMyElements);
        int NSIDE;
        bool DEBUG;
        int NSTOKES, BaselineLength;
        string outputFolder;
        double getNPIX(void){
            return 12. * NSIDE * NSIDE +1; //total pixel size, each pixel is an element which contains 3 floats which are IQU
        }

        int getIndexM(int row, int col) {
                return row * (2*NSTOKES-1 - row)/2 + col;
        }
};

#endif
