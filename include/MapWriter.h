#ifndef MAPWRITER_H
#define MAPWRITER_H

#include <string>
#include <vector>

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_DataAccess.h"

extern "C" {
    #include "read_fits_lib.h"
}

using namespace std;


class MapWriter

{
    private:
        int NPIX;

        Epetra_BlockMap * TargetMap;
        Epetra_BlockMap fromMap_;
        Epetra_Export * Exporter;
        Epetra_MpiComm Comm_;

    public:
        MapWriter(const Epetra_BlockMap & fromMap, const Epetra_MpiComm & Comm, int NPIXin)  : fromMap_(fromMap), Comm_(Comm) { 
            
            NPIX = NPIXin;

            int NumMyElements_target;

            if( Comm.MyPID() == 0 ) {
                NumMyElements_target = NPIX;
            } else {
                NumMyElements_target = 0;
            }

            TargetMap = new Epetra_BlockMap(NPIX,NumMyElements_target,3,0,Comm);
            Exporter = new Epetra_Export(fromMap,*TargetMap);
        }

        int write(Epetra_Vector vector, string fileName) {

            Epetra_Vector exp(*TargetMap);

            exp.Export(vector,*Exporter,Add);

            if( Comm_.MyPID() == 0 ) {
                cout << "Exported to proc 0" << endl;
                double * hitmap_vector;
                hitmap_vector = new double[NPIX*3];
                exp.ExtractCopy(hitmap_vector);
                int mapNPIX = NPIX -1;
                double * intensity;
                double * q;
                double * u;
                intensity = new double[mapNPIX];
                q = new double[mapNPIX];
                u = new double[mapNPIX];
                int i = 0;
                for (i = 0; i < mapNPIX;i=i+1){
                    intensity[i] = hitmap_vector[3*i];
                    q[i] = hitmap_vector[3*i +1];
                    u[i] = hitmap_vector[3*i +2];
                }
                cout << "Writing Map" << endl;
                write_map(fileName.c_str(), intensity, q, u, mapNPIX);
                cout << "Map written" << endl;
            } 
            return 0;
        }


};

#endif
