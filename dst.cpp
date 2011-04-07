#include <string>
#include <math.h>

#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>

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
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Time.h"
#include <EpetraExt_MatrixMatrix.h>
#include "Teuchos_CommandLineProcessor.hpp"

#include "DataManager.h"
#include "CreateMatrices.h"
#include "ReadParameterFile.h"

#include "AztecOO.h"

#include "H5Cpp.h"

using namespace std;
using namespace H5;
using boost::format;

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif  

Teuchos::CommandLineProcessor My_CLP;
string parameterFilename = "config.xml";
My_CLP.setOption("p", &parameterFilename, "Parameter filename");
My_CLP.parse( argc, argv );

int MyPID = Comm.MyPID();
int i_M, a, s_index, err;
double ScalarMultiplier = 1.;

int MAX_SAMPLES_PER_PROC = 5e6;

Epetra_Time time(Comm);
DataManager* dm;

log(MyPID, parameterFilename);
readParameterFile(parameterFilename, dm);

bool hasI = dm->NSTOKES>2;
log(MyPID, format("Number of elements [mil]: %.10d") % (dm->getLength()/1.e6));

Epetra_Map LinearMap(dm->getLength(), 0, Comm);

int NumMyElements = dm->adjustDistribution(LinearMap.MinMyGID(), LinearMap.NumMyElements());

//cout << "ID:" << MyPID << "Before:"<< LinearMap.NumMyElements()<< "My elements:" << NumMyElements <<endl;
Epetra_Map Map(dm->getLength(), NumMyElements, 0, Comm);

Epetra_Map PixMap(dm->getNPIX(),0,Comm);
Epetra_BlockMap PixBlockMap(dm->getNPIX(),dm->NSTOKES,0,Comm);

// declaring matrices
log(MyPID,"POINTING MATRIX");
Epetra_MultiVector summap(PixMap, dm->NSTOKES);
int Msize =  dm->NSTOKES * (dm->NSTOKES + 1) / 2;
Epetra_MultiVector M(PixMap, Msize);
Epetra_Vector tempmap(PixMap);
Epetra_Vector * hitmap = new Epetra_Vector(PixMap);
Epetra_CrsMatrix * P;
Epetra_CrsGraph * Graph; 

Epetra_MultiVector yqu = Epetra_MultiVector(Map, 3);

double ** yqu_view = new double *[3];
yqu.ExtractView(&yqu_view);

string LABEL[3] = {"I", "Q", "U"};

if (!hasI)  {
    LABEL[0]="Q"; LABEL[1]="U";
}

Epetra_Vector pix = Epetra_Vector(Map);
//Epetra_IntVector pix = Epetra_IntVector(Map);

double * pix_view;
pix.ExtractView(&pix_view);

Epetra_Vector tempvec = Epetra_Vector(Map);
Epetra_Vector tempvec2 = Epetra_Vector(Map);

double weight = 1.;
int ch_number = 0;

log(MyPID,"READ POINTING");
time.ResetStartTime();
dm->getData("pointing", Map.MinMyGID(), NumMyElements, pix_view);
dm->getData("qw", Map.MinMyGID(), NumMyElements, yqu_view[1]);
dm->getData("uw", Map.MinMyGID(), NumMyElements, yqu_view[2]);
log(MyPID, format("Read data timer %f") % time.ElapsedTime());

log(MyPID,"READ DATA");
time.ResetStartTime();
dm->getData("data", Map.MinMyGID(), NumMyElements, yqu_view[0]);
log(MyPID, format("Read data timer %f") % time.ElapsedTime());

time.ResetStartTime();
createGraph(Map, PixMap, pix, Graph);
log(MyPID, format("Create Graph timer %f") % time.ElapsedTime());

time.ResetStartTime();
P = new Epetra_CrsMatrix(Copy, *Graph);
P->PutScalar(1.);
P->FillComplete(PixMap, Map);
log(MyPID, format("Create P timer %f") % time.ElapsedTime());

log(MyPID,"SUM MAP");
time.ResetStartTime();

////// I
if (hasI) {
    log(MyPID,"I");
    P->Multiply1(true,*(yqu(0)),tempmap); //SUMMAP = Pt y
    summap(0)->Update(weight, tempmap, 1.);
    log(MyPID, format("%f") % time.ElapsedTime());
}

time.ResetStartTime();
log(MyPID,"HitMap");
tempvec.PutScalar(1.);
P->Multiply1(true,tempvec,tempmap);
hitmap->Update(1., tempmap, 1.);
log(MyPID, format("%f") % time.ElapsedTime());
time.ResetStartTime();

WriteH5Vec(hitmap, "hitmap", dm->outputFolder);
delete hitmap;

// Q U
log(MyPID,"QU");
int qu_i=0;
for (int i=1; i<3; ++i) { // Q=1 U=2
    tempvec.Multiply(1., *(yqu(0)), *(yqu(i)), 0.);
    P->Multiply1(true,tempvec,tempmap); //SUMMAP = Pt y
    if (hasI) {
        qu_i = i;
    } else {
        qu_i = i-1;
    }
    summap(qu_i)->Update(weight, tempmap, 1.);
}
log(MyPID, format("%f") % time.ElapsedTime());

log(MyPID, "Create M");
log(MyPID, "IQU loop");
time.ResetStartTime();
if (hasI) {
    for (int j=0; j<3; ++j) {
        for (int k=j; k<3; ++k) {
            log(MyPID,format("M %d %d") % j % k);
            if (k == 0) { //also j=0
                tempvec.PutScalar(1.);
            } else if (j == 0 ) {
                tempvec.Update(1., *(yqu(k)), 0.);
            } else {
                tempvec.Multiply(1., *(yqu(j)), *(yqu(k)), 0.);
            }
            P->Multiply1(true,tempvec,tempmap); //SUMMAP = Pt y
            i_M = dm->getIndexM(j, k);
            log(MyPID, format("Setting M %d") % i_M );
            M(i_M)->Update(weight, tempmap, 1.);
        }
    } // M loop IQU
} else {
    for (int j=0; j<2; ++j) {
        for (int k=j; k<2; ++k) {
            log(MyPID,format("M %d %d") % j % k);
            tempvec.Multiply(1., *(yqu(j+1)), *(yqu(k+1)), 0.);
            i_M = dm->getIndexM(j, k);
            P->Multiply1(true,tempvec,tempmap); //SUMMAP = Pt y
            log(MyPID, format("Setting M %d") % i_M );
            M(i_M)->Update(weight, tempmap, 1.);
        }
    } // M loop IQU
}
log(MyPID, format("IQU loop: %f") % time.ElapsedTime());

log(MyPID,"Writing SUMMAP");

time.ResetStartTime();
for (int j=0; j<dm->NSTOKES; ++j) {
    WriteH5Vec(summap(j), "summap_" + LABEL[j], dm->outputFolder);
}
log(MyPID, format("%f") % time.ElapsedTime());


log(MyPID,"Computing RCOND and Inverting");
time.ResetStartTime();
Epetra_Vector rcond(PixMap);

Epetra_Vector binmap(PixMap);
invertM(PixMap, dm, M, rcond, summap);
log(MyPID, format("%f") % time.ElapsedTime());

WriteH5Vec(&rcond, "rcondmap", dm->outputFolder);

log(MyPID,"Writing BINMAP");
time.ResetStartTime();
for (int j=0; j<dm->NSTOKES; ++j) {
    WriteH5Vec(summap(j), "binmap_" + LABEL[j], dm->outputFolder);
}
log(MyPID, format("%f") % time.ElapsedTime());

log(MyPID,"Writing M");
time.ResetStartTime();
for (int j=0; j<dm->NSTOKES; ++j) {
    for (int k=j; k<dm->NSTOKES; ++k) {
        i_M = dm->getIndexM(j, k);
        WriteH5Vec(M(i_M), "M_" + LABEL[j] + "_" + LABEL[k], dm->outputFolder);
    }
}
log(MyPID, format("%f") % time.ElapsedTime());

//Destriping
log(MyPID,"M_time");
Epetra_MultiVector M_time(Map, Msize);
P->Multiply(false, M, M_time);
log(MyPID, format("%f") % time.ElapsedTime());

log(MyPID,"Z");
Epetra_CrsMatrix Z(Copy, Map, 70);
EpetraExt::MatrixMatrix::Multiply(*P, false, *P, true, Z);
log(MyPID, format("%f") % time.ElapsedTime());
Epetra_CrsMatrix Z2(Z);

log(MyPID,"Z weighting");
time.ResetStartTime();
//if (hasI) {
//    cout << "NOTIMPLEMENTED" << endl;
//} else {
//    for (int i=0; i<NumMyElements;i++) {
//        NumMyEntries = Z.NumMyEntries(i);
//        Z.ExtractMyRowView(i, NumEntries, Values, Indices)
//        for (int c=0; c<NumEntries; c++) {
//            GlobalCol = Z.GCID(Indices[c]);
//            Values[c] += qw
//
//    }
//
//    tempvec2.PutScalar(0.);
//    Z.PutScalar(0.);
//    for (int j=0; j<2; ++j) {
//        for (int k=j; k<2; ++k) {
//            log(MyPID,format("M %d %d") % j % k);
//            if (k != j) { //cross terms must be doubled
//                ScalarMultiplier = 2.;
//            } else {
//                ScalarMultiplier = 1.;
//            }
//            tempvec.Multiply(ScalarMultiplier, *(yqu(j+1)), *(yqu(k+1)), 0.);
//            i_M = dm->getIndexM(j, k);
//            tempvec2.Multiply(1., tempvec, *(M_time(i_M)), 1.);
//            log(MyPID, format("Setting M %d") % i_M );
//        }
//    }
//    Z.LeftScale(tempvec2); //by row
//}
Z.LeftScale(*(yqu(1))); //by row
tempvec.Multiply(1., *(yqu(1)), *(M_time(0)), 0.);
tempvec.Multiply(1., *(yqu(2)), *(M_time(1)), 1.);
Z.RightScale(tempvec); //by col
Z2.LeftScale(*(yqu(2))); //by row
tempvec.Multiply(1., *(yqu(1)), *(M_time(1)), 0.);
tempvec.Multiply(1., *(yqu(2)), *(M_time(2)), 1.);
Z2.RightScale(tempvec); //by col
log(MyPID, format("%f") % time.ElapsedTime());

int NumEntries;
for (int i=0; i<NumMyElements;i++) {
    NumEntries = Z.NumMyEntries(i);
    for (int c=0; c<NumEntries; c++) {
        Z[i][c] += Z2[i][c];
    }
}

log(Comm.MyPID(),"1-PB");
Z.Scale(-1);
Epetra_Vector diagZ(Map);
Z.ExtractDiagonalCopy(diagZ);
for(unsigned int i=0; i<Map.NumMyElements(); ++i) {
        diagZ[i] = 1. + diagZ[i];
}
Z.ReplaceDiagonalValues(diagZ);

//cout << Z << endl;

//Create F
log(Comm.MyPID(),"----------------Creating F n_base * time");
time.ResetStartTime();
int NumLocalBaselines=0;
vector<int> BaselineLengths;
dm->numLocalBaselines(Map.MinMyGID(), Map.NumMyElements(), NumLocalBaselines, BaselineLengths);
Epetra_Map BaselinesMap(-1,NumLocalBaselines, 0, Comm);

Epetra_CrsMatrix F(Copy, Map,1, true);
int * MyGlobalElements = Map.MyGlobalElements();
int * MyBaselineGlobalElements = BaselinesMap.MyGlobalElements();
double one = 1.0;
int k=0;
int baselinenum;
for(unsigned int i=0 ; i<BaselinesMap.NumMyElements(); ++i ) { //loop on local baselines
    baselinenum = MyBaselineGlobalElements[i];
    //cout <<BaselineLengths[i]  << "   " << baselinenum << endl;
    for (unsigned int j=0; j<BaselineLengths[i]; j++) {
        //cout << "k " << k << "   " << baselinenum << "   " << MyGlobalElements[k] << endl;
        F.InsertGlobalValues(MyGlobalElements[k], 1, &one, &baselinenum);
        k++;
    }
}
F.FillComplete(BaselinesMap, Map);

log(MyPID, format("%f") % time.ElapsedTime());
log(Comm.MyPID(),"----------------Creating FZT n_base * time");
time.ResetStartTime();
Epetra_CrsMatrix FZT(Copy, Map,10);

log(Comm.MyPID(),"M-M");
err = EpetraExt::MatrixMatrix::Multiply(Z, true, F, false, FZT); 
log(Comm.MyPID(),"M-M END");
if (err != 0) {
    cout << "err "<<err<<" from MatrixMatrix::Multiply"<<endl;
    return(err);
}

log(MyPID, format("%f") % time.ElapsedTime());
log(Comm.MyPID(),"----------------Creating RHS n_base");

time.ResetStartTime();
Epetra_Vector RHS(BaselinesMap);
FZT.Multiply(true,*(yqu(0)),RHS);

log(MyPID, format("%f") % time.ElapsedTime());
log(Comm.MyPID(),"----------------Creating D n_base x n_base");

time.ResetStartTime();
Epetra_CrsMatrix D(Copy, BaselinesMap, 30);

log(Comm.MyPID(),"M-M");
err = EpetraExt::MatrixMatrix::Multiply(FZT, true, F, false, D); 
log(Comm.MyPID(),"M-M END");
if (err != 0) {
    cout << "err "<<err<<" from MatrixMatrix::Multiply"<<endl;
    return(err);
}

log(MyPID, format("%f") % time.ElapsedTime());

//// Solving
Epetra_Vector baselines(BaselinesMap);

time.ResetStartTime();
// Create Linear Problem
Epetra_LinearProblem problem(&D, &baselines, &RHS);
// Create AztecOO instance
AztecOO solver(problem);

int options[AZ_OPTIONS_SIZE];
double params[AZ_PARAMS_SIZE];
AZ_defaults(options, params);
solver.SetAllAztecOptions( options );
solver.Iterate(100, 1.0E-10);

cout << "Solver performed " << solver.NumIters() << " iterations." << endl
   << "Norm of true residual = " << solver.TrueResidual() << endl;

log(MyPID, format("%f") % time.ElapsedTime());

Epetra_Vector destripedTOD(Map);
F.Multiply(false,baselines,destripedTOD);

for(unsigned int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local elements
    destripedTOD[i] = yqu[0][i] - destripedTOD[i];
}

//cout << destripedTOD << endl;

log(MyPID,"Writing DESTRIPED");
time.ResetStartTime();
log(MyPID, format("%f") % time.ElapsedTime());
for (int j=0; j<dm->NSTOKES; ++j) {
    tempvec.Multiply(1., *(yqu(1)), *(M_time(j)), 0.);
    tempvec.Multiply(1., *(yqu(2)), *(M_time(j+1)), 1.);
    for(unsigned int i=0 ; i<Map.NumMyElements(); ++i ) { //loop on local elements
        tempvec[i] *= destripedTOD[i];
    }
    P->Multiply(true, tempvec, tempmap);
    //cout << tempmap << endl;
    WriteH5Vec(&tempmap, "destripedmap_" + LABEL[j], dm->outputFolder);
}

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
return(0);

};
