#ifndef CSTR_H 
#define CSTR_H

#include "userData.h"

const int Ns_ = 20;
const double tpred = C_tpred;
const double ts = tpred/Ns_; 
const int lsim = 10;

const int s_nsteps = 20;    
const int m_nsteps = 10;

const char* saveFileStates = "cstrStates.txt";
const char* saveFileInputs = "cstrInputs.txt";
const char* saveFileTimes = "cstrInstants.txt";
const char* saveFileDimensions = "cstrDimensions.txt";
const char* saveFileKKT = "cstrKKT.txt";
const char* saveFileFeas = "cstrFeas.txt";
const char* saveFileSolveTime = "cstrSolveTime.txt";
const char* saveFileCG = "cstrCG.txt";

#endif