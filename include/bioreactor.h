#ifndef BIOREC_H 
#define BIOREC_H

#include "userData.h"

const int Ns_ = 20;
const double tpred = 48.;
const double ts = tpred/Ns_; 
const int lsim = 200;

const int s_nsteps = 20;    
const int m_nsteps = 1;

/* Initial state of periodic trajectory (economic optimum) */
const double xini_0 = 5.479699008;
const double xini_1 = 17.18171851;
const double xini_2 = 19.28396712;
const double xini_3 = 0.;
const double xini_4 = 0.;

/* Files for storing states, inputs,... */
const char* saveFileStates = "bioreactorStates.txt";
const char* saveFileInputs = "bioreactorInputs.txt";
const char* saveFileTimes = "bioreactorInstants.txt";
const char* saveFileDimensions = "bioreactorDimensions.txt";
const char* saveFileKKT = "bioreactorKKT.txt";
const char* saveFileFeas = "bioreactorFeas.txt";
const char* saveFileSolveTime = "bioreactorSolveTime.txt";
const char* saveFileCG = "bioreactorCG.txt";

#endif