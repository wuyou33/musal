#ifndef CRAZY_H
#define CRAZY_H

/***
** Parameters for Walid's crazyflie 
**/

#include "userData.h"

const int Ns_ = 20;
const double tpred = 0.4; 
const double ts = 0.02;
const int lsim = 2;

const int s_nsteps = 10;
const int m_nsteps = 1;

const char* saveFileStates = "crazyStates.txt";
const char* saveFileInputs = "crazyInputs.txt";
const char* saveFileTimes = "crazyInstants.txt";
const char* saveFileDimensions = "crazyDimensions.txt";
const char* saveFileKKT = "crazyKKT.txt";
const char* saveFileFeas = "crazyFeas.txt";
const char* saveFileSolveTime = "crazySolveTime.txt";
const char* saveFileCG = "crazyCG.txt";

#endif