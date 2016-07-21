#ifndef INV_PEND_H
#define INV_PEND_H

/***
** Parameters for inverted pendulum 
**/

#include "userData.h"

int Ns_,lsim,s_nsteps,m_nsteps;
double tpred,ts;

Ns_ = 30;
tpred = 1.5;
ts = tpred/Ns_; 
lsim = 80;

s_nsteps = 10;    
m_nsteps = 1;

const char* saveFileStates = "pendulumStates.txt";
const char* saveFileInputs = "pendulumInputs.txt";
const char* saveFileTimes = "pendulumInstants.txt";
const char* saveFileDimensions = "pendulumDimensions.txt";
const char* saveFileKKT = "pendulumKKT.txt";
const char* saveFileFeas = "pendulumFeas.txt";
const char* saveFileSolveTime = "pendulumSolveTime.txt";
const char* saveFileCG = "pendulumCG.txt";

/** State bounds */
/* Lower bounds on pendulum's states */
double xmin[DX];
xmin[0] = MIN_X_ONE;
xmin[1] = NEGINF;
xmin[2] = NEGINF; 
xmin[3] = NEGINF;

/* Upper bounds on pendulum's states */
double xmax[DX];
xmax[0] = MAX_X_ONE;
xmax[1] = POSINF;
xmax[2] = POSINF;
xmax[3] = POSINF;

/** Input bounds */
double umin[DU];

//    umin[0] = MIN_U_ONE;
//    umin[1] = MIN_U_TWO; // Unicycle

/* Upper bounds on crazyflie's input */
double umax[DU];

//    umax[0] = MAX_U_ONE;
//    umax[1] = MAX_U_TWO;

/* Initial state */
//xini[2] = 0;
    //xini[0] = 0.8; // Unicycle
    //xini[1] = -0.2;
    //xini[2] = 3.*M_PI/4.;

#endif