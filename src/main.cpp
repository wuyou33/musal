//
//  main.cpp
//  MusAL
//
//  Created by Jean on 3/11/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//
#include <iostream>
#include <iomanip>
#include <string.h>
#include <math.h>
#include <time.h>
// ADOL-C includes (all headers in "adolc/adolc.h")
//#include <adolc/adolc.h>

#include "userData.h"
#include "MyDynamics.h"
#include "MyTrackCost.h"
#include "MyPathConstraint.h"
#include "MyMayer.h"
#include "Simulator.h"
#include "SolverOptions.h"
#include "ExplicitButcherTabs.h"

/* Physical set-up */
//#include "crazyflie.h"
#include "bioreactor.h"
//#include "cstr.h"

int main(int argc, const char * argv[])
{
    int i;
    ExplicitButcherTabs butch;
    
    /** 
    ** Set simulator
    **/
    Simulator<MyDynamics,MyTrackCost,MyPathConstraint,MyMayer> sim;
    SolverOptions opts;

    /** State bounds */
    /* Lower bounds */
    double xmin[DX];
    xmin[0] = NEGINF;
    xmin[1] = NEGINF;
    xmin[2] = NEGINF;
    xmin[3] = NEGINF; 
    xmin[4] = NEGINF;

    /* Upper bounds */
    double xmax[DX];
    xmax[0] = POSINF;
    xmax[1] = POSINF;
    xmax[2] = POSINF;
    xmax[3] = POSINF;
    xmax[4] = POSINF; 

    sim.setStateBounds(xmin,xmax);
    for (i=0;i<DX;++i)   
        std::cout<<"xmin="<<xmin[i]<<", xmax="<<xmax[i]<<std::endl;

    /** Terminal state bounds */
    double xmaxT[DX];
    xmaxT[0] = POSINF;
    xmaxT[1] = POSINF;
    xmaxT[2] = POSINF;
    xmaxT[3] = MAX_XT_4;
    xmaxT[4] = MAX_XT_5;

    double xminT[DX];
    xminT[0] = xmin[0];
    xminT[1] = xmin[1];
    xminT[2] = xmin[2];
    xminT[3] = xmin[3];
    xminT[4] = xmin[4];

    sim.setTerminalStateBounds(xminT,xmaxT);
    for (i=0;i<DX;++i)   
        std::cout<<"xminT="<<xminT[i]<<", xmaxT="<<xmaxT[i]<<std::endl; 

    /** Input bounds */
    double umin[DU];
    umin[0] = MIN_U_1;
//    umin[1] = MIN_U_2;
    double umax[DU];
    umax[0] = MAX_U_1;
//    umax[1] = MAX_U_2;

    sim.setInputBounds(umin,umax);
    for (i=0;i<DU;++i)   
        std::cout<<"umin="<<umin[i]<<", umax="<<umax[i]<<std::endl;

    sim.setTolKKTal(0.08);
    sim.setSimulationLength(lsim);
    sim.setPredictionTime(tpred);
    sim.setSamplingPeriod(ts);
    sim.setNumShootingNodes(Ns_+1);
    
    sim.setNumIntStepShooting(m_nsteps);
    sim.setButcherShooting(butch.gets_erk_4(),
                           butch.getA_erk_4(),
                           butch.getb_erk_4(),
                           butch.getc_erk_4());
//    sim.setButcherShooting(butch.gets_erk_trap(),butch.getA_erk_trap(),butch.getb_erk_trap(),butch.getc_erk_trap());
//    sim.setButcherShooting(butch.gets_erk_rung(),butch.getA_erk_rung(),butch.getb_erk_rung(),butch.getc_erk_rung());

    sim.setNumIntStepSimulation(s_nsteps);
    sim.setButcherSimulation(butch.gets_erk_4(),
                             butch.getA_erk_4(),
                             butch.getb_erk_4(),
                             butch.getc_erk_4());
    
    /* Set options for augmented Lagrangian solver */
    opts.setMaxPrimalIters(500);
    opts.setBacktrackMaxIter(50);
    opts.setMaxDualIters(1);
    opts.setBacktrackCoef(0.5);
    opts.setBacktrackConst(0.0001);
    opts.setBetaInter(0.1);
    opts.setBetaExtra(5.);
    opts.setInitialPenalty(50.);
    opts.setPenaltyMulCoef(5.);
    opts.setBandWidth(10);
    opts.setPivotTolerance(1.E-6);
//    opts.setEcTolIni(1.E-5);
//    opts.setKktTolIni(1.E-5);
    sim.setSolverOptions(opts);
    sim.init();

    /* Run */
    /** Initial state */
    double xini[DX];
    /* Bioreactor */
    xini[0] = xini_0;
    xini[1] = xini_1;
    xini[2] = xini_2;
    xini[3] = xini_3;
    xini[4] = xini_4;
    /* CSTR */
/*    xini[0] = 2.;
    xini[1] = 1.08;
    xini[2] = 108.;
    xini[3] = 107.7; */
    sim.setInitialState(xini);
    sim.run();
    return 0;
    /* Save */
    sim.saveDimensions(saveFileDimensions);
    sim.saveTimes(saveFileTimes);
    sim.saveInputs(saveFileInputs);
    sim.saveStates(saveFileStates);
    sim.saveKKT(saveFileKKT);
    sim.saveFeas(saveFileFeas);
    sim.saveSolveTime(saveFileSolveTime);
    sim.saveCG(saveFileCG);

    return 0;
}


