//
//  SolverOptions.cpp
//  MusAL
//
//  Created by Jean on 4/23/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//
#include <iostream>


#include "SolverOptions.h"
#include "defaultOptions.h"


SolverOptions::SolverOptions(){
    /* Options for refinement phase */
    maxRefIter = MAX_REF_ITER;
    maxItSrch = MAX_PRJ_ITER;
    betaSrch = REF_BETA_BCK;
    cstSrch = REF_CST_DEC;

    /* Options for Cauchy phase */
    sini = STEP_INI;
    betaBck = BETA_BCK;
    cstBck = CST_BCK;
    maxItBck = MAX_BCK_ITER;
    
    betaInter = BETA_INTER;
    betaExtra = BETA_EXTRA;
    maxInterIter = MAX_INTER_ITER;
    maxExtraIter = MAX_EXTRA_ITER;

    /* Options for trust-region iterations */
    trRadIni = TR_RAD_INI;
    
    eta0 = ETA_ZERO;
    eta1 = ETA_ONE;
    eta2 = ETA_TWO;
    
    sig1 = SIGMA_ONE;
    sig2 = SIGMA_TWO;
    sig3 = SIGMA_THREE;
    
    maxPit = MAX_PRIMAL_ITER;
    
    /* Options for SR1 hessian approximations */
    sr_shoot_ini = SR_SHOOT_INI;
    sr_path_ini = SR_PATH_INI;

    sr_skip = SR_SKIP;

    /* Options for regularization coefficient in diagonal preconditionner */
    pr_reg_coef = PREC_REG_COEF;

    /* Bandwidth of band preconditioner */
    bwd = BAND_WIDTH;

    /* Pivot tolerance for LDL' on band preconditioner */
    pivTol = PIV_TOL;

    /* Options on (dual) outer loop */
    maxDit = MAX_DUAL_ITER;
    
    mulPen = MUL_PEN_COEF;
    rhoIni = RHO_INI;

    kktTolIni = KKT_TOL_INI;
    ecTolIni = EC_TOL_INI;

    /* Default tolerances */
    activTol = ACTIV_TOL;

    trRadTol = TR_RAD_TOL;

    kktTolAbs = KKT_TOL;
    necTolAbs = NEC_TOL;

    diffTol = DIFF_TOL;
    objTol = OBJ_TOL;

    /* SPG options */
    gamSpg = GAM_SPG;

    histLenSpg = HIST_LEN_SPG;

    sigOneSpg = SIG_ONE_SPG;
    sigTwoSpg = SIG_TWO_SPG;

    aminSpg = AMIN_SPG;
    amaxSpg = AMAX_SPG;
    ainiSpg = AINI_SPG;
  
    maxBckItSpg = MAX_BCK_IT_SPG;
}


SolverOptions::~SolverOptions(){
    
}
