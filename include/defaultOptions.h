//
//  defaultOptions.h
//  MusAL
//
//  Created by Jean on 4/29/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//
#ifndef MusAL_defaultOptions_h
#define MusAL_defaultOptions_h

/**
** Default solver options
**/
/* Options for refinement phase: 
*	- REF_PCG_RE: PCG with restart when problem bound hit during search
*	- REF_MOR_TOR: Moré-Toraldo projections (projected searchs along PCG directions)
*/
#define REFIN_PCG_RES 1
#define REFIN_MOR_TOR 0

/* Options for Cauchy phase:
* - CCH_PRJ_SRCH: standard backtracking projected search
* - CCH_INTER_EXTRA: interpolation & extrapolation procedure
*/
#define CCH_PRJ_SRCH 1
#define CCH_INTER_EXTRA 0

/* Preconditioner */
#define NO_PRE 0 // No preconditioner 
#define BAND_PRE 1 // Band preconditioner 
#define JACO_PRE 0 // Diagonal preconditioner

/* Options for Jacobi preconditionner */
#define PREC_REG_COEF 1.

/* Default bandwidth for band preconditioner */
#define BAND_WIDTH 2

/* Pivot tolerance for LDL' on band preconditioner */
#define PIV_TOL 1E-5

/* Options for inexact Newton loop */
#define QUAD_CV 1 // For local quadratic convergence 
#define SUPL_CV 0 // For local superlinear convergence

/* Options for refinement search (Lin & Moré say that MAX_REF_ITER is equal to pb dimension) */
#define MAX_REF_ITER 200 // Max # refinement steps
#define MAX_PRJ_ITER 100 // Max # backtracking iter in proj search
#define REF_BETA_BCK 0.5 // Backtracking coef
#define REF_CST_DEC 0.0001 // Sufficient decrease constant 

/* Options for Cauchy search */
#define STEP_INI 1.0 // Initial step-size
#define BETA_BCK 0.7 // Backtracking coefficient
#define CST_BCK 0.0001 // Constant for sufficient decrease
#define MAX_BCK_ITER 100 // Max # backtracking iter
#define ACTIV_TOL 1E-20 // Tolerance on activity of bound constraints 

/* Options for interpolation-extrapolation */
#define BETA_INTER 0.1 
#define BETA_EXTRA 10.
#define MAX_INTER_ITER 50
#define MAX_EXTRA_ITER 50

/* Options for trust-region loop */
#define ETA_ZERO 1E-4 //0.3
#define ETA_ONE 0.25
#define ETA_TWO 0.75

#define SIGMA_ONE 0.25
#define SIGMA_TWO 0.5
#define SIGMA_THREE 4.0

#define MAX_PRIMAL_ITER 300

#define TR_RAD_INI 0.1
#define TR_RAD_TOL 1E-16

/* Options for SR1 updates */
#define SR_SHOOT_INI 1.
#define SR_PATH_INI 1. 

#define SR_SKIP 1E-8 // Constant for skipping condition in SR1 update

/* Options for dual loop */
#define RHO_INI 10.
#define MUL_PEN_COEF 10.

#define MAX_DUAL_ITER 1

#define KKT_TOL_INI 1.
#define EC_TOL_INI 1.

/* Default tolerance on inner KKT and 2-norm of equality constraints */
#define KKT_TOL 1E-4;
#define NEC_TOL 1E-5;

/* Default tolerance on infinity norm of momentum of primal sequence */
#define DIFF_TOL 1E-6  //1E-8

/* Default tolerance on magnitude of objective difference */
#define OBJ_TOL 1E-8    //1E-10

/* Default SPG options */
#define GAM_SPG 1E-4
#define HIST_LEN_SPG 10
#define SIG_ONE_SPG 0.1
#define SIG_TWO_SPG 0.9
#define AMIN_SPG 1E-30
#define AMAX_SPG 1E30
#define AINI_SPG 1.
#define MAX_BCK_IT_SPG 200

#endif