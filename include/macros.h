//
//  macros.h
//  MusAL
//
//  Created by Jean on 3/13/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_macros_h
#define MusAL_macros_h

/* Simulation mode: real-time or full */
#define SIM_RTI 1
#define SIM_FUL 0 

/* Options for basic linear algebra */
#define BLAS 0 // Use CBLAS
#define MANU 1 // Use manuel implementation

/* Options for scaling variables */
#define SCAL_VAR 1
#define UNSCAL_VAR 0

/* Error detection */
#define DETEC_CAUCH_ERR 0
#define DETEC_REFIN_ERR 0

/* Display options */
/* Printing flag for TPCG loop */
#define PRINT_CG 0
/* Printing flag for trust-region loop */	
#define PRINT_TR 1
/* Printing flag for Moré-Toraldo refinement */
#define PRINT_MOR_TOR 0
/* Printing flag for dual loop */
#define PRINT_DUA 0
/* Counting cumulative CG iters */
#define COUNT_CG 1

/* Sizes */
#define SZDBL sizeof(double)
#define SZINT sizeof(int)

/* Macros */
#define MIN(a,b) ((a<b) ? a : b)
#define MAX(a,b) ((a>b) ? a : b)
#define ABS(a) ((a>=0) ? a : -(a))
#define PROJ(x,l,u) ((x<=l) ? l : (x>=u) ? u : x)

/* Positive & negative infinite */
#define POSINF 1.0E20
#define NEGINF -1.0E20

/* Positive and negative zero */
#define POSZER 1.0E-20
#define NEGZER -1.0E-20

/**********************************************************/
/*************** DEFAULT PARAMETERS FOR INTEGRATOR *******/
/**********************************************************/

/* Maximum grid length */
#define LGMAX 100000

/* Max. # integration steps */
#define NSMAX 1E5;

/**
** PARAMETERS FOR MyTrajectory
**/
/* Maximum reference length */
#define  LREFMAX 1000

#endif
