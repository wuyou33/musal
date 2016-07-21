//
//  userDimensions.h
//  MusAL
//
//  Created by Jean on 4/29/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_userDimensions_h
#define MusAL_userDimensions_h

#include <math.h>

#include "recFuncs.h"

/**
** CSTR 
**/
/*#define DX 4
#define DU 2
#define DY 4
#define DPC 0
#define DPER 0
#define MAY 0 */
/* Parameters */
/*#define C_hab 4.2
#define C_hbc -11.
#define C_had -41.85
#define C_k10 1.287E12
#define C_k20 9.043E6
#define C_e1 9758.3
#define C_e2 8560.
#define C_t0 273.15
#define C_tin 104.9
#define C_cin 5.1E3
#define C_a 30.828
#define C_b 86.688
#define C_d 3.522E-4
#define C_g 0.1
#define C_tpred 3. */
/* State bounds */
/*#define MIN_X_1 0.
#define MAX_X_1 6.
#define MIN_X_2 0.
#define MAX_X_2 4.
#define MIN_X_3 70.
#define MAX_X_3 150.
#define MIN_X_4 0.
#define MAX_X_4 300. */
/* Input bounds */
/*#define MIN_U_1 3.
#define MAX_U_1 35.
#define MIN_U_2 0.
#define MAX_U_2 200. */

/**
** Bioreactor
**/
/* Macro for indicating if OCP with periodic constraints 
   is to be solved (from scratch) to compute economically optimal trajectory */
#define PEROCP 1
/* Macro indicating if bioreactor is used */
#define BIOREC 1
/* State, input and output dimensions */
#define DX 5
#define DU 1
#define DY 5
/* Path-constraint dimension */
#define DPC 0
/* Mayer term */
#define MAY 0
/* Periodicity constraints dimension 
	(First DPER components of initial & final shooting nodes) */
#define DPER 3
/* System's parameters */
#define B_D 0.15 // Dilution rate
#define B_T 48. // End time of cycle
#define B_mm 0.48 // 
#define B_Km 1.2
#define B_Pm 50.
#define B_Yxs 0.4
#define B_Ki 22.
#define B_a 2.2
#define B_b 0.2
/* Terminal state bound constraints */
#define MAX_XT_4 32.9
#define MAX_XT_5 5.8 
/* Input bound constraints */
#define MIN_U_1 28.7
#define MAX_U_1 40.0 

/**
** Walid's crazyflie
**/
/* State, input and output dimensions */
/*#define DX 15
#define DU 4
#define DY 12 */
/* Path-constraint dimension */
//#define DPC 0
/* Mayer term */
//#define MAY 0
/* Gravity */
//#define Q_g -9.8065
/* Mass and length */
/*#define Q_m 0.02946
#define Q_l 0.046 */
/* 6th parameter estimate in Walid's report (inertia, drag & lift coefficients) */
/*#define Q_Jx 0.000021757 
#define Q_Jy 0.000021757 
#define Q_Jz 0.000062614
#define Q_Cf 0.22026 
#define Q_Cm 0.012614 */
/* Inner controller parameters */
/*#define Q_KpPhi 3.5
#define Q_KiPhi 2.0
#define Q_KpThe 3.5 
#define Q_KiThe 2.0
#define Q_KpDps 70.0 
#define Q_KiDps 50.0
#define Q_KpDph 70.0 
#define Q_KpDth 70.0 */
/* Constants */
/*const double CPI4 = Cosine<60>(M_PI/4);
#define RAD2DEG 180./M_PI
#define CPWM 1./65535. */
/* State bound constraints */
/*#define MIN_X_1 -10.0
#define MAX_X_1 10.0
#define MIN_X_2 -10.0
#define MAX_X_2 10.0
#define MIN_X_3 -10.0
#define MAX_X_3 10.0
#define MIN_X_4 -10.0 
#define MAX_X_4 10.0
#define MIN_X_5 -10.0
#define MAX_X_5 10.0
#define MIN_X_6 -10.0
#define MAX_X_6 10.0
#define MIN_X_7 -M_PI/6 
#define MAX_X_7 M_PI/6
#define MIN_X_8 -M_PI/6
#define MAX_X_8 M_PI/6
#define MIN_X_9 -M_PI
#define MAX_X_9 M_PI
#define MIN_X_10 -2.*M_PI
#define MAX_X_10 2.*M_PI
#define MIN_X_11 -2.*M_PI
#define MAX_X_11 2.*M_PI
#define MIN_X_12 -2.*M_PI
#define MAX_X_12 2.*M_PI
#define MIN_X_13 -20.0
#define MAX_X_13 20.0
#define MIN_X_14 -20.0
#define MAX_X_14 20.0
#define MIN_X_15 -500.0
#define MAX_X_15 500.0 */
/* Input bound constraints */
/*#define MIN_U_1 0.0
#define MAX_U_1 65535.0
#define MIN_U_2 -30.0
#define MAX_U_2 30.0
#define MIN_U_3 -30.0
#define MAX_U_3 30.0
#define MIN_U_4 -200.0
#define MAX_U_4 200.0*/

/**
** Inverted pendulum
**/
/* State, input and output dimensions */
/*#define DX 4
#define DU 1
#define DY 4 */
/* Path-constraint dimension */
//#define DPC 0
/* Mayer term */
//#define MAY 0
/* System parameters */
/*#define PRM_l 0.5
#define PRM_g 9.81
#define PRM_m 0.1
#define PRM_M 1.0 */
/* Bound constraints on state */
/*#define MIN_X_ONE -2.0
#define MAX_X_ONE 2.0 */
/* Bound constraints on input */
/*#define MIN_U_ONE -20.0
#define MAX_U_ONE 20.0 */

/*****************************************************************************/

/**
* Inverted pendulum with polynomial dynamics (need to introduce path constraint c^2+s^2=1)
*/
/* State, input and output dimensions */
/*#define DX 5
#define DU 1
#define DY 5 */

/* Path-constraint dimension */
//#define DPC 1

/* System parameters */
/*#define PRM_l 0.5
#define PRM_g 9.81
#define PRM_m 0.1
#define PRM_M 1.0 */

/* Bound constraints on horizontal position */
/*#define MIN_X_ONE -2.0
#define MAX_X_ONE 2.0 */

/* Bound constraints on input */
/*define MIN_U_ONE -20.0
#define MAX_U_ONE 20.0 */

/*****************************************************************************/
/**
* Crane model
*/
/* State, input and output dimensions */
//#define DX 4
//#define DU 1
//#define DY 4

/* Path-constraint dimension */
//#define DPC 0

/* System parameters */
//#define PRM_g 9.81
//#define PRM_b 0.2

/* Bound constraints on state */

/* Bound constraints on input */
//#define MIN_U_ONE -1.0
//#define MAX_U_ONE 1.0

/*****************************************************************************/
/**
* DC motor
*/
/* State, input and output dimensions */
//#define DX 2
//#define DU 1
//#define DY 1

/* Path-constraint dimension */
//#define DPC 0

/* System parameters */
//#define PRM_Ra 12.548 
//#define PRM_La 0.307
//#define PRM_km 0.22567
//#define PRM_ua 60.
//#define PRM_B 0.00783
//#define PRM_J 0.00385
//#define PRM_tl 1.47

/* Bound constraints on state */
//#define MIN_X_ONE -10.0 
//#define MAX_X_ONE 10.0
//#define MIN_X_TWO -10.0
//#define MAX_X_TWO 10.0

/* Bound constraints on input */
//#define MIN_U_ONE -5.0
//#define MAX_U_ONE 5.0

/*****************************************************************************/

/** 
* Unicycle
*/
/* State, input and output dimensions */
//#define DX 3
//#define DU 2
//#define DY 3

/* Path-constraint dimension */
//#define DPC 0

/* Bound constraints on states */
//#define MIN_X_ONE -100. 
//#define MAX_X_ONE 100.
//#define MIN_X_TWO -100.
//#define MAX_X_TWO 100.
//#define MIN_X_THREE -100.*M_PI 
//#define MAX_X_THREE 100.*M_PI

/* Bound constraints on inputs */
//#define MIN_U_ONE 0.
//#define MAX_U_ONE 0.5
//#define MIN_U_TWO -0.5*M_PI
//#define MAX_U_TWO 0.5*M_PI


#endif
