#ifndef CRANE_H
#define CRANE_H

/***
** Crane parameters
**/

#include "userData.h"

/*    ts = 0.3; 
    tsim = ts;
    tpred = 3.;

    s_hfix = ts/10.;    
    m_hfix = ts/6.; */

    /** State bounds */
    /* Lower bounds on crazyflie's states */
    double xmin[DX];

//    xmin[0] = MIN_X_ONE; // For inverted pendulum
//    xmin[0] = NEGINF; // For crane
//    xmin[1] = NEGINF;
//    xmin[2] = NEGINF; //-1.;
//    xmin[3] = NEGINF; //-1.;
//    xmin[4] = NEGINF;
//    xmin[0] = MIN_X_ONE; // DC motor
//    xmin[1] = MIN_X_TWO;
//    xmin[0] = MIN_X_ONE; // Unicycle
//    xmin[1] = MIN_X_TWO;
//    xmin[2] = MIN_X_THREE;

    /* Upper bounds on crazyflie's states */
    double xmax[DX];

//    xmax[0] = MAX_X_ONE; // For inverted pendulum
//    xmax[0] = POSINF; // For crane
//    xmax[1] = POSINF;
//    xmax[2] = POSINF; //1.
//    xmax[3] = POSINF; //1.
//    xmax[4] = POSINF;
//    xmax[0] = MAX_X_ONE; // DC motor
//    xmax[1] = MAX_X_TWO;
//    xmax[0] = MAX_X_ONE; // Unicycle
//    xmax[1] = MAX_X_TWO;
//    xmax[2] = MAX_X_THREE;


#endif