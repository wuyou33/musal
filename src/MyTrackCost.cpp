//
//  MyTrackCost.cpp
//  MusAL
//
//  Created by Jean on 3/13/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//
#include <algorithm>
#include <iostream>
#include <string.h>

#include "macros.h"
#include "MyTrackCost.h"

/**
* Default constructor 
*/
MyTrackCost::MyTrackCost(){

    memset(Q,0,dy*SZDBL);
    memset(R,0,du*SZDBL);
    memset(gX,0,dx*SZDBL);
    memset(gU,0,du*SZDBL);

    /* Bioreactor, no quadratic term */


    /* Weights on states for Walid's crazyflie */
/*    Q[0] = 200.;
    Q[1] = 200.;
    Q[2] = 100.;
    Q[3] = 1.;
    Q[4] = 1.;
    Q[5] = 1.;
    Q[6] = 20.;
    Q[7] = 20.;
    Q[8] = 20.;
    Q[9] = 1.;
    Q[10] = 1.;
    Q[11] = 1.; */

    /* Weights on input for Walid's crazyflie */
 /*   R[0] = 5.E-8;
    R[1] = 5.E-3;
    R[2] = 5.E-3;
    R[3] = 1.E-3;
    */

    /* Weights on states for inverted pendulum and crane*/
/*    Q[0] = 10.; // Horizontal position
    Q[1] = 0.1;  
    Q[2] = 10.; // Angle 
    Q[3] = 0.1; */
    /* Weights  on input for inverted pendulum and crane */
//    R[0] = 0.01;

    /* Weights on states for pendulum with polynomial dynamics */
/*    Q[0] = 10.;
    Q[1] = 0.1;
    Q[2] = 10.;
    Q[3] = 10.;
    Q[4] = 0.1; */

    /* Weights on inputs for pendulum with polynomial dynamics */
//    R[0] = 0.01;

    /* Weights on states for dc motor*/
//    Q[0] = 0.;
//    Q[1] = 100.;

    /* Weights on input for dc motor */
//    R[0] = 0.1;

    /* Weights on states for unicycle */
//    Q[0] = 500.;
//    Q[1] = 500.;
//    Q[2] = 1000.;

    /* Weights on inputs for unicycle */
//    R[0] = 3.;
//    R[1] = 3.;
}

/** 
* Destructor 
*/
MyTrackCost::~MyTrackCost(){}

/** 
* Compute stage-cost 
*/
void MyTrackCost::evalStage(double* val,double* x,double* u,double* yref,double* uref){

    /* USER-DEFINED */
    double val_,df;

    /**
    ** CSTR
    **/
//    *val = -x[1]*u[0]/C_tpred;

    /**
    ** Bioreactor
    **/
    *val = -1.*B_D*x[2]/B_T;

    /**
    ** For crazyflie 
    **/
    /* Output tracking term */
/*    df = *x-*yref;
    val_ = *Q*df*df*0.5; 
    df = *(x+1)-*(yref+1);
    val_ += *(Q+1)*df*df*0.5;
    df = *(x+2)-*(yref+2);
    val_ += *(Q+2)*df*df*0.5;
    df = *(x+3)-*(yref+3);
    val_ += *(Q+3)*df*df*0.5;
    df = *(x+4)-*(yref+4);
    val_ += *(Q+4)*df*df*0.5;
    df = *(x+5)-*(yref+5);
    val_ += *(Q+5)*df*df*0.5;
    df = *(x+6)-*(yref+6);
    val_ += *(Q+6)*df*df*0.5;
    df = *(x+7)-*(yref+7);
    val_ += *(Q+7)*df*df*0.5;
    df = *(x+8)-*(yref+8);
    val_ += *(Q+8)*df*df*0.5;
    df = *(x+9)-*(yref+9);
    val_ += *(Q+9)*df*df*0.5;
    df = *(x+10)-*(yref+10);
    val_ += *(Q+10)*df*df*0.5;
    df = *(x+11)-*(yref+11);
    val_ += *(Q+11)*df*df*0.5; */

    /* Input term */
/*    df = *u-*uref;
    val_ += *R*df*df*0.5;
    df = *(u+1)-*(uref+1);
    val_ += *(R+1)*df*df*0.5;
    df = *(u+2)-*(uref+2);
    val_ += *(R+2)*df*df*0.5;
    df = *(u+3)-*(uref+3);
    val_ += *(R+3)*df*df*0.5;

    *val = val_; */

    /* Inverted pendulum and crane */
/*    df = *x-*yref;
    val_ = *Q*df*df; 
    df = *(x+1)-*(yref+1);
    val_ += *(Q+1)*df*df;
    df = *(x+2)-*(yref+2);
    val_ += *(Q+2)*df*df;
    df = *(x+3)-*(yref+3);
    val_ += *(Q+3)*df*df; */ 

    /* Inverted pendulum with polynomial dynamics */
/*    df = *x-*yref;
    val_ = *Q*df*df; 
    df = *(x+1)-*(yref+1);
    val_ += *(Q+1)*df*df;
    df = *(x+2)-*(yref+2);
    val_ += *(Q+2)*df*df;
    df = *(x+3)-*(yref+3);
    val_ += *(Q+3)*df*df;
    df = *(x+4)-*(yref+4);
    val_ += *(Q+4)*df*df; */

//    df = *u-*uref;
//    val_ += *R*df*df;


    /* DC motor */
/*    df = *(x+1)-*yref;
    val_ = *(Q+1)*df*df;

    df = *u-*uref;
    val_ += *R*df*df;

    *val = val_; */

    /* Unicycle */
/*    df = *x-*yref;
    val_ = *Q*df*df;
    df = *(x+1)-*(yref+1);
    val_ += *(Q+1)*df*df;
    df = *(x+2)-*(yref+2);
    val_ += *(Q+2)*df*df;

    df = *u-*uref;
    val_ += *R*df*df;
    df = *(u+1)-*(uref+1);
    val_ += *(R+1)*df*df; 
    *val = val_;*/
}


/** 
* Compute stage-cost gradient wrt state 
*/
void MyTrackCost::evalGradXstage(double* x,double* u,double* yref,double* uref){
    
    /* USER-DEFINED */
    double df;

    /**
    ** CSTR
    **/
//    gX[1] = -C_b*u[0]/C_tpred;

    /**
    ** Bioreactor
    **/
    gX[2] = -1.*B_D/B_T;

    /**
    ** For crazyflie
    **/
/*    df = *x-*yref;
    *gX = *Q*df;
    df = *(x+1)-*(yref+1);
    *(gX+1) = *(Q+1)*df;
    df = *(x+2)-*(yref+2);
    *(gX+2) = *(Q+2)*df;
    df = *(x+3)-*(yref+3);
    *(gX+3) = *(Q+3)*df;
    df = *(x+4)-*(yref+4);
    *(gX+4) = *(Q+4)*df;
    df = *(x+5)-*(yref+5);
    *(gX+5) = *(Q+5)*df;
    df = *(x+6)-*(yref+6);
    *(gX+6) = *(Q+6)*df;
    df = *(x+7)-*(yref+7);
    *(gX+7) = *(Q+7)*df;
    df = *(x+8)-*(yref+8);
    *(gX+8) = *(Q+8)*df;
    df = *(x+9)-*(yref+9);
    *(gX+9) = *(Q+9)*df;
    df = *(x+10)-*(yref+10);
    *(gX+10) = *(Q+10)*df;
    df = *(x+11)-*(yref+11);
    *(gX+11) = *(Q+11)*df; */

    /* Inverted pendulum and crane */
/*    df = *x-*yref;
    *gX = *Q*df*2.;
    df = *(x+1)-*(yref+1);
    *(gX+1) = *(Q+1)*df*2.;
    df = *(x+2)-*(yref+2);
    *(gX+2) = *(Q+2)*df*2.;
    df = *(x+3)-*(yref+3);
    *(gX+3) = *(Q+3)*df*2.; */

    /* Inverted pendulum with polynomical dynamics */
/*    df = *x-*yref;
    *gX = *Q*df*2.;  
    df = *(x+1)-*(yref+1);
    *(gX+1) = *(Q+1)*df*2.;
    df = *(x+2)-*(yref+2);
    *(gX+2) = *(Q+2)*df*2.;
    df = *(x+3)-*(yref+3);
    *(gX+3) = *(Q+3)*df*2.;
    df = *(x+4)-*(yref+4);
    *(gX+4) = *(Q+4)*df*2.; */

    /* DC motor */
 /*   *gX = 0.;
    df = *(x+1)-*yref;
    *(gX+1) = *(Q+1)*df*2.; */

    /*Uunicycle */
/*    df = *x-*yref;
    *gX = *Q*df*2.;
    df = *(x+1)-*(yref+1);
    *(gX+1) = *(Q+1)*df*2.;
    df = *(x+2)-*(yref+2);
    *(gX+2) = *(Q+2)*df*2.; */
}


/** 
* Compute stage-cost gradient wrt input 
*/
void MyTrackCost::evalGradUstage(double* x,double* u,double* yref,double* uref){
    
    /* USER-DEFINED */
    double df;

    /**
    ** CSTR
    **/
//    gU[0] = -C_b*x[1]/C_tpred;

    /**
    ** Bioreactor
    **/

    /**
    ** For crazyflie 
    **/
/*    df = *u-*uref;
    *gU = *R*df;
    df = *(u+1)-*(uref+1);
    *(gU+1) = *(R+1)*df;
    df = *(u+2)-*(uref+2);
    *(gU+2) = *(R+2)*df;
    df = *(u+3)-*(uref+3);
    *(gU+3) = *(R+3)*df; */

/*    df = *u-*uref;
    *gU = *R*df*2.; */
//    df = *(u+1)-*(uref+1);
//    *(gU+1) = *(R+1)*df*2.;

}
