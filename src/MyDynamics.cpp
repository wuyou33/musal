//
//  MyDynamics.cpp
//  MusAL
//
//  Created by Jean on 3/13/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//
#include <string.h>
#include <algorithm>
#include <math.h>
#include <iostream>

#include "MyDynamics.h"

/**
 * Default constructor 
 */
MyDynamics::MyDynamics(){
    
    int dxdx,dxdu;
    
    dxdx = dx*dx;
    dxdu = dx*du;
    
    memset(JacX,0,dxdx*SZDBL);
    memset(JacXt,0,dxdx*SZDBL);
    memset(JacU,0,dxdu*SZDBL);
    memset(JacUt,0,dxdu*SZDBL);
}

/** 
 * Destructor 
 */
MyDynamics::~MyDynamics(){}

/**
 * Method to evaluate ode rhs 
 */
void MyDynamics::rhs(double* xdot,double* x,double* u){

    /* USER-DEFINED */
    double *xxd,*xx,*uu;

    xx = x;
    uu = u;
    xxd = xdot;

    /**
    ** CSTR reactor
    **/
/*    double x1,x2,x3,x4;
    double u1,u2;

    x1 = *xx++;
    x2 = *xx++;
    x3 = *xx++;
    x4 = *xx++;
    u1 = *uu++;
    u2 = *uu++;

    double k1,k2;
    double x12;

    k1 = C_k10*exp(-C_e1/(x3+C_t0));
    k2 = C_k20*exp(-C_e2/(x3+C_t0));
    x12 = x1*x1;

    *xxd++ = -k1*x1-k2*x12+(C_cin-x1)*u1;
    *xxd++ = k1*(x1-x2)-x2*u1;
    *xxd++ = -C_d*(k1*(C_hab*x1+C_hbc*x2)+k2*x12*C_had)+C_a*(x4-x3)+(C_tin-x3)*u1;
    *xxd++ = C_b*(x3-x4)+C_g*u2;
*/
    /**
    ** Bioreactor
    **/
    double x1,x2,x3,x4,x5;
    double u1;

    x1 = *xx++; // X
    x2 = *xx++; // S
    x3 = *xx++; // P
    x4 = *xx++; // Sf^
    x5 = *xx++; // X^
    u1 = *uu++; // Sf

    /* Intermediate vars */
    double mu;
    mu = B_mm*(1.-x3/B_Pm)*x2/(B_Km+x2+x2*x2/B_Ki);

    *xxd++ = (-B_D+mu)*x1;
    *xxd++ = B_D*(u1-x2)-mu*x1/B_Yxs;
    *xxd++ = -B_D*x3+(B_a*mu+B_b)*x1;
    *xxd++ = u1/B_T;
    *xxd++ = x1/B_T;

    /**
    ** Walid's quadcopter: 15 states (12 quad states, 3 pid states), 4 inputs 
    **/
/*    
    double rx,ry,rz;
    double vx,vy,vz;
    double phi,the,psi;
    double wPhi,wThe,wPsi;
    double xc1,xc2,xc3;
    double cPsi,sPsi;
    double cPhi,sPhi;
    double cThe,sThe,tThe;
    double scThePhi;
    double u1,u2,u3,u4;
    double pw1,pw2,pw3,pw4;
    double w1,w2,w3,w4;
    double w1s,w2s,w3s,w4s;
    double uphi,uthe,udps;
    double ephi,ethe,edps;
    double duwPhi,duwThe;
 */
    /* Extract states */
    /* Position */
/*    rx = *xx;
    ry = *xx++;
    rz = *xx++;*/
    /* Speed */
/*    vx = *xx++;
    vy = *xx++;
    vz = *xx++; */
    /* Angular position */
/*    phi = *xx++;
    the = *xx++;
    psi = *xx++; */
    /* Angular speeds */
/*    wPhi = *xx++;
    wThe = *xx++;
    wPsi = *xx++; */
    /* PID states */
/*    xc1 = *xx++;
    xc2 = *xx++;
    xc3 = *xx++; */
    /* Some precomputations */
/*    cPsi = cos(psi);
    sPsi = sin(psi);
    cThe = cos(the);
    sThe = sin(the);
    tThe = sThe/cThe;
    cPhi = cos(phi);
    sPhi = sin(phi);
    scThePhi = sThe*cPhi; */
    /* Extract inputs */
//    u1 = *uu++; // Total thrust
//    u2 = *uu++; // Rolling moment
//    u3 = *uu++; // Pitching moment
//    u4 = *uu++; // Yawing moment
    /* Compute PWMs */
/*    ephi = u2-RAD2DEG*phi;
    ethe = -u3+RAD2DEG*the;
    edps = -u4-RAD2DEG*wPsi;
    uphi = Q_KpPhi*ephi+Q_KiPhi*xc1;
    uthe = Q_KpThe*ethe+Q_KiThe*xc2;
    udps = Q_KpDps*edps+Q_KiDps*xc3;
    duwPhi = Q_KpDph*(uphi-RAD2DEG*wPhi);
    duwThe = Q_KpDth*(uthe+RAD2DEG*wThe); 
    w1 = CPWM*(u1-0.5*(duwPhi-duwThe)-udps);
    w2 = CPWM*(u1-0.5*(duwPhi+duwThe)+udps);
    w3 = CPWM*(u1+0.5*(duwPhi-duwThe)-udps);  
    w4 = CPWM*(u1+0.5*(duwPhi+duwThe)+udps);
    w1s = w1*w1; w2s = w2*w2; w3s = w3*w3; w4s = w4*w4;
    pw1 = Q_Cf*(w1s+w2s+w3s+w4s);
    pw2 = Q_Cf*Q_l*CPI4*(-w1s-w2s+w3s+w4s);
    pw3 = Q_Cf*Q_l*CPI4*(-w1s+w2s+w3s-w4s);
    pw4 = Q_Cm*(-w1s+w2s-w3s+w4s); */
    /* Quad dynamics */
/*    *xxd++ = vx;
    *xxd++ = vy;
    *xxd++ = vz;
    *xxd++ = (cPsi*scThePhi+sPsi*sPhi)*u1/Q_m;
    *xxd++ = (sPsi*scThePhi-cPsi*sPhi)*u1/Q_m;
    *xxd++ = Q_g+cThe*cPhi*u1/Q_m;
    *xxd++ = wPhi+tThe*sPhi*wThe+cPhi*tThe*wPsi;
    *xxd++ = cPhi*wThe-sPhi*wPsi;
    *xxd++ = (sPhi*wThe+cPhi*wPsi)/cThe;
    *xxd++ = (pw2-(Q_Jz-Q_Jy)*wThe*wPsi)/Q_Jx;
    *xxd++ = (pw3-(Q_Jx-Q_Jz)*wPhi*wPsi)/Q_Jy; 
    *xxd++ = (pw4-(Q_Jy-Q_Jx)*wPhi*wThe)/Q_Jz; */
    /* PID dynamics */
/*    *xxd++ = ephi;
    *xxd++ = ethe;
    *xxd++ = edps; */

    /* Inverted pendulum
    (TO DO: turn into polynomial dynamics to make evaluation of dynamics & jacobians potentially cheaper)
     */
/*    double x_1,x_2,x_3,x_4,x_4_2;
    double u_1;
    double cx_3,sx_3,csx_3;
    double den;
    
    x_1 = *x; // Horizontal pos
    x_2 = *(x+1); // Horizontal speed
    x_3 = *(x+2);  // Angle
    x_4 = *(x+3);   // Angular speed
    x_4_2 = x_4*x_4;
    
    cx_3 = cos(x_3);
    sx_3 = sin(x_3);
    csx_3 = cx_3*sx_3;

    u_1 = *u;
    
    den = PRM_M+PRM_m-PRM_m*cx_3*cx_3;
    
    *xdot = x_2;
    *(xdot+1) = (PRM_m*PRM_l*sx_3*x_4_2+PRM_m*PRM_g*csx_3+u_1)/den;
    *(xdot+2) = x_4;
    *(xdot+3) = -(PRM_m*PRM_l*csx_3*x_4_2+u_1*cx_3+(PRM_M+PRM_m)*PRM_g*sx_3)/(PRM_l*den); */

    /* 
    * Inverted pendulum with polynomial dynamics
    */
/*    double x_1,x_2,x_3,x_4,x_5;
    double x_4_5s;
    double u_1;
    double den;

    x_1 = *x; // Horiz. pos.
    x_2 = *(x+1); // Horiz. spd.
    x_3 = *(x+2); // Angle cosine
    x_4 = *(x+3); // Angle sine
    x_5 = *(x+4); // Angular speed

    u_1 = *u;

    den = PRM_M+PRM_m-PRM_m*x_3*x_3;
    x_4_5s = x_4*x_5*x_5;

    *xdot = x_2;
    *(xdot+1) = (PRM_m*PRM_l*x_4_5s+PRM_m*PRM_g*x_3*x_4+u_1)/den; 
    *(xdot+2) = -x_4*x_5;
    *(xdot+3) = x_3*x_5;
    *(xdot+4) = -(PRM_m*PRM_l*x_3*x_4_5s+u_1*x_3+(PRM_M+PRM_m)*PRM_g*x_4)/(PRM_l*den); */

    /* Crane */
/*    double x_1,x_2,x_3,x_4;
    double u_1;

    x_1 = *x; // p: horizontal trolley position
    x_2 = *(x+1); // v: trolley velocity
    x_3 = *(x+2); // phi: excitation angle
    x_4 = *(x+3); // om: angular velocity of mass point

    u_1 = *u;
    
    *xdot = x_2;
    *(xdot+1) = u_1;
    *(xdot+2) = x_4;
    *(xdot+3) = -1.*PRM_g*sin(x_3)-u_1*cos(x_3)-PRM_b*x_4; */

    /* DC motor */
/*    double x_1,x_2;
    double u_1;

    x_1 = *x;
    x_2 = *(x+1);

    u_1 = *u;

    *xdot = (PRM_ua-PRM_Ra*x_1-PRM_km*x_2*u_1)/PRM_La;
    *(xdot+1) = (PRM_km*x_1*u_1-PRM_B*x_2-PRM_tl)/PRM_J; */

    /* Unicycle */
/*    double x_3,u_1,u_2;

    x_3 = *(x+2);

    u_1 = *u;
    u_2 = *(u+1);

    *xdot = u_1*cos(x_3);
    *(xdot+1) = u_1*sin(x_3);
    *(xdot+2) = u_2; */
}

/** 
* Display dynamics rhs
*/
void MyDynamics::displayRhs(double* x,double* u){

    int i;
    double xdot[dx];

    rhs(xdot,x,u);

    std::cout<<"ODE rhs="<<std::endl;
    for (i=0;i<dx;++i)
        std::cout<<xdot[i]<<std::endl;
}

/** 
 * Method to evaluate state-jacobian of ode rhs 
 */
void MyDynamics::jxRhs(double* x,double* u){

    /* USER-DEFINED */
    double *xx,*uu;
    xx = x;
    uu = u;

    /**
    ** CSTR
    **/
/*    double x1,x2,x3,x4;
    double u1,u2;

    x1 = *xx++;
    x2 = *xx++;
    x3 = *xx++;
    x4 = *xx++;
    u1 = *uu++;
    u2 = *uu++;

    double den,den2;
    double k1,k2;
    double dk1,dk2;
    double x12;

    den = x3+C_t0;
    den2 = den*den;
    k1 = C_k10*exp(-C_e1/den);
    k2 = C_k20*exp(-C_e2/den);
    dk1 = k1*C_e1/den2;
    dk2 = k2*C_e2/den2;
    x12 = x1*x1;

    JacX[0] = -k1-2.*k2*x1-u1;
    JacX[2] = -dk1*x1-dk2*x12;

    JacX[4] = k1;
    JacX[5] = -k1-u1;
    JacX[6] = dk1*(x1-x2);

    JacX[8] = -C_d*(k1*C_hab+2.*k2*C_had*x1);
    JacX[9] = -C_d*k1*C_hbc;
    JacX[10] = -C_d*(dk1*(C_hab*x1+C_hbc*x2)+dk2*x12*C_had)-C_a-u1;
    JacX[11] = C_a;

    JacX[14] = C_b;
    JacX[15] = -C_b;
*/

    /**
    ** Bioreactor
    **/
    double x1,x2,x3,x4,x5,x12;
    double u1;
    double mu,num,den,den2,dmu;

    x1 = *xx++; // X
    x2 = *xx++; // S
    x3 = *xx++; // P
    x4 = *xx++;
    x5 = *xx++;
    u1 = *uu++;

    x12 = x1*x2;
    den = B_Km+x2+x2*x2/B_Ki;
    num = (1.-x3/B_Pm)*x2;
    mu = B_mm*num/den;
    den2 = den*den;
    dmu = B_mm*((1.-x3/B_Pm)*den-num*(2.*x2/B_Ki+1))/den2;

    JacX[0] = -B_D+mu;
    JacX[1] = x1*dmu;
    JacX[2] = -x1*B_mm*x2/(B_Pm*den);

    JacX[5] = -mu/B_Yxs;
    JacX[6] = -B_D-x1*dmu/B_Yxs;
    JacX[7] = B_mm*x12/(B_Pm*B_Yxs*den);

    JacX[10] = B_a*mu+B_b;
    JacX[11] = B_a*x1*dmu;
    JacX[12] = -B_D-B_a*B_mm*x12/(B_Pm*den);

    JacX[20] = 1./B_T;

    /**
    ** Walid's quadcopter
    **/
/*    double x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15;
    double u1,u2,u3,u4;
    double cx7,sx7,sx72,cx8,sx8,tx8,cx9,sx9;
    double *xx,*uu;

    xx = x;
    x1 = *xx++;
    x2 = *xx++;
    x3 = *xx++;
    x4 = *xx++;
    x5 = *xx++;
    x6 = *xx++;
    x7 = *xx++;
    x8 = *xx++;
    x9 = *xx++;
    x10 = *xx++;
    x11 = *xx++;
    x12 = *xx++;
    x13 = *xx++;
    x14 = *xx++;
    x15 = *xx++;

    uu = u;
    u1 = *uu++;
    u2 = *uu++;
    u3 = *uu++;
    u4 = *uu++; */

    /* Precomputations */
/*    cx7 = cos(x7);
    sx7 = sin(x7);
    sx72 = sin(0.5*x7);
    cx8 = cos(x8);
    sx8 = sin(x8);
    tx8 = sx8/cx8;
    cx9 = cos(x9);
    sx9 = sin(x9); */

    /* Jacobian matrix */
/*    JacX[3] = 1.0;
    JacX[19] = 1.0;
    JacX[35] = 1.0;
    JacX[51] = u1*(cx7*sx9-cx9*sx7*sx8)/Q_m;
    JacX[52] = u1*cx7*cx8*cx9/Q_m;
    JacX[53] = u1*(cx9*sx7-cx7*sx8*sx9)/Q_m;
    JacX[66] = -u1*(cx7*cx9+sx7*sx8*sx9)/Q_m;
    JacX[67]= u1*cx7*cx8*sx9/Q_m;
    JacX[68] = u1*(sx7*sx9+cx7*cx9*sx8)/Q_m;
    JacX[81] = -u1*cx8*sx7/Q_m;
    JacX[82] = -u1*cx7*sx8/Q_m;
    JacX[96] = (sx8*(x11*cx7-x12*sx7))/cx8;
    JacX[97] = (x12*cx7+x11*sx7)/(cx8*cx8);
    JacX[99] = 1.0;
    JacX[100] = sx7*tx8;
    JacX[101] = cx7*tx8;
    JacX[111] = -x12*cx7-x11*sx7;
    JacX[115] = cx7;
    JacX[116] = -sx7;
    JacX[126] = (x11*cx7-x12*sx7)/cx8;
    JacX[127] = sx8*(x12*(2.0*sx72-1.0)-x11*sx7)/(sx8*sx8-1.0);
    JacX[130] = sx7/cx8;
    JacX[131] = cx7/cx8;
    JacX[141] = u1*(-4.305083284519519E-3);
    JacX[142] = u4*(-3.013558299163663E-1)-x12*1.726641718587006E1+x15*2.152541642259759E-1;
    JacX[144] = u1*(-1.230023795577005E-3);
    JacX[145] = u4*(-8.610166569039038E-2)-x12*6.81114043707411+x15*6.150118977885027E-2;
    JacX[146] = u3*3.013558299163663E-1-x8*1.726641718587006E1-x11*6.81114043707411-x14*1.722033313807808E-1;
    JacX[147] = u1*4.293593022139282E-5;
    JacX[148] = u4*(-3.005515115497497E-3)-x12*1.722033313807808E-1+x15*2.146796511069641E-3;
    JacX[149] = u3*(-3.756893894371871E-3)+x8*2.152541642259759E-1+x11*6.150118977885027E-2+x14*2.146796511069641E-3;
    JacX[156] = u4*(-3.013558299163663E-1)-x12*1.726641718587006E1+x15*2.152541642259759E-1;
    JacX[157] = u1*(-4.305083284519519E-3);
    JacX[159] = u4*(-8.610166569039038E-2)-x12*3.055383669137354+x15*6.150118977885027E-2;
    JacX[160] = u1*(-1.230023795577005E-3);
    JacX[161] = u2*3.013558299163663E-1-x7*1.726641718587006E1-x10*3.055383669137354+x13*1.722033313807808E-1;
    JacX[162] = u4*3.005515115497497E-3+x12*1.722033313807808E-1-x15*2.146796511069641E-3;
    JacX[163] = u1*(-4.293593022139282E-5);
    JacX[164] = u2*(-3.756893894371871E-3)+x7*2.152541642259759E-1+x10*6.150118977885027E-2-x13*2.146796511069641E-3;
    JacX[171] = u3*3.226410088804494E-1-x8*1.848596810669266E1-x11*5.281705173340761-x14*1.843662907888282E-1;
    JacX[172] = u2*3.226410088804494E-1-x7*1.848596810669266E1-x10*5.281705173340761+x13*1.843662907888282E-1;
    JacX[174] = u3*9.218314539441411E-2-x8*5.281705173340761-x11*1.509058620954503-x14*5.267608308252234E-2;
    JacX[175] = u2*9.218314539441411E-2-x7*5.281705173340761-x10*1.509058620954503+x13*5.267608308252234E-2;
    JacX[176] = u1*(-1.505030945214924E-3);
    JacX[177] = u3*(-3.217798803954346E-3)+x8*1.843662907888282E-1+x11*5.267608308252234E-2+x14*1.838742173688198E-3;
    JacX[178] = u2*3.217798803954346E-3-x7*1.843662907888282E-1-x10*5.267608308252234E-2+x13*1.838742173688198E-3;
    JacX[179] = u1*1.87626752417163E-5;
    JacX[186] = -5.729577951308232E1;
    JacX[202] = 5.729577951308232E1;
    JacX[221] = -5.729577951308232E1;
*/

    /* Inverted pendulum */
/*    double x_1,x_2,x_3,x_4,x_4_2;
    double u_1;
    double cx_3,sx_3,csx_3,ccx_3,ssx_3;
    double den,den2,dden;
    double ml,mg,Mmg;
    
    ml = PRM_m*PRM_l;
    mg = PRM_m*PRM_g;
    Mmg = (PRM_M+PRM_m)*PRM_g;

    x_1 = *x; // Horizontal position
    x_2 = *(x+1); // Horizontal speed
    x_3 = *(x+2);  // Angle
    x_4 = *(x+3);   // Angular speed
    x_4_2 = x_4*x_4;
    
    cx_3 = cos(x_3);
    sx_3 = sin(x_3);
    csx_3 = cx_3*sx_3;
    ccx_3 = cx_3*cx_3;
    ssx_3 = sx_3*sx_3;
    
    u_1 = *u;
    
    den = PRM_M+PRM_m-PRM_m*ccx_3;
    den2 = den*den;
    
    dden = 2*PRM_m*csx_3;
    
    JacX[1] = 1.;
    
    JacX[6] = (PRM_m*(PRM_l*cx_3*x_4_2+PRM_g*(ccx_3-ssx_3))*den-(ml*sx_3*x_4_2+mg*csx_3+u_1)*dden)/den2;
    JacX[7] = 2*ml*sx_3*x_4/den;
    
    JacX[11] = 1.;
    
    JacX[14] = (2.*PRM_m*csx_3*(ml*csx_3*x_4_2+u_1*cx_3+Mmg*sx_3)-den*(ml*x_4_2*(ccx_3-ssx_3)-u_1*sx_3+Mmg*cx_3))/(PRM_l*den2);
    JacX[15] = -2.*PRM_m*csx_3*x_4/den;  */

    /* Inverted pendulum with polynmial dynamics */
/*    double x_1,x_2,x_3,x_4,x_5;
    double x_3s,x_5s;
    double u_1;
    double den,dens;

    x_1 = *x; // Horiz. pos. 
    x_2 = *(x+1);  // Horiz. spd.
    x_3 = *(x+2);  // Angle cosine
    x_4 = *(x+3);  // Angle sine
    x_5 = *(x+4); // Angular spd

    u_1 = *u;

    x_3s = x_3*x_3;
    x_5s = x_5*x_5;
    den = PRM_M+PRM_m-PRM_m*x_3s;
    dens = den*den;

    JacX[1] = 1.;

    JacX[7] = ((PRM_M+PRM_m)*PRM_m*PRM_g*x_4+PRM_m*PRM_m*PRM_g*x_3s*x_4+2.*PRM_m*x_3*(PRM_m*PRM_l*x_4*x_5s+u_1))/dens;
    JacX[8] = PRM_m*(PRM_l*x_5s+PRM_g*x_3)/den;
    JacX[9] = 2.*PRM_m*PRM_l*x_4*x_5/den;

    JacX[13] = -x_5;
    JacX[14] = -x_4;

    JacX[17] = x_5;
    JacX[19] = x_3;

    JacX[22] = -((PRM_m*PRM_l*x_4*x_5s+u_1)*(PRM_M+PRM_m-PRM_m*x_3s)+2.*PRM_m*x_3*(PRM_m*PRM_l*x_3*x_4*x_5s+u_1*x_3+(PRM_M+PRM_m)*PRM_g*x_4))/(PRM_l*dens);
    JacX[23] = -(PRM_m*PRM_l*x_3*x_5s+(PRM_M+PRM_m)*PRM_g)/(PRM_l*den);
    JacX[24] = -2.*PRM_m*PRM_l*x_3*x_4*x_5/(PRM_l*den); */

    /* Crane */
/*    double x_3;
    double u_1;

    x_3 = *(x+2); // phi: excitation angle
    u_1 = *u;

    JacX[1] = 1.;    

    JacX[11] = 1.;

    JacX[14] = -1.*PRM_g*cos(x_3)+u_1*sin(x_3);
    JacX[15] = -1.*PRM_b;
*/

    /* DC motor */
/*    double u_1;

    u_1 = *u;

    JacX[0] = -1.*PRM_Ra/PRM_La;
    JacX[1] = -1.*PRM_km*u_1/PRM_La;
    JacX[2] = PRM_km*u_1/PRM_J;
    JacX[3] = -1.*PRM_B/PRM_J; */

    /* Unicycle */
/*    double x_3,u_1,u_2;

    x_3 = *(x+2);

    u_1 = *u;
    u_2 = *(u+1);

    JacX[2] = -1.*u_1*sin(x_3);
    JacX[5] = u_1*cos(x_3); */
}

/**
* Display state-jacobian of dynamics
*/
void MyDynamics::displayJacX(double* x,double* u){

    int i;

    jxRhs(x,u);

    std::cout<<"ODE jacX="<<std::endl;
    for (i=0;i<dx*dx;++i){
        if (i%dx==dx-1){
            std::cout<<JacX[i]<<std::endl;
        }
        else {
            std::cout<<JacX[i]<<" ";
        }
    }
}

/** 
 * Method to evaluate input-jacobian of ode rhs 
 */
void MyDynamics::juRhs(double* x,double* u){
    
    /* USER-DEFINED */
    
    /**
    ** CSTR
    **/
/*    double x1,x2,x3;
    x1 = x[0]; x2 = x[1]; x3 = x[2];

    JacU[0] = C_cin-x1;
    JacU[2] = -x2;
    JacU[4] = C_tin-x3;
    JacU[7] = C_g;
    */

    /**
    ** Bioreactor
    **/
    JacU[1] = B_D;
    JacU[3] = 1./B_T;

    /** 
    ** Walid's quadcopter
    **/
/*    double x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15;
    double u1,u2,u3,u4;
    double cx7,sx7,sx72,cx8,sx8,cx9,sx9;
    double *xx,*uu;

    xx = x;
    x1 = *xx++;
    x2 = *xx++;
    x3 = *xx++;
    x4 = *xx++;
    x5 = *xx++;
    x6 = *xx++;
    x7 = *xx++;
    x8 = *xx++;
    x9 = *xx++;
    x10 = *xx++;
    x11 = *xx++;
    x12 = *xx++;
    x13 = *xx++;
    x14 = *xx++;
    x15 = *xx++;

    uu = u;
    u1 = *uu++;
    u2 = *uu++;
    u3 = *uu++;
    u4 = *uu++;
*/
    /* Precomputations */
/*    cx7 = cos(x7);
    sx7 = sin(x7);
    cx8 = cos(x8);
    sx8 = sin(x8);
    cx9 = cos(x9);
    sx9 = sin(x9);

    JacU[12] = sx7*sx9*3.394433129667346E1+cx7*cx9*sx8*3.394433129667346E1;
    JacU[16] = cx9*sx7*(-3.394433129667346E1)+cx7*sx8*sx9*3.394433129667346E1;
    JacU[20] = cx7*cx8*3.394433129667346E1;
    JacU[36] = u2*7.513787788743744E-5-x7*4.305083284519519E-3-x10*1.230023795577005E-3+x13*4.293593022139282E-5;
    JacU[37] = u1*7.513787788743744E-5;
    JacU[38] = u4*5.25965145212062E-3+x12*3.013558299163663E-1-x15*3.756893894371871E-3;
    JacU[39] = u3*5.25965145212062E-3-x8*3.013558299163663E-1-x11*8.610166569039038E-2-x14*3.005515115497497E-3;
    JacU[40] = u3*7.513787788743744E-5-x8*4.305083284519519E-3-x11*1.230023795577005E-3-x14*4.293593022139282E-5;
    JacU[41] = u4*5.25965145212062E-3+x12*3.013558299163663E-1-x15*3.756893894371871E-3;
    JacU[42] = u1*7.513787788743744E-5;
    JacU[43] = u2*5.25965145212062E-3-x7*3.013558299163663E-1-x10*8.610166569039038E-2+x13*3.005515115497497E-3;
    JacU[44] = u4*(-2.626774533840282E-5)-x12*1.505030945214924E-3+x15*1.87626752417163E-5;
    JacU[45] = u3*(-5.631147906920105E-3)+x8*3.226410088804494E-1+x11*9.218314539441411E-2+x14*3.217798803954346E-3;
    JacU[46] = u2*(-5.631147906920105E-3)+x7*3.226410088804494E-1+x10*9.218314539441411E-2-x13*3.217798803954346E-3;
    JacU[47] = u1*(-2.626774533840282E-5);
    JacU[49] = 1.0;
    JacU[54] = -1.0;
    JacU[59] = -1.0; */

   /* Inverted pendulum */ 
/*    double x_1,x_2,x_3,x_4;
    double u_1;
    double cx_3;
    double den,iden;
    
    x_1 = *x; // Horizontal position
    x_2 = *(x+1); // Horizontal speed
    x_3 = *(x+2); // Angle
    x_4 = *(x+3); // Angular speed
    
    cx_3 = cos(x_3);
    u_1 = *u;
    den = PRM_M+PRM_m-PRM_m*cx_3*cx_3;
    iden = 1./den;

    JacU[1] = iden;
    JacU[3] = -iden*cx_3/PRM_l; 
*/

    /* Inverted pendulum with polynomial dynamics */
/*    double x_1,x_2,x_3,x_4,x_5;
    double u_1;
    double den;

    x_1 = *x; // Horiz. pos.
    x_2 = *(x+1); // Horz. spd.
    x_3 = *(x+2); // Angle cosine
    x_4 = *(x+3); // Angle sine
    x_5 = *(x+4); // Angular speed

    den = PRM_M+PRM_m-PRM_m*x_3*x_3;

    JacU[1] = 1./den;
    JacU[4] = -x_3/(PRM_l*den); */

    /* Crane */
/*    double x_3;
    x_3 = *(x+2); // Excitation angle
    JacU[1] = 1.;
    JacU[3] = -1.*cos(x_3); */

    /* DC motor */
/*    double x_1,x_2;
    x_1 = *x;
    x_2 = *(x+1);
    JacU[0] = -1.*PRM_km*x_2/PRM_La;
    JacU[1] = PRM_km*x_1/PRM_J; */

    /* Unicycle */
/*    double x_3;
    x_3 = *(x+2);
    JacU[0] = cos(x_3);
    JacU[2] = sin(x_3);
    JacU[5] = 1.; */
}

/** 
* Display input-jacobian of dynamics
*/
void MyDynamics::displayJacU(double* x,double* u){

    int i;

    juRhs(x,u);

    std::cout<<"ODE jacU="<<std::endl;
    for (i=0;i<dx*du;++i){
        if (i%du==du-1){
            std::cout<<JacU[i]<<std::endl;
        }
        else {
            std::cout<<JacU[i]<<" ";
        }
    }

}

/**
 * Method to evaluate transpose of state-jacobian of ode rhs 
 */
void MyDynamics::jxRhsT(double* x,double* u){
    
    /* USER-DEFINED */
    double *xx,*uu;
    xx = x;
    uu = u;

    /**
    ** CSTR
    **/
/*    double x1,x2,x3;
    double u1,u2,u3;

    x1 = *xx++;
    x2 = *xx++;
    x3 = *xx++;
    u1 = *uu++;
    u2 = *uu++;

    double den,den2;
    double k1,k2;
    double dk1,dk2;
    double x12;

    den = x3+C_t0;
    den2 = den*den;
    k1 = C_k10*exp(-C_e1/den);
    k2 = C_k20*exp(-C_e2/den);
    dk1 = k1*C_e1/den2;
    dk2 = k2*C_e2/den2;
    x12 = x1*x1;

    JacXt[0] = -k1-2.*k2*x1-u1;
    JacXt[1] = k1;
    JacXt[2] = -C_d*(k1*C_hab+2.*k2*C_had*x1);

    JacXt[5] = -k1-u1;
    JacXt[6] = -C_d*k1*C_hbc;

    JacXt[8] = -dk1*x1-dk2*x12;
    JacXt[9] = dk1*(x1-x2);
    JacXt[10] = -C_d*(dk1*(C_hab*x1+C_hbc*x2)+dk2*x12*C_had)-C_a-u1;
    JacXt[11] = C_b;

    JacXt[14] = C_a;
    JacXt[15] = -C_b;
*/

    /**
    ** Bioreactor
    **/
    double x1,x2,x3,x4,x5,x12;
    double u1;
    double mu,num,den,den2,dmu;

    x1 = *xx++; // X
    x2 = *xx++; // S
    x3 = *xx++; // P
    x4 = *xx++;
    x5 = *xx++;
    u1 = *uu++; // Sf

    x12 = x1*x2;
    den = B_Km+x2+x2*x2/B_Ki;
    num = (1.-x3/B_Pm)*x2;
    mu = B_mm*num/den;
    den2 = den*den;
    dmu = B_mm*((1.-x3/B_Pm)*den-num*(2.*x2/B_Ki+1))/den2;

    JacXt[0] = -B_D+mu;
    JacXt[1] = -mu/B_Yxs;
    JacXt[2] = B_a*mu+B_b;
    JacXt[4] = 1./B_T;

    JacXt[5] = x1*dmu;
    JacXt[6] = -B_D-x1*dmu/B_Yxs;
    JacXt[7] = B_a*x1*dmu;

    JacXt[10] = -x1*B_mm*x2/(B_Pm*den);
    JacXt[11] = B_mm*x12/(B_Pm*B_Yxs*den);
    JacXt[12] = -B_D-B_a*B_mm*x12/(B_Pm*den);

    /** 
    ** Walid's quadcopter
    *
/*    double x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15;
    double u1,u2,u3,u4;
    double cx7,sx7,sx72,cx8,sx8,tx8,cx9,sx9;
    double *xx,*uu;

    xx = x;
    x1 = *xx++;
    x2 = *xx++;
    x3 = *xx++;
    x4 = *xx++;
    x5 = *xx++;
    x6 = *xx++;
    x7 = *xx++;
    x8 = *xx++;
    x9 = *xx++;
    x10 = *xx++;
    x11 = *xx++;
    x12 = *xx++;
    x13 = *xx++;
    x14 = *xx++;
    x15 = *xx++;

    uu = u;
    u1 = *uu++;
    u2 = *uu++;
    u3 = *uu++;
    u4 = *uu++;
*/
    /* Precomputations */
/*    cx7 = cos(x7);
    sx7 = sin(x7);
    sx72 = sin(0.5*x7);
    cx8 = cos(x8);
    sx8 = sin(x8);
    tx8 = sx8/cx8;
    cx9 = cos(x9);
    sx9 = sin(x9);

    JacXt[45] = 1.0;
    JacXt[61] = 1.0;
    JacXt[77] = 1.0;
    JacXt[93] = u1*cx7*sx9*3.394433129667346E1-u1*cx9*sx7*sx8*3.394433129667346E1;
    JacXt[94] = u1*cx7*cx9*(-3.394433129667346E1)-u1*sx7*sx8*sx9*3.394433129667346E1;
    JacXt[95] = u1*cx8*sx7*(-3.394433129667346E1);
    JacXt[96] = tx8*(x11*cx7-x12*sx7);
    JacXt[97] = -x12*cx7-x11*sx7;
    JacXt[98] = (x11*cx7-x12*sx7)/cx8;
    JacXt[99] = u1*(-4.305083284519519E-3);
    JacXt[100] = u4*(-3.013558299163663E-1)-x12*1.726641718587006E1+x15*2.152541642259759E-1;
    JacXt[101] = u3*3.226410088804494E-1-x8*1.848596810669266E1-x11*5.281705173340761-x14*1.843662907888282E-1;
    JacXt[102] = -5.729577951308232E1;
    JacXt[108] = u1*cx7*cx8*cx9*3.394433129667346E1;
    JacXt[109] = u1*cx7*cx8*sx9*3.394433129667346E1;
    JacXt[110] = u1*cx7*sx8*(-3.394433129667346E1);
    JacXt[111] = (x12*cx7+x11*sx7)/(cx8*cx8);
    JacXt[113] = (sx8*(x12*(2.0*sx72*sx72-1.0)-x11*sx7))/(sx8*sx8-1.0);
    JacXt[114] = u4*(-3.013558299163663E-1)-x12*1.726641718587006E1+x15*2.152541642259759E-1;
    JacXt[115] = u1*(-4.305083284519519E-3);
    JacXt[116] = u2*3.226410088804494E-1-x7*1.848596810669266E1-x10*5.281705173340761+x13*1.843662907888282E-1;
    JacXt[118] = 5.729577951308232E1;
    JacXt[123] = u1*cx9*sx7*3.394433129667346E1-u1*cx7*sx8*sx9*3.394433129667346E1;
    JacXt[124] = u1*sx7*sx9*3.394433129667346E1+u1*cx7*cx9*sx8*3.394433129667346E1;
    JacXt[141] = 1.0;
    JacXt[144] = u1*(-1.230023795577005E-3);
    JacXt[145] = u4*(-8.610166569039038E-2)-x12*3.055383669137354+x15*6.150118977885027E-2;
    JacXt[146] = u3*9.218314539441411E-2-x8*5.281705173340761-x11*1.509058620954503-x14*5.267608308252234E-2;
    JacXt[156] = sx7*tx8;
    JacXt[157] = cx7;
    JacXt[158] = sx7/cx8;
    JacXt[159] = u4*(-8.610166569039038E-2)-x12*6.81114043707411+x15*6.150118977885027E-2;
    JacXt[160] = u1*(-1.230023795577005E-3);
    JacXt[161] = u2*9.218314539441411E-2-x7*5.281705173340761-x10*1.509058620954503+x13*5.267608308252234E-2;
    JacXt[171] = cx7*tx8;
    JacXt[172] = -sx7;
    JacXt[173] = cx7/cx8;
    JacXt[174] = u3*3.013558299163663E-1-x8*1.726641718587006E1-x11*6.81114043707411-x14*1.722033313807808E-1;
    JacXt[175] = u2*3.013558299163663E-1-x7*1.726641718587006E1-x10*3.055383669137354+x13*1.722033313807808E-1;
    JacXt[176] = u1*(-1.505030945214924E-3);
    JacXt[179] = -5.729577951308232E1;
    JacXt[189] = u1*4.293593022139282E-5;
    JacXt[190] = u4*3.005515115497497E-3+x12*1.722033313807808E-1-x15*2.146796511069641E-3;
    JacXt[191] = u3*(-3.217798803954346E-3)+x8*1.843662907888282E-1+x11*5.267608308252234E-2+x14*1.838742173688198E-3;
    JacXt[204] = u4*(-3.005515115497497E-3)-x12*1.722033313807808E-1+x15*2.146796511069641E-3;
    JacXt[205] = u1*(-4.293593022139282E-5);
    JacXt[206] = u2*3.217798803954346E-3-x7*1.843662907888282E-1-x10*5.267608308252234E-2+x13*1.838742173688198E-3;
    JacXt[219] = u3*(-3.756893894371871E-3)+x8*2.152541642259759E-1+x11*6.150118977885027E-2+x14*2.146796511069641E-3;
    JacXt[220] = u2*(-3.756893894371871E-3)+x7*2.152541642259759E-1+x10*6.150118977885027E-2-x13*2.146796511069641E-3;
    JacXt[221] = u1*1.87626752417163E-5;
*/
    /* Inverted pendulum */
/*    double x_1,x_2,x_3,x_4,x_4_2;
    double u_1;
    double cx_3,sx_3,csx_3,ccx_3,ssx_3;
    double den,den2,dden;
    double ml,mg,Mmg;
    
    ml = PRM_m*PRM_l;
    mg = PRM_m*PRM_g;
    Mmg = (PRM_M+PRM_m)*PRM_g;

    x_1 = *x;
    x_2 = *(x+1);
    x_3 = *(x+2);
    x_4 = *(x+3);
    x_4_2 = x_4*x_4;
    
    cx_3 = cos(x_3);
    sx_3 = sin(x_3);
    csx_3 = cx_3*sx_3;
    ccx_3 = cx_3*cx_3;
    ssx_3 = sx_3*sx_3;
    
    u_1 = *u;
    den = PRM_M+PRM_m-PRM_m*ccx_3;
    den2 = den*den;
    dden = 2.*PRM_m*csx_3;

    JacXt[4] = 1.;    

    JacXt[9] = (PRM_m*(PRM_l*cx_3*x_4_2+PRM_g*(ccx_3-ssx_3))*den-(ml*sx_3*x_4_2+mg*csx_3+u_1)*dden)/den2;
    JacXt[11] = (2.*PRM_m*csx_3*(ml*csx_3*x_4_2+u_1*cx_3+Mmg*sx_3)-den*(ml*x_4_2*(ccx_3-ssx_3)-u_1*sx_3+Mmg*cx_3))/(PRM_l*den2);

    JacXt[13] = 2.*ml*sx_3*x_4/den;
    JacXt[14] = 1.;
    JacXt[15] = -2.*PRM_m*csx_3*x_4/den; */

    /* Inverted pendulum with polynomial dynamics */
/*    double x_1,x_2,x_3,x_4,x_5;
    double x_3s,x_5s;
    double u_1;
    double den,dens;

    x_1 = *x; // Horiz. pos. 
    x_2 = *(x+1);  // Horiz. spd.
    x_3 = *(x+2);  // Angle cosine
    x_4 = *(x+3);  // Angle sine
    x_5 = *(x+4); // Angular spd

    u_1 = *u; 

    x_3s = x_3*x_3;
    x_5s = x_5*x_5;
    den = PRM_M+PRM_m-PRM_m*x_3s;
    dens = den*den;

    JacXt[5] = 1.;

    JacXt[11] = ((PRM_M+PRM_m)*PRM_m*PRM_g*x_4+PRM_m*PRM_m*PRM_g*x_3s*x_4+2.*PRM_m*x_3*(PRM_m*PRM_l*x_4*x_5s+u_1))/dens;
    JacXt[14] = -((PRM_m*PRM_l*x_4*x_5s+u_1)*(PRM_M+PRM_m-PRM_m*x_3s)+2.*PRM_m*x_3*(PRM_m*PRM_l*x_3*x_4*x_5s+u_1*x_3+(PRM_M+PRM_m)*PRM_g*x_4))/(PRM_l*dens);

    JacXt[16] = PRM_m*(PRM_l*x_5s+PRM_g*x_3)/den;
    JacXt[17] = -x_5;
    JacXt[18] = x_5;
    JacXt[19] = -(PRM_m*PRM_l*x_3*x_5s+(PRM_M+PRM_m)*PRM_g)/(PRM_l*den);

    JacXt[21] = 2.*PRM_m*PRM_l*x_4*x_5/den;
    JacXt[22] = -x_4;
    JacXt[23] = x_3;
    JacXt[24] = -2.*PRM_m*PRM_l*x_3*x_4*x_5/(PRM_l*den); */


    /* Crane */
/*    double x_3;
    double u_1;

    x_3 = *(x+2); // phi: excitation angle

    u_1 = *u;

    JacXt[4] = 1.;

    JacXt[11] = -1.*PRM_g*cos(x_3)+u_1*sin(x_3);

    JacXt[14] = 1.;
    JacXt[15] = -1.*PRM_b; */


    /* DC motor */
/*    double u_1;

    u_1 = *u;

    JacXt[0] = -1.*PRM_Ra/PRM_La;
    JacXt[1] = PRM_km*u_1/PRM_J;
    JacXt[2] = -1.*PRM_km*u_1/PRM_La;
    JacXt[3] = -1.*PRM_B/PRM_J; */

    /* Unicycle */
/*    double x_3,u_1,u_2;

    x_3 = *(x+2);

    u_1 = *u;
    u_2 = *(u+1);

    JacXt[6] = -1.*u_1*sin(x_3);
    JacXt[7] = u_1*cos(x_3); */
}

void MyDynamics::displayJacXt(double* x,double* u){

    int i;

    jxRhsT(x,u);

    std::cout<<"ODE jacXt="<<std::endl;
    for (i=0;i<dx*dx;++i){
        if (i%dx==dx-1){
            std::cout<<JacXt[i]<<std::endl;
        }
        else {
            std::cout<<JacXt[i]<<" ";
        }
    }

}

/** 
 * Method to evaluate tranpose of input-jacobian of ode rhs 
 */
void MyDynamics::juRhsT(double* x,double* u){
    
    /* USER-DEFINED */
    /**
    ** CSTR
    **/
/*    double x1,x2,x3;
    x1 = x[0]; x2 = x[1]; x3 = x[2];

    JacUt[0] = C_cin-x1;
    JacUt[1] = -x2;
    JacUt[2] = C_tin-x3;
    JacUt[7] = C_g;
    */

    /**
    ** Bioreactor
    **/
    JacUt[1] = B_D;
    JacUt[3] = 1./B_T;

    /**
    ** Walid's quadcopter
    **/
/*    double x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15;
    double u1,u2,u3,u4;
    double cx7,sx7,sx72,cx8,sx8,cx9,sx9;
    double *xx,*uu;

    xx = x;
    x1 = *xx++;
    x2 = *xx++;
    x3 = *xx++;
    x4 = *xx++;
    x5 = *xx++;
    x6 = *xx++;
    x7 = *xx++;
    x8 = *xx++;
    x9 = *xx++;
    x10 = *xx++;
    x11 = *xx++;
    x12 = *xx++;
    x13 = *xx++;
    x14 = *xx++;
    x15 = *xx++;

    uu = u;
    u1 = *uu++;
    u2 = *uu++;
    u3 = *uu++;
    u4 = *uu++; */

    /* Precomputations */
/*    cx7 = cos(x7);
    sx7 = sin(x7);
    cx8 = cos(x8);
    sx8 = sin(x8);
    cx9 = cos(x9);
    sx9 = sin(x9);

    JacUt[3] = sx7*sx9*3.394433129667346E1+cx7*cx9*sx8*3.394433129667346E1;
    JacUt[4] = cx9*sx7*(-3.394433129667346E1)+cx7*sx8*sx9*3.394433129667346E1;
    JacUt[5] = cx7*cx8*3.394433129667346E1;
    JacUt[9] = u2*7.513787788743744E-5-x7*4.305083284519519E-3-x10*1.230023795577005E-3+x13*4.293593022139282E-5;
    JacUt[10] = u3*7.513787788743744E-5-x8*4.305083284519519E-3-x11*1.230023795577005E-3-x14*4.293593022139282E-5;
    JacUt[11] = u4*(-2.626774533840282E-5)-x12*1.505030945214924E-3+x15*1.87626752417163E-5;
    JacUt[24] = u1*7.513787788743744E-5;
    JacUt[25] = u4*5.25965145212062E-3+x12*3.013558299163663E-1-x15*3.756893894371871E-3;
    JacUt[26] = u3*(-5.631147906920105E-3)+x8*3.226410088804494E-1+x11*9.218314539441411E-2+x14*3.217798803954346E-3;
    JacUt[27] = 1.0;
    JacUt[39] = u4*5.25965145212062E-3+x12*3.013558299163663E-1-x15*3.756893894371871E-3;
    JacUt[40] = u1*7.513787788743744E-5;
    JacUt[41] = u2*(-5.631147906920105E-3)+x7*3.226410088804494E-1+x10*9.218314539441411E-2-x13*3.217798803954346E-3;
    JacUt[43] = -1.0;
    JacUt[54] = u3*5.25965145212062E-3-x8*3.013558299163663E-1-x11*8.610166569039038E-2-x14*3.005515115497497E-3;
    JacUt[55] = u2*5.25965145212062E-3-x7*3.013558299163663E-1-x10*8.610166569039038E-2+x13*3.005515115497497E-3;
    JacUt[56] = u1*(-2.626774533840282E-5);
    JacUt[59] = -1.0;
    */

    /* Inverted pendulum */
/*    double x_1,x_2,x_3,x_4;
    double u_1;
    double cx_3;
    double den,iden;
    
    x_1 = *x;
    x_2 = *(x+1);
    x_3 = *(x+2);
    x_4 = *(x+3);
    
    cx_3 = cos(x_3);
    u_1 = *u;
    den = PRM_M+PRM_m-PRM_m*cx_3*cx_3;
    iden = 1./den;
    
    JacUt[1] = iden;
    JacUt[3] = -iden*cx_3/PRM_l; */

    /* Inverted pendulum with polynomial dynamics */
/*    double x_1,x_2,x_3,x_4,x_5;
    double u_1;
    double den;

    x_1 = *x; // Horiz. pos.
    x_2 = *(x+1); // Horz. spd.
    x_3 = *(x+2); // Angle cosine
    x_4 = *(x+3); // Angle sine
    x_5 = *(x+4); // Angular speed

    den = PRM_M+PRM_m-PRM_m*x_3*x_3;

    JacUt[1] = 1./den;
    JacUt[4] = -x_3/(PRM_l*den); */

    /* Crane */
/*    double x_3;

    x_3 = *(x+2); // Excitation angle

    JacUt[1] = 1.;
    JacUt[3] = -1.*cos(x_3); */

    /* DC motor */
/*    double x_1,x_2;

    x_1 = *x;
    x_2 = *(x+1);
    JacUt[0] = -1.*PRM_km*x_2/PRM_La;
    JacUt[1] = PRM_km*x_1/PRM_J; */

    /* Unicycle */
/*    double x_3;

    x_3 = *(x+2);

    JacUt[0] = cos(x_3);
    JacUt[1] = sin(x_3);
    JacUt[5] = 1.; */
}

void MyDynamics::displayJacUt(double* x,double* u){

    int i;

    juRhsT(x,u);

    std::cout<<"ODE jacUt="<<std::endl;
    for (i=0;i<dx*du;++i){
        if (i%dx==dx-1){
            std::cout<<JacUt[i]<<std::endl;
        }
        else {
            std::cout<<JacUt[i]<<" ";
        }
    }
}

/**
** Method to evaluate adjoint ODE rhs
**   (\dot{\lmb}_x',\dot{\lmb}_u')' = -(Jx(x,u)',Ju(x,u)')'*lmb_x
**/
void MyDynamics::adjRhs(double* arhs,double* x,double* u,double* adj){

    int i,j,jj;
    double aux,*prd;
 
    /* Eval. state-jacobian transpose */
    jxRhsT(x,u);   
    /* Eval. input-jacobian transpose */
    juRhsT(x,u);
       
    /**
    * Eval. adjoint dynamics
    */
    prd = arhs;
    /* State-adjoint */
    jj = 0;
    for (i=0;i<dx;++i){
        /* State-jacobian transpose */
        aux = 0.;
        for (j=0;j<dx;++j){
            aux -= *(JacXt+jj)**(adj+j);
            ++jj;
        }
        *prd++ = aux;
    }
    /* Input-adjoint */
    jj = 0;
    for (i=0;i<du;++i){
        /* Input-jacobian transpose */
        aux = 0.;
        for (j=0;j<dx;++j){
            aux -= *(JacUt+jj)**(adj+j);
            ++jj;
        }
        *prd++ = aux;
    }
}