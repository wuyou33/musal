//
//  MyDynamics.h
//  MusAL
//
//  Created by Jean on 3/12/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_MyDynamics_h
#define MusAL_MyDynamics_h

#include <iostream>
#include <math.h>

#include "macros.h"
#include "userData.h"

/**
* Implements ODE rhs (dynamics)
*/
class MyDynamics {

public:
    /* Default constructor */
    MyDynamics();
    
    /* Destructor */
    ~MyDynamics();
    
    /* Method to evaluate ode rhs */
    void rhs(double*,double*,double*);
    
    /* Method to evaluate state-jacobian of ode rhs */
    void jxRhs(double*,double*);

    /* Method to evaluate input-jacobian of ode rhs */
    void juRhs(double*,double*);
    
    /* Method to evaluate transpose of state-jacobian of ode rhs */
    void jxRhsT(double*,double*);
    
    /* Method to evaluate tranpose of input-jacobian of ode rhs */
    void juRhsT(double*,double*);
    
    /* Method to evaluate adjoint ODE rhs */
    void adjRhs(double*,double*,double*,double*);
    
    /**
     * Set methods 
    */
    inline void setInputScaling(double* Dq){
        // Empty, just for template matching with AugODEquation
    }

    inline void setOutputRef(double* yref){
    	// Empty, just for template matching with AugODEquation
    }

	inline void setInputRef(double* uref){
    	// Empty, just for template matching with AugODEquation
    }

    /** 
    * Get methods 
    */
    inline int getDx(){
        return dx;
    }
    
    inline int getDu(){
        return du;
    }
    
    inline int getDa(){
        return da;
    }
    
    inline double* getJacX(){
        return JacX;
    }
    
    inline double* getJacU(){
        return JacU;
    }
    
    inline double* getJacXt(){
        return JacXt;
    }
    
    inline double* getJacUt(){
        return JacUt;
    }

    /** 
    * Display methods
    */
    void displayRhs(double*,double*);

    void displayJacX(double*,double*);

    void displayJacU(double*,double*);

    void displayJacXt(double*,double*);

    void displayJacUt(double*,double*);

private:
    /* State dimension, USER DEFINED */
    const static int dx = DX;
    
    /* Input dimension, USER DEFINED */
    const static int du = DU;
    
    /* Adjoint dimension */
    const static int da = dx+du;
    
    /* State-jacobian and its transpose */
    double JacX[DX*DX];
    double JacXt[DX*DX];
    
    /* Input-jacobian and its transpose */
    double JacU[DX*DU];
    double JacUt[DU*DX];
};

#endif
