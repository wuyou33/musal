//
//  MyPathConstraint.h
//  MusAL
//
//  Created by Jean on 3/14/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_MyPathConstraint_h
#define MusAL_MyPathConstraint_h

#include <iostream>
#include <algorithm>
#include <math.h>

#include "macros.h"
#include "userData.h"

/****
* Implements path-constraint g(x,u)<=0, transformed into g(x,u)+r, r is a nonnegative slack variable
**/
class MyPathConstraint {


public:
    /* Default constructor */
    MyPathConstraint();
    
    /* Destructor */
    ~MyPathConstraint();
    
    /* Evaluate path-constraint */
    void eval(double*,double*,double*,double*);
    
    /* Evaluate transpose of state-jacobian of path-constraint */
    void evalJacXt(double*,double*);
    
    /* Evaluate transpose of input-jacobian of path-constraint */
    void evalJacUt(double*,double*);
    
    /* Evaluate product of vector against state-jacobian transpose */
    void prdJacXt(double*,double*);
    
    /* Evaluate product of vector against input-jacobian transpose */
    void prdJacUt(double*,double*);
    
    /* Gets */
    inline int getDx(){
        return dx;
    }
    
    inline int getDu(){
        return du;
    }
    
    inline int getDpc(){
        return dpc;
    }
    
    inline double* getJacXt(){
        return jacXt;
    }
    
    inline double* getJacUt(){
        return jacUt;
    }
    
    
private:
    /* State dimension, USER DEFINED */
    const static int dx = DX;
    
    /* Input dimension, USER DEFINED */
    const static int du = DU;
    
    /* Constraint dimension */
    const static int dpc = DPC;
    
    /* State jacobian */
    double jacXt[dx*dpc];
    
    /* Input jacobian */
    double jacUt[du*dpc];
    
};

    

#endif
