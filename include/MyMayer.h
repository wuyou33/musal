//
//  MyMayer.h
//  MusAL
//
//  Created by Jean on 3/14/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_MyMayer_h
#define MusAL_MyMayer_h

#include <iostream>
#include <math.h>

#include "userData.h"

/***********************************
** Implements Mayer term in OCP (terminal penalty)
********************************/
class MyMayer {

public:
    /* Default constructor */
    MyMayer();
    
    /* Destructor */
    ~MyMayer();
    
    /* Evaluate Mayer term */
    double eval(double*);
    
    /* Evaluate gradient of Mayer term */
    void grad(double*);
    
    /* Set methods */
    inline void setStateScaling(double* Ds){
        this->Ds = Ds;
    }

    inline void setStateCentering(double* cs){
        this->cs = cs;
    }

    /* Get methods */
    inline double* getGrad(){
        return gX;
    }
    
    inline int getDx(){
        return dx;
    }
    
private:
    
    /* State dimension, USER DEFINED */
    const static int dx = DX;
    
    /* Gradient */
    double gX[dx];

    /* Scaled terminal state */
    double xsc[dx];

    /* State scaling (from [0,1] to true set) */
    double* Ds;

    /* State centering */
    double* cs;
};


#endif
