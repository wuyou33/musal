//
//  MyTrackCost.h
//  MusAL
//
//  Created by Jean on 3/13/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_MyTrackCost_h
#define MusAL_MyTrackCost_h


#include <iostream>
#include <math.h>

#include "userData.h"


class MyTrackCost {
    
    
public:
    /* Default constructor */
    MyTrackCost();
    
    /* Destructor */
    ~MyTrackCost();
   
    /**
    ** USER-DEFINED methods
    **/
    /* Compute stage-cost */
    void evalStage(double*,double*,double*,double*,double*);
    
    /* Compute stage-cost gradient wrt state */
    void evalGradXstage(double*,double*,double*,double*);
    
    /* Compute stage-cost gradient wrt input */
    void evalGradUstage(double*,double*,double*,double*);
    
    /** 
    ** Get methods 
    **/
    inline int getDx(){
        return dx;
    }
    
    inline int getDu(){
        return du;
    }
    
    inline int getDy(){
    	return dy;
    }
    
    inline double* getQ(){
        return Q;
    }
    
    inline double* getR(){
        return R;
    }
    
    inline double* getGx(){
        return gX;
    }
    
    inline double* getGu(){
        return gU;
    }
    

private:
    
    /* State dimension, USER-DEFINED */
    const static int dx = DX;
    
    /* Output dimension, USER-DEFINED */
    const static int dy = DY;
    
    /* Input dimension, USER-DEFINED */
    const static int du = DU;
    
    /* USER-DEFINED output-weighting matrix (diagonal) */
    double Q[DY];
    
    /* USER-DEFINED input-weighting matrix (diagonal) */
    double R[DU];
    
    /* State gradient */
    double gX[DX];
    
    /* Input gradient */
    double gU[DU];
    
};
    

#endif
