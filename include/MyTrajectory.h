//
//  MyTrajectory.h
//  MusAL
//
//  Created by Jean on 4/9/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_MyTrajectory_h
#define MusAL_MyTrajectory_h


#include "userData.h"


class MyTrajectory {


public:
    /* Default constructor */
    MyTrajectory();
    
    /* Copy constructor */
    MyTrajectory(const MyTrajectory&);
    
    /* Assignment operator */
    MyTrajectory& operator= (const MyTrajectory&);
    
    /* Destructor */
    ~MyTrajectory();
    
    
    /**
     * Set methods
    */
     /* Set initial time */
     inline void setInitialTime(double ti){
        this->ti = ti;
     }

     /* Set final time */
     inline void setFinalTime(double tf){
        this->tf = tf;
     }

    /* Set sampling period */
    inline void setSamplingPeriod(double ts){
        this->ts = ts;
    }
    
    /* Set # samples */
    inline void setNumberSamples(int Ns){
        this->Ns = Ns;
    }


    /**
    * Get methods 
     */
    inline double* getOutputReference(){
        return yRef;
    }
    
    inline double* getInputReference(){
        return uRef;
    }
    
    inline int getNs(){
        return Ns;
    }

    inline int getDy(){
        return dy;
    }

    inline int getDu(){
        return du;
    }
    

    /**
    * Processing methods
    */
    /* Compute output reference between ti & tf */
    void computeOutputReference();
    
    /* Compute input reference between ti & tf */
    void computeInputReference();
    
    /* Display output & input reference */
    void displayReference();

    /**
    * User-defined methods 
    */
    /* Compute output reference point at time instant */
    void computeOutputPoint(double,double*);

    /* Compute input reference point at time instant */
    void computeInputPoint(double,double*);
    

private:

    /* Dimensions */
    const static int du = DU;
    const static int dy = DY;

    /* Final time */
    double tf;
    
    /* Initial time */
    double ti;

    /* Sampling period */
    double ts;
    
    /* Number of reference points */
    int Ns;

    /* Output reference */
    double* yRef;
    
    /* Input reference */
    double* uRef;
    
};


#endif
