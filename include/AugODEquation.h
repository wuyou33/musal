//
//  AugODEquation.h
//  MusAL
//
//  Created by Jean on 3/12/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_AugODEquation_h
#define MusAL_AugODEquation_h

#include <iostream>
#include <algorithm>


template <class dynT,class cosT> class AugODEquation {
    
public:    
    /* Default constructor */
    AugODEquation();
    
    /* Destructor */
    ~AugODEquation();

    /* Copy constructor */
    AugODEquation(const AugODEquation&);

    /* Assignment operator */
    AugODEquation& operator= (const AugODEquation&);
   
    /**
    ** Processing methods
    **/
    /* Method to evaluate augmented ode rhs */
    void rhs(double*,double*,double*);
    
    /* Method to evaluate adjoint ode rhs */
    void adjRhs(double*,double*,double*,double*);

	/**
	** Set methods 
	**/
    /* Set ode dynamics */
    inline void setDynamics(const dynT& ode){
    	this->ode = ode;
        ds = ode.getDx();
    	dx = ds+1;
    	du = ode.getDu();
    	da = dx+du;
        /* Set pointers to jacobians transpose */
        dynJxT = ode.getJacXt();
        dynJuT = ode.getJacUt();
    }
    
    /* Set stage-cost */
    inline void setStageCost(const cosT& cost){
        this->cost = cost;
        dy = cost.getDy();
        /* Set pointers to gradients */
        cosGx = cost.getGx();
        cosGu = cost.getGu();
    }
    
    /* Set pointer to output reference */
    inline void setOutputRef(double* yref){
    	this->yref = yref;
    }
    
    /* Set pointer to input reference */
    inline void setInputRef(double* uref){
    	this->uref = uref;
    }

    /* Set inverse input scaling */
    inline void setInputScaling(double* Dq){
        this->Dq = Dq;
    }

    /** 
    ** Gets 
    **/
    inline int getDx(){
        return dx;
    }
    
    inline int getDs(){
        return ds;
    }
    
    inline int getDy(){
    	return dy;
   	}	 
    
    inline int getDa(){
        return da;
    }
    
    inline int getDu(){
        return du;
    }
    
    inline dynT& getDynamics(){
        return ode;
    }
    
    inline cosT& getStageCost(){
        return cost;
    }
    
private:
	/* Dynamics */
    dynT ode;
    
    /* Stage-cost */
    cosT cost;
    
    /**
    * Dimensions 
    */
    /* Augmented-state dimension (ds+1) */
    int dx;
    
    /* Augmented-adjoint dimension (dx+du) */
    int da;
    
    /* Parameter (input) dimension */
    int du;
    
    /* State dimension */
    int ds;
    
    /* Output dimension */
    int dy;
    
    /** 
    * Pointers 
    */
    /* Pointer to output reference */
    double* yref;
    
    /* Pointer to input reference */
    double* uref;
    
    /* Pointer to state jacobian transpose (dynamics) */
    double* dynJxT;
    
    /* Pointer to state gradient (stage-cost) */
    double* cosGx;
    
    /* Pointer to input jacobian transpose (dynamics) */
    double* dynJuT;
    
    /* Pointer to input gradient (stage-cost) */
    double* cosGu;

    /* Input scaling matrix */
    double* Dq;
};
#include "AugODEquation.tpp"


#endif
