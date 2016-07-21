//
//  ExplicitRKintegrator.h
//  MusAL
//
//  Created by Jean on 3/12/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_ExplicitRKintegrator_h
#define MusAL_ExplicitRKintegrator_h

#include <iostream>
#include <algorithm>
#include <assert.h>


template <class dynT> class ExplicitRKintegrator {
    
public:
    /* Default constructor */
    ExplicitRKintegrator();
    
    /* Constructor from # stages, Butcher tableau, state & prm dimensions, max. # steps & bounds on step-size */
//    ExplicitRKintegrator(int,double*,double*,double*,int,int,int,double,double);
    
    /* Copy constructor */
    ExplicitRKintegrator(const ExplicitRKintegrator<dynT>&);
    
    /* Destructor */
    ~ExplicitRKintegrator();
    
    /* Overloaded assignment operator */
    ExplicitRKintegrator<dynT> & operator= (const ExplicitRKintegrator<dynT>&);
    
    /**
    * Set methods 
    */
    /* Set pointer to initial state */
    inline void setStateIni(double* xini){
        this->xini = xini;
    }
    
    /* Set pointer to final state */
    inline void setStateFin(double* xfin){
        this->xfin = xfin;
    }
    
    /* Set pointer to intermediate states */
    inline void setStates(double* states){
        this->states = states; 
    }

    /* Set pointer to interpolated states */
    inline void setStatesInterp(double* statesInterp){
        this->statesInterp = statesInterp;
    }

    /* Set pointer to final adjoint */
    inline void setAdjointFin(double* adfin){
        this->adfin = adfin;
    }
    
    /* Set pointer to initial adjoint */
    inline void setAdjointIni(double* adini){
        this->adini = adini;
    }
    
    /* Set pointer to input */
    inline void setInput(double* u){
        this->u = u;
    }
    
    /* Set pointer to output-reference, called when templated on AugODEquation */
    inline void setOutputRef(double* yref){
        ode.setOutputRef(yref);
    }
    
    /* Set pointer to input reference, called when template on AugODEquation */
    inline void setInputRef(double* uref){
        ode.setInputRef(uref);
    }
    
   	/* Set integration time-interval */
	inline void setTimeInterval(double timeInterval){
		this->timeInterval = timeInterval;
	}
    
    /* Set fixed step-size */
    inline void setHfix(double hfix){
        assert(hfix>0.);
    	this->hfix = hfix;
    }
    
    /* Set initial step-size for adaptive step-size integration */
    inline void setHini(double hini){
        this->hini = hini;
    }
    
    /* Set minimum step-size for adaptive integration */
    inline void setHmin(double hmin){
    	this->hmin = hmin;
    }
    
    /* Set maximum step-size for adaptive integration */
    inline void setHmax(double hmax){
    	this->hmax = hmax;
    }
    
    /* Set # steps */
    inline void setNsteps(int nsteps){
        assert(nsteps>=1);
        this->nsteps = nsteps;
    }
    
    /* Set ODE */
    inline void setODE(const dynT& ode){
        this->ode = ode;
    }

    /* Set inverse input scaling */
    inline void setInputScaling(double* Dq){
        ode.setInputScaling(Dq);
    }
    
    /* Set Butcher Tableau */
    void setButcher(int,double*,double*,double*);
    
    /* Set integration step-sizes,... */
    void setIntegSteps(int,int,double,double,double,double);

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

    inline int getIntegrationOrder(){
        return s;
    }
    
    inline double* getA(){
        return A;
    }
    
    inline double* getb(){
        return b;
    }
    
    inline double* getc(){
        return c;
    }
    
    inline double* getStateNodes(){
        return stateNodes;
    }
    
    inline double* getAdjointNodes(){
        return adjointNodes;
    }
    
    inline double* getInterpolatedStates(){
        return statesInterp;
    }
    
    inline double getTimeInterval(){
        return timeInterval;
    }
    
    inline int getNsteps(){
        return nsteps;
    }

    /***************** Processing methods *************************/
    
    /* Compute scaled input */
    void applyInvInputScaling();

    /* Initialization routine */
    void init();
    
    /* Perform one state integration step */
    void performStepState(double*,double*,double*);
    
    /* Perform one adjoint integration step */
    void performStepAdjoint(double*,double*,double*);
    
    /* Run state integration loop with fixed step-size */
    void integrateStateFixed();
    
    /* Run adjoint integration loop */
    void integrateAdjointFixed();
    
    /* Compute state interpolates between all integration steps */
    void interpState();
    
    /* Run state integration loop with step-size control */
    void integrateStateAdaptive();
    
    /* Run adjoint integration loop */
    void integrateAdjointAdaptive();
    
    /***************** Display methods ******************/
    
    /* Display integration parameters */
    void displayParameters();

    /* Display Butcher tableau */
    void displayButcher();

    /* Display interpolated states */
    void displayInterpolates();

    /* Display states */
    void displayStates();
    
private:
	/******************** Dimensions **************************/
    /* State dimension */
    int dx;
    int dxbl;
    
    /* Input dimension */
    int du;
    int dubl;
    
    /* Adjoint dimension */
    int da;
    int dabl;
    
    /* Useful dimensions */
    int sdx;
    int sdxnsteps;
    int dxnsteps;
    
    /********************* Integrator parameters **************************/
    /* Fixed step-size */
    double hfix;
    
    /* Initial step-size for adaptive step-size integration */
    double hini;
    
    /* Max. step-size */
    double hmax;
    
    /* Min. step-size */
    double hmin;
    
    /* State step-size history (for adaptive step integration) ALLOCATED */
    double* hsta;
    
    /* Adjoint step-size history (for adaptive step integration) ALLOCATED */
    double* hadj;
    
    /* # integration steps */
    int nsteps;
    
    /* Max. # integration steps */
    int nstepsMax;
    
    /* # stages */
    int s;
    
    /* Butcher-A matrix ALLOCATED*/
    double* A;
    
    /* Butcher-b matrix ALLOCATED */
    double* b;
    
    /* Butcher-c matrix ALLOCATED */
    double* c;
    
    /************************* Integrator data (to be allocated locally) ****************************/
    /* ODE to be integrated */
    dynT ode;
    
    /* Integration interval */
    double timeInterval;
        
    /* Adjoints over integration interval (backwards in time) ALLOCATED */    
    double* adjoints;
        
    /* Intermediate states for forward integration ALLOCATED */
    double* stateNodes;
    
    /* Intermediate adjoints for backward integration ALLOCATED */
    double* adjointNodes;
    
    /* Auxiliary variables for product agains Butcher A-matrix ALLOCATED */
    double* auxState;
    double* auxAdjoint;
    
    /******************* Pointers ***********************/
    /* Pointer to initial state */
    double* xini;
    
    /* Pointer to final state */
    double* xfin;
    
    /* Pointer to final adjoint */
    double* adfin;
    
    /* Pointer to initial adjoint */
    double* adini;
    
    /* Pointer to input */
    double* u;

    /* Pointer to states over integration interval */
    double* states;

    /* Pointer to interpolated states on integration interval */
    double* statesInterp;
};
#include "ExplicitRKintegrator.tpp"


#endif
