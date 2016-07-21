//
//  MuShoot.h
//  MusAL
//
//  Created by Jean on 3/12/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_MuShoot_h
#define MusAL_MuShoot_h

#include <iostream>
#include <cassert>

#include "AugODEquation.h"
#include "ExplicitRKintegrator.h"


/*** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
* Implements multiple-shooting discretization of OCP:
* 	- outputs augmented Lagrangian value
*	- outputs partially separable gradient of augmented Lagrangian
*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***/
template <class dynT,class costT,class pconT,class mayerT> class MuShoot {
    
    /* Augmented ODE type */
    typedef AugODEquation<dynT,costT> AugOde;
    
public:
    /* Default constructor */
    MuShoot();
    
    /* Copy constructor */
    MuShoot(const MuShoot<dynT,costT,pconT,mayerT>&);
    
    /* Destructor */
    ~MuShoot();
    
    /* Overloaded assigment operator */
    MuShoot<dynT,costT,pconT,mayerT>& operator= (const MuShoot<dynT,costT,pconT,mayerT>&);
    
    /**
     * Set methods
     */
    /* Set prediction time */
    inline void setPredictionTime(double Tpred){
        assert(Tpred>0.);
        this->Tpred = Tpred;
    }
    
    /* Set sampling period */
    inline void setSamplingPeriod(double ts){
        assert(ts>0.);
        erkInteg.setTimeInterval(ts);
    }

    /* Set # shooting nodes */
    inline void setNumberShootingNodes(int Ns){
        assert(Ns>=2);
        this->Ns = Ns;
        this->Ns_ = Ns-1;
    }

    /* Set # shooting intervals */
    inline void setNumberShootingIntervals(int Ns_){
        assert(Ns_>=1);
        this->Ns_ = Ns_;
        this->Ns = Ns_+1;
    }

    /* Set fixed integration step-size */
    inline void setIntegrationStepSize(double hfix){
        assert(hfix>0.);
    	erkInteg.setHfix(hfix);
    }
        
    /* Set # integration steps */
    inline void setNumberIntegrationSteps(int nsteps){
        assert(nsteps>=1);
        erkInteg.setNsteps(nsteps);
    }

    /* Set Butcher tableau given by s,A,b,c */
    inline void setButcher(int s,double* A,double* b,double* c){
        assert(s>0);
        erkInteg.setButcher(s,A,b,c);
    }
    
	/* Set ERK integrator */
	inline void setERKinteg(const ExplicitRKintegrator<AugOde>& erkInteg){
        this->erkInteg = erkInteg;
    }

    /* Set pointer to initial state */
    inline void setXest(double* xest){
        this->xest = xest;
#ifdef BIOREC
        xest[3] = 0.;
        xest[4] = 0.;
#endif
    }
    
    /* Set output & input references */
    void setReference(double*,double*);
    
    /* Set pointer to primal optimizer */
    inline void setPrimal(double* z){
        this->z = z;
    }
    
    /* Set pointer to candidate */
    inline void setCandidate(double* zS){
        this->zS = zS;
    }
    
    /* Set pointer to dual optimizer */
    inline void setDual(double* mu){
        this->mu = mu;
    }
    
    /* Set pointer to equality constraints */
    inline void setEcon(double* econ){
        this->econ = econ;
    }
    
    /* Set pointer to full AL gradient */
    inline void setGradient(double* gAL){
        this->gAL = gAL;
    }
    
    /* Set pointer to shooting part of AL gradient at iterate */
    inline void setGradientShooting(double* gALshoot){
        this->gALshoot = gALshoot;
    }
    
    /* Set pointer to shooting part of AL gradient at candidate */
    inline void setGradientShootingCandidate(double* gALshootS){
        this->gALshootS = gALshootS;
    }

    /* Set pointer to path-slacks part of AL gradient at iterate */
    inline void setGradientPath(double* gALpath){
        this->gALpath = gALpath;
    }
    
    /* Set pointer to shooting part of AL gradient at candidate */
    inline void setGradientPathCandidate(double* gALpathS){
        this->gALpathS = gALpathS;
    }

    /* Set pointer to difference of path-constraint part of AL gradient */
    inline void setDiffGradientPath(double* dffGpath){
        this->dffGpath = dffGpath;
    }
    
    /* Set lower bound on input */
    void setUmin(double*);
    
    /* Set upper bound on input */
    void setUmax(double*);
    
    /* Set lower bound on state */
    void setXmin(double*);
    
    /* Set upper bound on state */
    void setXmax(double*);

    /* Set lower bound on terminal state */
    void setXminT(double*);

    /* Set upper bound on terminal state */
    void setXmaxT(double*);

    /**
    * Get methods
    */
    inline double getPredictionTime(){
        return Tpred;
    }
    
    inline int getDprim(){
        return dZ;
    }
    
    inline int getDprimShoot(){
    	return dZshoot;
    }
    
    inline int getDprimPath(){
    	return dZpath;
    }
    
    inline int getDdua(){
        return dL;
    }
    
    inline int getDduaShoot(){
        return dLshoot;
    }
    
    inline int getDduaPath(){
        return dLpath;
    }

    inline int getNs_(){
        return Ns_;
    }
    
    inline int getNs(){
        return Ns;
    }
    
    inline int getIzTerm(){
    	return izTerm;
    }
    
    inline int getDx(){
        return dx;
    }
    
    inline int getDu(){
        return du;
    }
    
    inline int getDy(){
   		return dy;
   	}
    
    inline int getDa(){
        return da;
    }
    
    inline int getDxu(){
        return dxu;
    }
    
    inline int getDw(){
        return dw;
    }
    
    inline int getDpc(){
        return dpc;
    }
    
    inline int getDper(){
        return dper;
    }

    inline double* getAugStates(){
        return augStates;
    }
    
    inline double* getUmin(){
    	return umin;
    }
    
    inline double* getUmax(){
    	return umax;
    }
    
    inline double* getXmin(){
    	return xmin;
    }
    
    inline double* getXmax(){
    	return xmax;
    }

    inline double* getXminT(){
        return xminT;
    }

    inline double* getXmaxT(){
        return xmaxT;
    }

    inline double* getIntegAugStates(){
        return integAugStates;
    }
    
    inline double* getInterpAugStates(){
        return interpAugStates;
    }

    inline double* getXest(){
        return xest;
    }

    inline int getIndexLastNode(){
        return izTerm;
    }

    /**
    ** Processing methods  
    **/
    /* Initialize members */
    void init();

    /* State variable: compute scaling transformation & adapt problem bounds */
    void computeStateScaling();
            
    /* Apply inverse state scaling D^-1*(.-c)*/
    void applyInvStateScaling();

    /* Apply inverse input scaling */
    void applyInvInputScaling(double*);

    /* Input variable: compute scaling transformation & adapt problem bounds */
    void computeInputScaling();

    /* Compute augmented states over horizon from iterate ==> val(z) & gALc(z) from */
    double shootStatesFromIterate(double);
    
    /* Compute augmented states over horizon from candidate ==> vcan(zS) */
    double shootStatesFromCandidate(double);
    
    /* Compute augmented adjoints over horizon from iterate ==> gAL(z) */
    void shootAdjointsFromIterate();

    /* Compute augmented adjoints over horizon from candidate ==> gAL(zS) */
    void shootAdjointsFromCandidate();

    /* Compute augmented states & adjoints from iterate ==> vAL(z), gAL(z) */
    double shootFromIterate(double);

    /* Compute augmented states and adjoints from candidate ==> vAL(zS), gAL(zS) */
    double shootFromCandidate(double);

    /* Compute augmented states and adjoints from iterate with periodicity constraint ==> vAL(z), gAL(z) */
    double shootFromIteratePeriodic(double);

    /* Compute augmented states and adjoints from candidate with periodicity constraint ==> vAL(zS), gAL(zS) */
    double shootFromCandidatePeriodic(double);

    /**
    ** Methods for debugging 
    **/
    /* Display multiple shooting parameters */
    void displayShootingParameters();

    /* Display infinity-norm of augmented states */
    void displayAsInfNorm();

    /* Display infinity-norm of augmented adjoints */
    void displayAaInfNorm();

    /* Display state scaling transformation */
    void displayStateTransfo();

    /* Display input scaling transformation */
    void displayInputTransfo();

private:
    /* Integrator (local instantiation) */
    ExplicitRKintegrator<AugOde> erkInteg;
 
    /* Dynamics, USER-DEFINED */
    dynT ode;
    
    /* Stage-cost, USER-DEFINED */
    costT cost;
    
    /* Path-constraint, USER-DEFINED */
    pconT pcon;
    
    /* Mayer term, USER-DEFINED */
    mayerT mayer;
    
    /**
    * Dimensions
    */
    /* State dimension */
    int dx;
    int dxbl;
    
    /* Input dimension */
    int du;
    int dubl;
    
    /* Output dimension */
    int dy;
    int dybl;
    
    /* Augmented state dimension (dx+1) */
    int dw;
    
    /* Augmented adjoint dimension (dw+du) */
    int da;
    
    /* Path-constraints dimension */
    int dpc;
    int dpcbl;
    
    /* Shooting variable (s_n',q_n')' dimension */
    int dxu;
    
    /* Dimension of separable AL gradient */
    int dxux;

    /* Dimension of periodicity constraint on state */
    int dper;
    
    /* Prediction time */
    double Tpred;
    
    /* Number of shooting nodes: s_0,s_1,...s_N */
    int Ns;
    /* # shooting intervals (Ns-1): [t_0,t_1], [t_1,t_2], ...,[t_{N-1},t_N] */
    int Ns_;
    
    /* Primal optimizer dimension
      (s_0',q_0',s_1',q_1',...,s_{N-1}',q_{N-1}',s_N',r_0',...r_{N-1}')' */
    int dZ;
    
    /* Shooting optimizer dimension 
       (s_0',q_0',s_1',q_1',...,s_{N-1}',q_{N-1}',s_N') */
    int dZshoot;
    
    /* Path-slacks optimizer dimension 
    	(r_0',...,r_{N-1}')' */
    int dZpath;
    
    /* Index of terminal state in optimizer z */
    int izTerm;
    
    /* Index of terminal state in partially separable gradient of augmented Lagrangian */
    int igTerm;
    
    /* Index of path-constraints slack in partially separable AL gradient */
    int igPath;
    
    /* Dual optimizer dimension */
    int dL;
    
    /* Shooting dual dimension */
    int dLshoot;
    
    /* Path dual dimension */
    int dLpath;
    
    /* Dimension of shooting part of full AL gradient */
    int dGshoot;
    int dGshootDBL;
    int dGpathDBL;

    /* Dimension of array of augmented states */
    int dAugS;
    
    /* Dimension of array of augmented adjoints */
    int dAugA;

    /* Dimension of local array of intermediate augmented states */
    int dwnsteps;

    /* Dimension of local array of interpolated augmented states */
    int sdwnsteps;

    /**
    * Locally allocated pointers
    */
    /* Augmented states (x(t;s_n,q_n)',z(t;s_n,q_n)')', n = 0...N-1, ALLOCATED */
    double* augStates;
    
    /* Augmented adjoints, ALLOCATED */
    double* augAdjoints;
    
    /* Output & input reference trajectory, ALLOCATED */
    double* yRef;
    double* uRef;
    
    /* Auxiliary vector for path constraints g(s_n,q_n)+r_n, ALLOCATED */
    double* econPath;
    
    /* Auxiliary vector for final augmented-adjoint variables \lambda_n, (\lambda_n^w,\lambda_n^u),
     n=0,...Ns-1, ALLOCATED */
    double* adfin;
    
    /* Lower bound on input, ALLOCATED */
    double* umin;
    
    /* Upper bound on input, ALLOCATED */
    double* umax;
    
    /* Lower bound on state, ALLOCATED */
    double* xmin;
    
    /* Upper bound on state, ALLOCATED */
    double* xmax;

    /* Lower bound on terminal state, ALLOCATED */
    double* xminT;

    /* Upper bound on terminal state, ALLOCATED */
    double* xmaxT;

    /* Diagonal state scaling matrix, ALLOCATED */
    double* Ds;

    /* Inverse diagonal state scaling matrix, ALLOCATED */
    double* iDs;

    /* Diagonal input scaling matrix, ALLOCATED */
    double* Dq;

    /* Inverse diagonal input scaling matrix, ALLOCATED */
    double* iDq;

    /* State centering (for scaling transfo that brings state within [0,1]), ALLOCATED */
    double* cs;

    /* Input centering (for scaling transfo that brings input within [0,1]), ALLOCATED */
    double* cq;
    
    /* Initial shooting augmented state, ALLOCATED */
    double* augStateIni;

    /* Scaled input, ALLOCATED */
    double* usc;

    /* Intermediate integration augmented-states on all shooting intervals, ALLOCATED */
    double* integAugStates; 

    /* Interpolated augmented-states on all shooting intervals, ALLOCATED */
    double* interpAugStates;
    
    /**
    * Pointers 
    */
    /* Pointer to initial state (parameter) */
    double* xest;

    /* Pointer to primal optimizer */
    double* z;
    
    /* Pointer to candidate (output of TPCG computed in SubspaceRefiner class) */
    double* zS;

    /* Pointer to dual optimizer */
    double* mu;
    
    /* Pointer to equality constraints (init.-value embedding + shooting constraint + slacked path-constraints) */
    double* econ;
    
    /* Pointer to full AL gradient wrt (s_0',q_0',s_1',q_1',...,s_N',r_0',...,r_{N-1}')' */
    double* gAL;
    
    /* Pointer to separated gradient of AL shooting part at iterate */
    double* gALshoot;

    /* Pointer to separated gradient of AL shooting part at candidate */
    double* gALshootS;
    
    /* Pointer to path-slack part of full AL gradient at iterate */
    double* gALpath;

    /* Pointer to path-slack part of AL gradient at candidate */
    double* gALpathS;

    /* Pointer to difference of gradient of AL path-constraint part */
    double* dffGpath;
    
};
#include "MuShoot.tpp"

#endif
