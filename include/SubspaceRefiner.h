//
//  SubspaceRefiner.h
//  MusAL
//
//  Created by Jean on 4/9/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_SubspaceRefiner_h
#define MusAL_SubspaceRefiner_h

#include <iostream>
#include <assert.h>

#include "macros.h"
#include "QuadraticModel.h"

/**
* Implements subspace refinement on null-space of active constraints
*/
class SubspaceRefiner {

public:
    /* Default construtor */
    SubspaceRefiner();

    /* Copy constructor */
    SubspaceRefiner(const SubspaceRefiner&);
    
    /* Assigment operator */
    SubspaceRefiner& operator= (const SubspaceRefiner&);
    
    /* Destructor */
    ~SubspaceRefiner();
    
    /**
    ** Sets
    **/
    /* Set quadratic model */
    inline void setQuadraticModel(const QuadraticModel& qm){
        qmod = qm;
    }
    
    /* Set pointer to full gradient */
    inline void setGradient(double* g){
        gAL = g;
    }
    
    /* Set pointers to hessians */
    inline void setHessians(double* hessShoot,double* hessPath){
        hALshoot = hessShoot;
        hALpath = hessPath;
    }
    
    /* Set shooting dimensions */
    inline void setDimensionShoot(int dx,int du,int dZshoot,int dHshoot){
        this->dZshoot = dZshoot;
        dZshootBL = dZshoot*SZDBL;
        this->dHshoot = dHshoot;
        this->dx = dx;
        this->du = du;
        dxu = dx+du;
        dxux = dxu+dx;
    }
    
    /* Set path-constraint dimensions */
    inline void setDimensionPath(int dpc,int dZpath,int dHpath){
        this->dZpath = dZpath;
        dZpathBL = dZpath*SZDBL;
        this->dHpath = dHpath;
        this->dpc = dpc;
    }
    
    /* Set some other useful dimensions */
    inline void setDimension(int dZ){
        this->dZ = dZ;
        dZBL = dZ*SZDBL;
    }

    /* Set pointer to candidate */
    inline void setCandidate(double* zN){
        zNex = zN;
    }
    
    /* Set pointer to Cauchy point */
    inline void setCauchy(double* zC){
        this->zC = zC;
    }
    
    /* Set pointer to current iterate */
    inline void setIterate(double* z){
        this->z = z;
    }
    
    /* Set # free variables */
    inline void setNumFree(int nfree_){
    	nfree = nfree_;
    	dimRed = nfree_;
    	dimRedBL = dimRed*SZDBL;
        maxIt = nfree_;
    }
    
    /* Set pointer to number of free variables per group */
    inline void setNumFreeGroup(){
        
    }

    /* Set pointer to free/active variables  */
    inline void setFreeVars(int* freeVars){
        this->freeVars = freeVars;
    }
    
    /* Set pointers to full bounds (NLP bounds) */
    inline void setNLPbounds(double* zmin,double* zmax){
        this->zmin = zmin;   
        this->zmax = zmax;
    }

    /* Set # shooting nodes */
    inline void setNumIntervals(int Ns_){
        this->Ns_ = Ns_;
    }
    
    /* Set max # TPCG iterations (dimension of free subspace) */
    inline void setMaxIt(int maxIt){
        this->maxIt = maxIt;
    }
    
    /* Set squared tolerance on norm square of CG residual */
    inline void setTol2(double tol2){
        this->tol2 = tol2;
    }
    
    /* Set regularization coefficient for preconditionner */
    inline void setRegCoefPrec(double pr_reg_coef){
        this->pr_reg_coef = pr_reg_coef;
    }

    /* Set tolerance on activity */
    inline void setActivityTolerance(double acTol){ 
        this->activTol = acTol;
    }

    /**
     ** Gets
     **/
    /* Outputs number of free variables */
    inline int getNumFree(){
        return nfree;
    }
    
    /* Outputs pointer to free/active variables */
    inline int* getFreeVars(){
        return freeVars;
    }
    
    /**
    ** Numerical processing methods 
    **/ 
    /* Subspace refinement (stopped when tolerance reached) */
    void findCandidate(const double);

    /* Subspace refinement with restart if problem constraint hit */
    void findCandidateRestart(const double);

    /* Real-time subspace refinement (stopped after fixed # iterations) */
    void findCandidateRealTime(const double);
    
    /* Evaluate model at current reduced CG iterate */
    double evalModelReducedCG(double*);

    /**
    ** Other methods
    **/
    /* Initialize and allocate local data */
    void init();

    /**
    ** Displaying methods
    **/
    /* Display shooting part of AL hessian */
    void displayHalShoot();
    
    /* Display diagonal of shooting part of AL hessian */
    void displayDiagHalShoot();

    /* Display diagonal preconditionner */
    void displayDiagPrec();

private:
    /* Quadratic model */
    QuadraticModel qmod;
    
    /**
    ** Dimensions
    **/
    /* Full primal dimension */
    int dZ;
    int dZBL;
    
    /* Dimension of shooting variable */
    int dZshoot;
    int dZshootBL;
    
    /* Dimension of path variable */
    int dZpath;
    int dZpathBL;
    
    /* Dimension of shooting AL hessian */
    int dHshoot;
    
    /* Dimension of path-constraint AL hessian */
    int dHpath;
    
    /* # shooting intervals */
    int Ns_;
    
    /* State dimension */
    int dx;
    
    /* Input dimension */
    int du;
    
    /* State+input dimension */
    int dxu;
    
    /* State+input+state dimension */
    int dxux;
    
    /* Path-constraint dimension */
    int dpc;
    
    /* # free variables */
    int nfree;
    
    /* Dimension of reduced subspace */
    int dimRed;
    int dimRedBL;

    /**
    ** Pointers to external data
    **/
    /* Pointer to current iterate */
    double* z;
    
    /* Pointer to candidate for next iterate */
    double* zNex;
    
    /* Pointer to Cauchy point */
    double* zC;
    
    /* Pointer to free/active variables */
    int* freeVars;

    /* Pointer to number of free variables per group */
    int* nfreeGroup;
    
    /* Pointer to AL gradient */
    double* gAL;
    
    /* Pointers to shooting hessian and path-constraints hessian */
    double* hALshoot;
    double* hALpath;
    
    /* Pointers to full NLP bounds */
    double *zmin;
    double *zmax;
    
    /**
    ** TPCG options 
    **/
    /* Max # TPCG iterations */
    int maxIt;
    
    /* Termination tolerance on TPCG iterations */
    double tol2;

    /* Regularization coefficient for making preconditionner positive definite */
    double pr_reg_coef;

    /* Tolerance on activity */
    double activTol;

    /**
    ** TPCG local allocated data
    **/
    /* Diagonal of AL hessian (full dimension), ALLOCATED */
    double* diagHal;
    
    /* Diagonal pre-conditionner for reduced AL hessian (reduced dimension), ALLOCATED */
    double* diagPrec;

	/* Auxiliary vector for H*(z-zC), ALLOCATED */
	double* auxFul;

    /* Directions for subspace search, ALLOCATED */
    double* zRed; // Subspace optimizer
    double* rRed; // Reduced residual
    double* yRed; // Preconditionned reduced residual
    double* pRed; // Conjugate direction
    double* xRed; // Curvature direction (xRed <- Hred*pRed)
    
    /* Auxiliary full vectors for product with reduced hessian, ALLOCATED */
    double* vecFulIn;
    double* vecFulOut;

    /* Auxiliary full-space vector to evaluate model, ALLOCATED */
    double* zeval;

    /* Reduced NLP bounds, ALLOCATED */
    double* zRedmin;
    double* zRedmax;
};


#endif
