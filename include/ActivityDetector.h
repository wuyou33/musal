//
//  ActivityDetector.h
//  MusAL
//
//  Created by Jean on 4/9/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_ActivityDetector_h
#define MusAL_ActivityDetector_h

#include <iostream>
#include <assert.h>

#include "QuadraticModel.h"
#include "macros.h"

/*******
* Implements an activity detection mechanism taken from one of the following:
*   - Gradient projection with generalized Armijo condition
*   - Proximal gradient step
*   - Lancelot's cauchy search
*   - Proximal Alternating Linearized Minimization
***/
class ActivityDetector {

public:
    /* Default constructor */
    ActivityDetector();
    
    /* Destructor */
    ~ActivityDetector();
    
    /* Copy constructor */
    ActivityDetector(const ActivityDetector&);
    
    /* Assignment operator */
    ActivityDetector& operator= (const ActivityDetector&);

    /**
    ** Set methods
    **/    
    /* Set quadratic model */
    inline void setQuadraticModel(const QuadraticModel& qmod){
        this->qmod = qmod;
    }
    
    /* Set dimension */
    inline void setDimensions(int dZ,int dZshoot,int dGshoot){
        this->dZ = dZ;
        this->dZshoot = dZshoot;
        this->dGshoot = dGshoot;
    }

    /* Set state dimension */
    inline void setStateDimension(int dx){
        assert(dx>=1);
        this->dx = dx;
        this->dx_ = dx-1; 
    }

    /* Set input dimension */
    inline void setInputDimension(int du){
        assert(du>=1);
        this->du = du;
    }

    /* Set groups dimensions */
    inline void setGroupDimensions(int dxu,int dxux,int dpc){
        this->dxu = dxu;
        this->dxu_ = dxu-1; 
        this->dxux = dxux;
        this->dpc = dpc;
    }
    
    /* Set lower and upper bounds */
    inline void setNLPbounds(double* zL,double* zU){
        this->zL = zL;
        this->zU = zU;
    }
    
    /* Set iterate */
    inline void setIterate(double* z){
        this->z = z;
    }
    
    /* Set Cauchy point */
    inline void setCauchy(double* zC){
        this->zC = zC;
    }
    
    /* Set gradients */
    inline void setGradient(double* grad){
    	this->grad = grad;
    }
    
    /* Set backtracking coefficient */
    inline void setBacktrackCoef(double beta){
        assert(beta>0);
        assert(beta<1);
        betaBck = beta;
    }

    /* Set backtracking constant */
    inline void setBacktrackConst(double c){
        assert(c>0);
        assert(c<1);
        cstBck = c;
    }
    
    /* Set maximum # backtracking iterations */
    inline void setBacktrackMaxIter(int maxIt){
        maxItBck = maxIt;
    }
    
    /* Set interpolation coefficient */
    inline void setBetaInter(double betaInter){
        this->betaInter = betaInter;
    }

    /* Set extrapolation coefficient */
    inline void setBetaExtra(double betaExtra){
        this->betaExtra = betaExtra;
    }

    /* Set max interpolation iter */
    inline void setMaxInterIter(int maxInterIter){
        this->maxInterIter = maxInterIter;
    }

    /* Set max extrapolation iter */
    inline void setMaxExtraIter(int maxExtraIter){
        this->maxExtraIter = maxExtraIter;
    }

    /* Set tolerance on activity */
    inline void setActivityTolerance(double acTol){
        activTol = acTol;
    }

    /* Set tolerance on inf-norm of primal difference */
    inline void setDiffTolerance(double diffTol){
        this->diffTol = diffTol; 
    }

    /* Set # shooting intervals */
    inline void setNs_(int Ns_){
        this->Ns_ = Ns_;
        this->Ns_BL = Ns_*SZDBL;
        this->NsBL = (Ns_+1)*SZDBL;
        this->Ns_NT = Ns_*SZINT;
        this->NsNT = (Ns_+1)*SZINT;
    }

    /* Set pointer to reduced residual */
    inline void setReducedResidual(double* rRed){
        this->rRed = rRed;
    }

    /* Set pointer to product prdHss */
    inline void setProductHessian(double* prdHss){
        this->prdHss = prdHss;
    }  

    /* Set pointer to free variables */
    inline void setFreeVars(int* freeVars){
        this->freeVars = freeVars;
    }

    /* Set pointer to indices of free variables */
    inline void setIndFreeVars(int* indFreeVars){
        this->indFreeVars = indFreeVars;
    }

    /* Set pointer to variables status */
    inline void setVarStatus(int* varStatus){
        this->varStatus = varStatus;
    }

    /* Set pointers to indices of free variables per state, input, slack and group */
    inline void setIndexFreeVars(int* ivarFreeState,int* ivarFreeInput,int* ivarFreeSlack,int* ivarFreeGroupShoot){
        this->ivarFreeState = ivarFreeState;
        this->ivarFreeInput = ivarFreeInput;
        this->ivarFreeSlack = ivarFreeSlack;
        this->ivarFreeGroupShoot = ivarFreeGroupShoot;
    }

    /* Set pointers to number of free variables per state, input, slack and group */
    inline void setNumberFreeVars(int* nfreeState,int* nfreeInput,int* nfreeSlack,int* nfreeGroupShoot){
        this->nfreeState = nfreeState;
        this->nfreeInput = nfreeInput;
        this->nfreeSlack = nfreeSlack;
        this->nfreeGroupShoot = nfreeGroupShoot;
    }

    /**
    ** Get methods
    **/
    /* Returns error indicator */
    inline int getError(){
        return err;
    }

    /* Returns # free variables */
    inline int getNumberFree(){
        return nfree;
    }

    /* Returns # free shooting variables */
    inline int getNumberFreeShoot(){
        return nfreeS;
    }
    
    /* Returns # free slack variables */
    inline int getNumberFreeSlack(){
        return nfreeP;
    }
    
    /**
    ** Numerical processing methods 
    **/
    /* Activity detection via projected-search (Armijo-like backtracking) */
    void findProjectedSearch(double*,double*,double*,const double);

    /* Activity detection via interpolation & extrapolation procedure (see TRON) */
    void findInterpExtrap(double*,double*,double*,const double);

    /* Active-set extraction at Cauchy point */
    void extractActiveSet(double*,double*);
    
    /**
    ** Methods for debugging
    **/
    /* Check model decrease along projected arc */
    void checkDescent();

    /* Allocate auxiliary memory space for debugging */
    void allocateDebug();

    /* Display Cauchy point */
    void displayCauchy();

    /* Display infinity-norm of Cauchy point */
    void displayInfNormCauchy();

    /* Display infinity-norm of model gradient */
    void displayInfNormGrad();

    /* Displays prdHss */
    void displayProductHessian();

    /* Displays gradient */
    void displayGradient();

    /* Displays # free states and indices */
    void displayFreeStates();

    /* Displays # free inputs and indices */
    void displayFreeInputs();

    /* Displays # free groups and indices */
    void displayFreeGroups();

    /* Displays indices of free variables */
    void displayFreeIndices();

private:
    /* Quadratic model */
    QuadraticModel qmod;

    /* Auxiliary vector for debugging, ALLOCATED */
    double* zdbg;

    /* Auxiliary vector for debugging, ALLOCATED */
    double* prdHssDbg;

    /**
    ** Dimensions 
    **/
    /* Optimizer dimension */
    int dZ;
    int dZshoot;

    /* Gradient dimensions */
    int dGshoot;

    /* # shooting intervals */
    int Ns_;
    int Ns_BL;
    int NsBL;
    int Ns_NT;
    int NsNT;

    /* State dimension */
    int dx;
    int dx_;

    /* Input dimension */
    int du;

    /* State + input dimension */
    int dxu;
    int dxu_;

    /* Shooting group dimension */
    int dxux;

    /* Path-slack dimension */
    int dpc;

    /**
    ** Parameters for Cauchy point search
    **/
    /* Backtracking multiplicative coefficient (<1) */
    double betaBck;

    /* Backtracking constant (sufficient decrease) */
    double cstBck;
    
    /* Maximum # backtracking iterations */
    int maxItBck;

    /* Interpolation coefficient */
    double betaInter;

    /* Extrapolation coefficient */
    double betaExtra;

    /* Max # interpolation iter */
    int maxInterIter;

    /* Max # extrapolation iter */
    int maxExtraIter;

    /* Tolerance on activity */
    double activTol;

    /* Tolerance on infinity norm of difference between Cauchy point and iterate */
    double diffTol;

    /**
    ** Updated variables during activity detection
    **/
    /* Error indicator: 1 if error in projected search, 0 otherwise */
    int err;

    /* Number of free variables */
    int nfree;

    /* Number of free shooting variables (s or q)*/
    int nfreeS;

    /* Number of free slack variables (r) */
    int nfreeP;

    /**
    ** Pointers to external data 
    **/
    /* Pointers to upper and lower bounds of NLP */
    double* zL;
    double* zU;
    
    /* Pointer to Cauchy point */
    double* zC;
    
    /* Pointer to primal iterate z */
    double* z;
    
    /* Pointer to full AL gradient */
    double* grad;
    
    /* Pointer to set of free/active variables */
    int* freeVars;

    /* Pointer to indices of free variables */
    int* indFreeVars;

    /* Pointer to variables status */
    int* varStatus;

    /* Pointers to indices of free variables in state, input, slack and shooting group */
    int* ivarFreeState;
    int* ivarFreeInput;
    int* ivarFreeSlack;
    int* ivarFreeGroupShoot;

    /* Pointers to number of free variables in state, input, slack and shooting group */
    int* nfreeState;
    int* nfreeInput;
    int* nfreeSlack;    
    int* nfreeGroupShoot;

    /* Pointer to reduced residual */
    double* rRed;

    /* Pointer to product B*(z-x) */
    double* prdHss;
};
    

#endif
