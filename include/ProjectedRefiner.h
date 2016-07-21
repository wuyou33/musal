#ifndef MusAL_ProjectedRefiner_h
#define MusAL_ProjectedRefiner_h

#include <iostream>
#include <assert.h>

#include "macros.h"
#include "QuadraticModel.h"
#include "Preconditioner.h"

/*****
** Implements Moré-Toraldo projected searchs in the refinement step 
***/
class ProjectedRefiner {

public: 
	/* Default constructor */
	ProjectedRefiner();

	/* Destructor */
	~ProjectedRefiner();

	/* Copy constructor */
	ProjectedRefiner(const ProjectedRefiner &);

	/* Assignment operator */
	ProjectedRefiner& operator= (const ProjectedRefiner&);

	/**
	** Set methods 
	**/
	/* Set quadratic model */
    inline void setQuadraticModel(const QuadraticModel& qmod){
        this->qmod = qmod;
    }

    /* Set preconditioner */
    inline void setPreconditioner(const Preconditioner& precon){
        this->precon = precon;
    }
    
    /* Set pointer to full gradient */
    inline void setGradient(double* gAL){
        this->gAL = gAL;
    }
    
    /* Set pointer to product against AL hessian approx */
    inline void setProductHessian(double* prdHss){
        this->prdHss = prdHss;
    }

    /* Set pointers to hessians */
    inline void setHessians(double* hALshoot,double* hALpath){
        this->hALshoot = hALshoot;
        this->hALpath = hALpath;
    }

    /* Set pointers to reduced hessians */
    inline void setReducedHessians(double* hALshootRed,double* hALpathRed){
        this->hALshootRed = hALshootRed;
        this->hALpathRed = hALpathRed;
    }
    
    /* Set shooting dimensions */
    inline void setDimensionShoot(int dx,int du,int dZshoot,int dHshoot,int dhS,int dGshoot){
        this->dZshoot = dZshoot;
        dZshootBL = dZshoot*SZDBL;
        this->dHshoot = dHshoot;
        this->dGshoot = dGshoot;
        this->dhS = dhS;
        this->dx = dx;
        this->du = du;
        dxu = dx+du;
        dxux = dxu+dx;
        precon.setShootDim(dZshoot,dxux,dhS,dHshoot,dGshoot);
    }
    
    /* Set path-constraint dimensions */
    inline void setDimensionPath(int dpc,int dZpath,int dHpath,int dhP,int dGpath){ 
        this->dZpath = dZpath;
        dZpathBL = dZpath*SZDBL; 
        this->dHpath = dHpath;
        this->dGpath = dGpath;
        this->dhP = dhP;
        this->dpc = dpc;
        precon.setSlackDim(dZpath,dpc,dhP,dHpath,dGpath);
    }
    
    /* Set some other useful dimensions */
    inline void setDimension(int dZ){
        this->dZ = dZ;
        dZBL = dZ*SZDBL;
    }

    /* Set pointer to candidate */
    inline void setMinorIterate(double* zminor){
        this->zminor = zminor;
    }
        
    /* Set pointer to current trust region iterate */
    inline void setIterate(double* z){
        this->z = z;
    }
    
    /* Set # free variables */
    inline void setNumFree(int nfree,int nfreeS,int nfreeP){
        this->dimRed = nfree;
        this->dimRedBL = dimRed*SZDBL;
        this->dimRedShoot = nfreeS;
        this->dimRedPath = nfreeP;
        precon.setDimRed(nfree,nfreeS,nfreeP);
    }

    /* Set pointers to NLP bounds */
    inline void setNLPbounds(double* zmin,double* zmax){
        this->zmin = zmin;
        this->zmax = zmax;
    }

    /* Set # shooting nodes */
    inline void setNumIntervals(int Ns_){
        this->Ns_ = Ns_;
        this->Ns_BL = Ns_*SZDBL;
        this->NsBL = (Ns_+1)*SZDBL;
        this->Ns_NT = Ns_*SZINT;
        this->NsNT = (Ns_+1)*SZINT;
    }
        
    /* Set regularization coefficient for preconditionner */
    inline void setRegCoefPrec(double pr_reg_coef){
        this->pr_reg_coef = pr_reg_coef;
        precon.setDiagRegCoef(pr_reg_coef);
    }

    /* Set tolerance on activity */
    inline void setActivityTolerance(double acTol){ 
        this->activTol = acTol;
    }

    /* Set pointer to free variables */
    inline void setFreeVars(int* freeVars){
        this->freeVars = freeVars;
    }

    /* Set pointer to indices of free variables */
    inline void setIndFreeVars(int* indFreeVars){
        this->indFreeVars = indFreeVars;
    }

    /* Set pointer to status of each variable */
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

    /* Set preconditioner bandwidth */
    inline void setBandWidth(int bwd){
        this->bwd = bwd;
        precon.setBandWidth(bwd);
    }

    /* Set pivot tolerance for LDL' on band preconditioner */
    inline void setPivotTolerance(double pivTol){
        precon.setPivotTolerance(pivTol);
    }

    /* Set backtracking coefficient for projected search */
    inline void setBetaSrch(double betaSrch){
        this->betaSrch = betaSrch;
    }

    /* Set constant for sufficient decrease */
    inline void setCstSrch(double cstSrch){
        this->cstSrch = cstSrch;
    }

    /* Set max # refinement iterations */
    inline void setMaxRefIter(int maxRefIter){
        this->maxRefIter = maxRefIter;
    }

    /* Set max # backtracking iter in projected search */
    inline void setMaxItSrch(int maxItSrch){
        this->maxItSrch = maxItSrch;
    }

    /* Set model value */
    inline void setValMod(double vmod){
        this->vmod = vmod;
    }

    /**
    ** Gets
    **/
    /* Outputs dimension of subspace of free variables at current minor iteration */
    inline int getDimRed(){
        return dimRed;
    }

    /* Outputs bandwidth */
    inline int getBandWidth(){
        return precon.getBandWidth();
    }

    /* Outputs pointer to reduced residual */
    inline double* getReducedResidual(){
        return rRed;
    }

    /* Outputs total number of CG iterations */
    inline int getCumCG(){
        return cumCG;
    }

    /* Outputs PCG status */
    inline int getPCGstatus(){
        return pcgStatus;
    }

    /* Outputs model value at current minor iterate */
    inline double getValMod(){
        return vmod;
    }

	/**
	** Numerical processing methods
	**/
	/* Refinement with projected searches along directions
	   generated by PCG iterations */
	void findCandidate(const double,const double,const double);

	/* Generate descent direction by PCG iterations */
	void genDesDir(const double,const double);

	/* Projected search along descent direction */
	void projSearch();

    /* Extract active-set at minor iterate */
    void extractActiveSet(double*,double*);

    /* Assemble reduced hessian in a block-wise separated fashion */
    void assembleReducedHessian();

    /**
    ** Other methods
    **/
    /* Initialize and allocate local data */
    void init();

    /**
    ** For debugging
    **/
    /* Displays shooting part of AL hessian */
    void displayHalShoot();

    /* Displays residuals rRed & yRed */
    void displayResiduals();

    /* Displays conjugate direction pRed */
    void displayConjugateDirection();

    /* Displays reduced descent direction wRed */
    void displayReducedDesDir();

    /* Displays reduced product xRed */
    void displayReducedProduct();

    /* Evaluate model q(z+w) where w is generated during PCG iter */
    void checkPCGiter(double*);

private: 
	/* Quadratic model */
	QuadraticModel qmod;

	/* Preconditioner */
	Preconditioner precon;

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
    
    /* Dimension of shooting AL gradient */
    int dGshoot;

    /* Dimension of shooting AL hessian */
    int dHshoot;
    
    /* Dimension of hessian shooting block */
    int dhS;

    /* Dimension of slack AL gradient */
    int dGpath;

    /* Dimension of slack AL hessian */
    int dHpath;
    
    /* Dimension of hessian slack block */
    int dhP;

    /* # shooting intervals */
    int Ns_;
    int Ns_BL;
    int NsBL;
    int Ns_NT;
    int NsNT;

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
    
    /* Dimension of reduced subspace */
    int dimRed;
    int dimRedBL;

    /* Dimension of shooting reduced subspace */
    int dimRedShoot;

    /* Dimension of slack reduced subspace */
    int dimRedPath;

    /**
    ** Other members
    **/
    /* Quadratic model value at minor iterate */
    double vmod;

    /* Total number of CG iterations */
    int cumCG;

    /**
    ** Pointers to external data
    **/
    /* Pointer to current iterate */
    double* z;
    
    /* Pointer to candidate point */
    double* zminor;

    /* Auxiliary working array for projected search, ALLOCATED */
    double* zwork;

    /* Pointers to NLP bounds */
    double *zmin;
    double *zmax;

    /* Pointer to AL gradient */
    double* gAL;
    
    /* Pointer to product with AL hessian approximation B(z)*(zMinorCur-z) */
    double* prdHss;

    /* Auxiliary working array for projected search, ALLOCATED */
    double* prdHssWrk;

    /* Pointers to shooting hessian and path-constraints hessian */
    double* hALshoot;
    double* hALpath;

    /* Pointers to reduced shooting hessian and slack hessian*/
    double* hALshootRed;
    double* hALpathRed;

    /* Pointer to free variables */
    int* freeVars;

    /* Pointer to indices of free variables */
    int* indFreeVars;

    /* Pointer to status of each primal variable */
    int* varStatus;

    /* Pointers to indices of free vars in states, inputs, slacks and groups */
    int* ivarFreeState;
    int* ivarFreeInput;
    int* ivarFreeSlack;
    int* ivarFreeGroupShoot;

    /* Pointers to # free vars per state, input, shooting group or slack */
    int* nfreeState;
    int* nfreeInput;
    int* nfreeGroupShoot;
    int* nfreeSlack;

    /**
    ** PCG parameters 
    **/
    /* PCG status
    * - -1: negative curvature direction or iterate goes outside trust region
    * - 1: iterated till stopping criterion met
    * - 2: reached max. # PCG iterations
    */
    int pcgStatus;
    
    /* Regularization coefficient for making Jacobi preconditionner positive definite */
    double pr_reg_coef;

    /* Bandwidth of band preconditioner */
    int bwd;

    /* Tolerance on activity */
    double activTol;

    /**
    ** PCG local allocated data
    **/
    /* Directions for subspace search, ALLOCATED */
    double* wRed; // Reduced descent direction
    double* rRed; // Reduced residual
    double* yRed; // Preconditionned reduced residual
    double* pRed; // Conjugate direction
    double* xRed; // Curvature direction (xRed <- Hred*pRed)

    /**
    ** Options for projected search
    **/
    /* Step-size */
    double a;

    /* Backtracking multiplicative coefficient (<1) */
    double betaSrch;

    /* Constant for sufficient decrease */
    double cstSrch;

    /* Maximum # backtracking iterations */
    int maxItSrch;

    /* Max # refinement iterations (proj. search) */
    int maxRefIter;
};


#endif