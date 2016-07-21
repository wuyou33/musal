#ifndef MusAL_PreconRefine_h
#define MusAL_PreconRefine_h

#include <iostream>
#include <assert.h>

#include "macros.h"
#include "QuadraticModel.h"
#include "Preconditioner.h"

/**
** Implements preconditionned refinement step
**/
class PreconRefine {

public:
	/* Default preconditionner */
	PreconRefine();

	/* Destructor */
	~PreconRefine();

	/* Copy constructor */
	PreconRefine(const PreconRefine&);

	/* Assignment operator */
	PreconRefine& operator= (const PreconRefine&);

    /**
    ** Set methods
    **/
    /* Set quadratic model */
    inline void setQuadraticModel(const QuadraticModel& qm){
        qmod = qm;
    }

    /* Set preconditioner */
    inline void setPreconditioner(const Preconditioner& precon_){
        precon = precon_;
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
    inline void setHessians(double* hessShoot,double* hessPath){
        hALshoot = hessShoot;
        hALpath = hessPath;
    }

    /* Set pointers to reduced hessians */
    inline void setReducedHessians(double* hessShootRed,double* hessPathRed){
        hALshootRed = hessShootRed;
        hALpathRed = hessPathRed;
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
    inline void setCandidate(double* zNex){
        this->zNex = zNex;
    }
    
    /* Set pointer to current iterate */
    inline void setIterate(double* z){
        this->z = z;
    }
    
    /* Set # free variables */
    inline void setNumFree(int nfree,int nfreeS,int nfreeP){
        this->dimRed = nfree;
        this->dimRedBL = dimRed*SZDBL;
        this->dimRedShoot = nfreeS;
        this->dimRedPath = nfreeP;
        this->maxIt = nfree;
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
    
    /* Set max # TPCG iterations (dimension of free subspace) */
    inline void setMaxIt(int maxIt){
        this->maxIt = maxIt;
    }
        
    /* Set max # cumulative PCG iter */
    inline void setMaxCumCG(const int maxCumCG){
        this->maxCumCG = maxCumCG;
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

    /**
    ** Gets
    **/
    /* Outputs dimension of subspace of free variables */
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

    /* Outputs total number of restarts */
    inline int getCumRestart(){
        return cumRestart;
    }

	/**
	** Numerical processing methods
	**/
	/* Subspace refinement (stopped when tolerance reached) */
	void findCandidate();

    /* Subspace refinement with restart if problem bound hit */
    void findCandidateRestart(const double,const double);

    /* Subspace refinement with restart and real-time refinement */
    void findCandidateRestartRT(const double,const double,const int);

    /* Extract active-set at Cauchy point */
    void extractActiveSet();

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

    /* Displays reduced optimizer zRed */
    void displayReducedOptimizer();

    /* Displays reduced product xRed */
    void displayReducedProduct();

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
    /* Total number of CG iterations (in one PCG loop) */
    int cumCG;

    /* Total number of restarts (in one PCG loop) */
    int cumRestart;

    /**
    ** Pointers to external data
    **/
    /* Pointer to current iterate */
    double* z;
    
    /* Pointer to candidate for next iterate */
    double* zNex;
    
    /* Pointers to NLP bounds */
    double *zmin;
    double *zmax;

    /* Pointer to AL gradient */
    double* gAL;
    
    /* Pointer to product with AL hessian approximation */
    double* prdHss;

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

    /* Pointers to # of free vars per state, input, slack and group */
    int* nfreeState;
    int* nfreeInput;
    int* nfreeSlack;
    int* nfreeGroupShoot;

    /**
    ** TPCG options 
    **/
    /* For RT: max. # cumulative PCG iter in primal loop */
    int maxCumCG;

    /* Max # TPCG iterations */
    int maxIt;

    /* Regularization coefficient for making Jacobi preconditionner positive definite */
    double pr_reg_coef;

    /* Bandwidth of band preconditioner */
    int bwd;

    /* Tolerance on activity */
    double activTol;

    /**
    ** TPCG local allocated data
    **/
    /* Directions for subspace search, ALLOCATED */
    double* zRed; // Subspace optimizer
    double* rRed; // Reduced residual
    double* yRed; // Preconditionned reduced residual
    double* pRed; // Conjugate direction
    double* xRed; // Curvature direction (xRed <- Hred*pRed)

    /* Reduced NLP bounds, ALLOCATED */
    double* zRedmin;
    double* zRedmax;
};


#endif MusAL_PreconRefine_h