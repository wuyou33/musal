#ifndef PRECOND_H
#define PRECOND_H

#include <iostream>
#include <assert.h>

#include "macros.h"

/**
** Implements preconditioner used in PreconRefine class 
**/
class Preconditioner {

public:
	/* Default constructor */
	Preconditioner();

	/* Destructor */
	~Preconditioner();

	/* Copy constructor */
	Preconditioner(const Preconditioner&);

	/* Assignment operator */
	Preconditioner& operator= (const Preconditioner&);

	/**
	** Set methods
	**/
	/* Set pointer to reduced hessian approx */
	inline void setReducedHessians(double* hALshootRed,double* hALpathRed){
		this->hALshootRed = hALshootRed;
		this->hALpathRed = hALpathRed;
	}

	/* Set full-space dimension */
	inline void setDimFull(int dZ){
		this->dZ = dZ;
	}

	/* Set reduced-space dimension */
	inline void setDimRed(int nfree,int nfreeS,int nfreeP){
		dRed = nfree;
		dRedBL = dRed*SZDBL;
		dRedShoot = nfreeS;
		dRedShootBL = dRedShoot*SZDBL;
		dRedPath = nfreeP;
		dRedPathBL = dRedPath*SZDBL;
	}	

	/* Set number of shooting intervals */
	inline void setNs_(int Ns_){
		this->Ns_ = Ns_;
	}

	/* Set shooting-block dimension */
	inline void setShootDim(int dZshoot,int dxux,int dhS,int dHshoot,int dGshoot){
		this->dxux = dxux;
		this->dhS = dhS;
		this->dHshoot = dHshoot;
		this->dHshootBL = dHshoot*SZDBL;
		this->dGshoot = dGshoot;
		this->dZshoot = dZshoot;
	}

	/* Set slack-block dimension */
	inline void setSlackDim(int dZpath,int dpc,int dhP,int dHpath,int dGpath){
		this->dpc = dpc;
		this->dhP = dhP;
		this->dHpath = dHpath;
		this->dHpathBL = dHpath*SZDBL;
		this->dGpath = dGpath;
		this->dZpath = dZpath;
	}

	/* Sets regularization coef for diagonal preconditioners */
	inline void setDiagRegCoef(double jacRegCoef){
		this->jacRegCoef = jacRegCoef;
	}

	/* Set pivot tolerance for factorizing preconditioner */
	inline void setPivotTolerance(double pivTol){
		this->pivTol = pivTol;
	}

	/* Set bandwidth for band preconditioner */
	inline void setBandWidth(int bw){
		this->bw = bw;
		this->bbw = bw+1;
	}

	/* Set pointer to free variables */
    inline void setFreeVars(int* freeVars){
        this->freeVars = freeVars;
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
	/* Returns pointer to lower triangular Cholesky factor L (shooting part) */
	inline double* getBandShoot(){
		return pMbShoot;
	}

	/* Returns pointer to Jacobi (diagonal) preconditioner D */
	inline double* getPrecD(){
		return precD;
	}

	/* Returns reduced-space dimension */
	inline int getDimRed(){
		return dRed;
	}

	/* Returns bandwidth */
	inline int getBandWidth(){
		return bw;
	}

	/**
	** Initialization method
	**/
	void init();

	/**
	** Numerical processing methods
	**/
	/* Builds Jacobi (diagonal) preconditioner */
	void buildJacobi();

	/* Applies Jacobi (diagonal) preconditioner */
	void applyJacobi(double*,double*); 

	/* Assemble and factorize banded preconditioner */
	void buildBand();

	/* Apply banded preconditioner */
	void applyBand(double*,double*);

	/**
	** For debugging
	**/
	/* Displays Jacobi preconditioner */
	void displayJacobi();

	/* Displays reduced hessian */
	void displayReducedHessian();

	/* Displays band preconditioner */
	void displayBand();

private:
	/* Number of shooting intervals */
	int Ns_;

	/* Shooting block dimension */
	int dxux;

	/* Path block dimension */
	int dpc;

	/* Full-space dimension */
	int dZ;

	/* Full-space dimension (shooting) */
	int dZshoot;

	/* Full-space dimension (path) */
	int dZpath;

	/* Dimension of null-space of active constraints */
	int dRed;
	int dRedBL;

	/* Dimension of shooting part of null-space of active constraints */
	int dRedShoot;
	int dRedShootBL;

	/* Dimension of slack part of null-space of active constraints */
	int dRedPath;
	int dRedPathBL;

	/* Dimension of separated groupwise shooting AL gradient */
	int dGshoot;

	/* Dimension of shooting hessian approx */
	int dHshoot;
	int dHshootBL;

	/* Dimension of a full shooting hessian block */
	int dhS;

	/* Dimension of separated groupwise slack AL gradient */
	int dGpath;

	/* Dimension of slack hessian approx */
	int dHpath;
	int dHpathBL;

	/* Dimension of a full path hessian block */
	int dhP;

	/**
	** Parameters
	**/
	/* Regularization coef for building positive-definite Jacobi (diagonal) preconditioner */
	double jacRegCoef;

	/* Bandwidth for band preconditioner */
	int bw;
	int bbw; // bw+1 

	/* Pivot tolerance for perturbed LDL' of band preconditioner */
	double pivTol;

	/**
	** Local data, ALLOCATED
	**/
	/* Diagonal D (by blocks) */
	double* precD;

	/* Band preconditioner (band format) */
	double* pMbShoot;
	double* pMbPath;

	/**
	** Pointers to external data
	**/
	/* Pointer to shooting part of reduced AL hessian approx */
	double* hALshootRed;

	/* Pointer to slack part of reduced A hessian approx */
	double* hALpathRed;

	/* Pointer to set of free/active variables */
    int* freeVars;

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
};

#endif PRECOND_H