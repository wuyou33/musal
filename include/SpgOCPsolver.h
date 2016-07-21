#ifndef SPG_SOLVER_H
#define SPG_SOLVER_H


#include <iostream>
#include <math.h>

#include "MuShoot.h"
#include "SolverOptions.h"


/**
* Implements real-time primal augmented Lagrangian iterations based on spectral projected gradient
*/
template <class dynT,class costT,class pconT,class mayerT> class SpgOCPsolver {

	typedef MuShoot<dynT,costT,pconT,mayerT> shooterT;
	

public:
	/* Default constructor */
	SpgOCPsolver();

	/* Destructor */
	~SpgOCPsolver();

	/* Copy constructor */
	SpgOCPsolver(const SpgOCPsolver<dynT,costT,pconT,mayerT>&);

	/* Assignment operator */
	SpgOCPsolver<dynT,costT,pconT,mayerT>& operator= (const SpgOCPsolver<dynT,costT,pconT,mayerT>&);

	/**
	* Set methods 
	*/
	/* Initialize primal optimizer */
    void setPrimal(double*);
    
    /* Initialize dual optimizer */
    void setDual(double*);
        
	/* Set initial state (parameter) */
	inline void setXest(double* xest){
		shooter.setXest(xest);
	}
    
    /* Set output & input reference */
    inline void setReference(double* yRef,double* uRef){
        shooter.setReference(yRef,uRef);
    }
    
    /* Set shooter */
    inline void setShooter(const shooterT& shooter_){
        shooter = shooter_;
    }

    /* Set solver options */
    inline void setSolverOptions(SolverOptions& opts){
        
        /* Max. amount of primal iterations */
        maxPit = opts.getMaxPrimalIters();

        /* Initial halved penalty parameter */
        halfRhoIni = 0.5*opts.getInitialPenalty();

        /* Squared tolerances on inner KKT & 2-norm of equality constraints */
        kktTol2 = opts.getKktTolAbs()*opts.getKktTolAbs();
        necTol2 = opts.getNecTolAbs()*opts.getNecTolAbs();

        /* Tolerance on infinity-norm of momentum of primal sequence */
        diffTol = opts.getDiffTol();

        /* Constant for sufficient decrease condition in non-monotone line-search */
        gamSpg = opts.getSpgGam();
        
        /* History length for non-monotone line-search */
        histLenSpg = opts.getSpgHistLen();

        /* Safeguarding constants for non-monotone linea-search */
        sigOneSpg = opts.getSpgSigOne();
        sigTwoSpg = opts.getSpgSigTwo();

        /* Bounds on step-size */
        aminSpg = opts.getSpgAmin();
        amaxSpg = opts.getSpgAmax();

        /* Initial step-size for non-monotone line-search */
        ainiSpg = opts.getSpgAini();

        /* Maximum # backtracking iterations in non-monotone line-search */
        maxBckIt = opts.getSpgMaxBckIt();
    }

    /* Set initial halved penalty parameter */
    inline void setHalvedRhoIni(double halfRhoIni){
        this->halfRhoIni = halfRhoIni;
    }

    /* Set maximum # primal iterations */
    inline void setMaxPrimalIter(int maxPit){
        this->maxPit = maxPit;
    }

    /**
    * Get methods 
    */
    /* Outputs laved initial penalty parameter */
    inline double getHalvedRhoIni(){
        return halfRhoIni;
    }

    /* Outputs halved penalty parameter */
    inline double getHalvedRho(){
        return halfRho;
    }

    /* Returns primal optimizer */
    inline double* getPrimal(){
    	return z;
    }

    /* Returns next iterate */
    inline double* getPrimalNext(){
        return znex;
    }
    
    /* Returns dual optimizer */
    inline double* getDual(){
    	return mu;
    }
    
    /* Returns first optimal input */
    inline double* getNmpcLaw(){
    	return uopt;
    }

    /* Outputs # shooting nodes in shooter */
    inline int getNs(){
        return shooter.getNs();
    }

    /* Outputs # shooting intervals in shooter */
    inline int getNs_(){
        return shooter.getNs_();
    }

    /* Returns KKT satisfaction on inner (AL) NLP */
    inline double getKKT(){
        return kkt;
    }

    /* Returns feasibility */
    inline double getFeas(){
        return nec;
    }

    /**
    * Numerical processing methods 
    */
    /* Initialization method (set primal bounds, performs some precomputations,...)*/
    void init();
    
    /* Solver real-time routine (fixed # primal iterations & 1 dual update) */
    void solveRealTime();

    /* Solver real-time homotopy method */
    void solveHomotopy();
    
    /* Solver real-time primal loop (fixed # primal iterations) */
    void primalRealTime();
    
    /* Compute squared KKT on augmented Lagrangian (primal) problem */
    void evalKKTprim();
    
    /* Compute squared 2-norm of equality constraints */
    void evalEconNorm();

    /* Shift primal & dual variables (in an NMPC fashion) */
    void shiftPrimalDual();

    /* Cold-start primal-dual optimizer */
    void coldStartPrimalDual();

    /** 
    * Other methods 
    */
    /* Display members */
    void displayMembers();


private:
	/* Multiple-shooting object for evaluating augmented Lagrangian and its gradient */
	shooterT shooter;

	/**
	* Algorithm options 
	*/
	/* Max. # primal iterations */
	int maxPit;

    /* Coefficient for sufficient decrease in backtracking SPG search */
    double gamSpg;

    /* History length for non-monotone line-search */
    int histLenSpg;

    /* Safeguarding parameters */
    double sigOneSpg;
    double sigTwoSpg;

    /* Upper and lower bounds on step-size */
    double aminSpg;
    double amaxSpg;

    /* Initial step-size between aminSpg & amaxSpg */
    double ainiSpg;

    /* Maximum # backtracking iterations */
    int maxBckIt;


	/**
	* Algorithm control parameters
	*/
	/* Initial halved penalty parameter */
	double halfRhoIni;

	/* Halved penalty parameter */
	double halfRho;

	/* Squared KKT conditions on augmented Lagrangian problem */
	double kkt2;

    /* KKT conditions on augmented Lagrangian problem */
    double kkt;

	/* Squared 2-norm of equality constraints */
	double nec2;

    /* 2-norm of equality constraints */
    double nec;

	/* Squared KKT tolerance on primal loop */
    double kktTol2;

    /* Squared tolerance on 2-norm of equality constraints */
    double necTol2;
   
    /* Tolerance on infinity-norm of momentum of primal sequence */
    double diffTol;


    /**
    * Solver dimensions
    */
    /* # shooting intervals */
    int Ns_;

    /* Shooting primal optimizer dimension */
    int dZshoot;
    
    /* Slack primal optimizer dimension */
    int dZpath;
    
    /* Primal optimizer dimension */
    int dZ;
    int dZBL;
    
    /* Shooting dual optimizer dimension */
    int dLshoot;
    
    /* Slack dual optimizer dimension */
    int dLpath;
    
    /* Dual optimizer dimension */
    int dL;
    int dLBL;
    
    /* Input dimension */
    int du;
    int dubl;
    
    /* State dimension */
    int dx;

	/* State + input dimension */
	int dxu;
  
  	/* State + input + state dimesion */
    int dxux;
    
    /* Path-constraint dimension */
    int dpc;
    
    /* Dimension of separated shooting AL gradient */
    int dGshoot;
    int dGshootBL;
    
    /* Dimension of slacks gradient */
    int dGpath;
    int dGpathBL;
    
    /* Dimension of shooting part of AL hessian and path-constraints part of AL hessian */
    int dhS; // Dimension of shooting block
    int dhP; // Dimension of path-constraint block
    int dHshoot; // Shooting AL hessian dimension
    int dHpath; // Path-constraint AL hessian dimension


   	/**
    * Solver ALLOCATED variables 
    */
    /* Primal optimizer z = (s_0',q_0',s_1',q_1',...,s_{N-1}',q_{N-1}',s_N',r_0',...,r_{N-1}')'
    * 	- s_n: shooting state
    *	- q_n: shooting input
    *	- r_n: path slack
    */
    double* z;

    /* Next primal iterate */
    double* znex;
    
    /* Difference between primal iterates */
    double* dffz;

    /* Dual optimizer mu = (mu_0',...,mu_N',lmb_0',...,lmb_{N-1})' 
    *	- mu_0: dual on initial-value embedding (s_0-xest = 0)
    *	- mu_n: dual on shooting constraint (s_{n+1}-x(t_{n+1};s_n,q_n) = 0)
    *	- lmb_n: dual on path constraint (g(s_n,q_n)+r_n = 0)
    */
    double* mu;
    
    /* Equality constraints (init.-value embedding + shooting + slacked path-constraints) */
    double* econ;
    
    /* Optimal control law (first shooting input) */
    double* uopt;
    
    /* Bounds on primal variable */ 
    double* zmin;
    double* zmax;
    
    /* Full gradient of augmented Lagrangian */
    double* gAL;

    /* Difference between full augmented Lagrangian gradient */
    double* dffGal;

    /* Shooting part of AL gradient */
    double* gALshoot;

    /* Slacks part of AL gradient */
    double* gALpath;

    /* Difference of slacks part of AL gradient (used as auxiliary vector in shooter) */
    double* dffGpath;

    /* History of augmented Lagrangian values */
    double* vALhist;
};
#include "SpgOCPsolver.tpp"

#endif
