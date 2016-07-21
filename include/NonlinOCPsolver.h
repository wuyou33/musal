//
//  NonlinOCPsolver.h
//  MusAL
//
//  Created by Jean on 3/18/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//
#ifndef MusAL_NonlinOCPsolver_h
#define MusAL_NonlinOCPsolver_h

#include <iostream>
#include <math.h>

#include "defaultOptions.h"
#include "MuShoot.h"
#include "ActivityDetector.h"
#include "QuadraticModel.h"
#include "SolverOptions.h"
#include "PreconRefine.h"
#include "ProjectedRefiner.h"

/*** *** *** *** *** *** *** ***
* Nonlinear OCP solver class
* templated over multiple-shooting discretization type
** ** ** ** ** ** **/
template <class dynT,class costT,class pconT,class mayerT> class NonlinOCPsolver {

    typedef MuShoot<dynT,costT,pconT,mayerT> shooterT;

public:
    /* Default constructor */
    NonlinOCPsolver();

    /* Destructor */
    ~NonlinOCPsolver();

    /* Overloaded copy constructor */
    NonlinOCPsolver(const NonlinOCPsolver<dynT,costT,pconT,mayerT>&);

    /* Overloaded assignment operator */
    NonlinOCPsolver<dynT,costT,pconT,mayerT>& operator= (const NonlinOCPsolver<dynT,costT,pconT,mayerT>&);

    /**
    ** Set methods
    **/
    /* Initialize primal optimizer */
    inline void initPrimal(double*);
    
    /* Initialize dual optimizer */
    inline void initDual(double*);
        
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
    
    /* Set activity detector */
    inline void setActivityDetector(const ActivityDetector& activDetector){
        this->activDetector = activDetector;
    }

    /* Set preconditoned subspace refiner */
    inline void setPrecSubspaceRefiner(const PreconRefine& precSubRefiner){
        this->precSubRefiner = precSubRefiner;
    }

    /* Set preconditioned projected subspace refiner */
    inline void setPrecProjRefiner(const ProjectedRefiner& precProjRefiner){
        this->precProjRefiner = precProjRefiner;
    }

    /* Set quadratic model */
    inline void setQuadraticModel(const QuadraticModel& qmod){
        this->qmod = qmod;
    }
    
    /* Set MAT files */
    inline void setMATfiles(char* zFile,char* muFile){
        this->zFile = zFile;
        this->muFile = muFile;
    }
    
    /* Set solver options */
    inline void setSolverOptions(SolverOptions& opts){

        /* Options for Cauchy search */
        sini = opts.getInitStepSize();

        /* Options for trust-region loop */
        eta0 = opts.getEtaZero();
        eta1 = opts.getEtaOne();
        eta2 = opts.getEtaTwo();
        
        sig1 = opts.getSigmaOne();
        sig2 = opts.getSigmaTwo();
        sig3 = opts.getSigmaThree();
        
        maxPit = opts.getMaxPrimalIters();

        trRadTol = opts.getToleranceTrustRadius();
        
        /* Options for activity detection */
        activDetector.setBacktrackCoef(opts.getBacktrackCoef());
        activDetector.setBacktrackConst(opts.getBacktrackConst());
        activDetector.setBacktrackMaxIter(opts.getBacktrackMaxIter());
        activDetector.setActivityTolerance(opts.getActivityTolerance());
        activDetector.setDiffTolerance(opts.getDiffTol());
        activDetector.setBetaInter(opts.getBetaInter());
        activDetector.setBetaExtra(opts.getBetaExtra());
        activDetector.setMaxInterIter(opts.getMaxInterIter());
        activDetector.setMaxExtraIter(opts.getMaxExtraIter());

#if REFIN_PCG_RES
        /* Options for preconditioned subspace refiner */
        precSubRefiner.setBandWidth(opts.getBandWidth());
        precSubRefiner.setPivotTolerance(opts.getPivotTolerance());
#endif
#if REFIN_MOR_TOR
        /* Options for preconditioned projected refiner */
        precProjRefiner.setBandWidth(opts.getBandWidth());
        precProjRefiner.setPivotTolerance(opts.getPivotTolerance());
        precProjRefiner.setMaxRefIter(opts.getMaxRefIter());
        precProjRefiner.setMaxItSrch(opts.getMaxItSrch());
        precProjRefiner.setBetaSrch(opts.getBetaSrch());
        precProjRefiner.setCstSrch(opts.getCstSrch());
#endif
        /* Options for SR1 updates */
        sr_shoot_ini = opts.getSrShootIni();
        sr_path_ini = opts.getSrPathIni();

        sr_skip = opts.getSrSkip();
        sr_skip2 = sr_skip*sr_skip;

        /* Options for positive-definite preconditionner design */
        pr_reg_coef = opts.getRegCoefPrec();
        
        /* Pivot tolerance for LDL' on preconditioner */
        pivTol = opts.getPivotTolerance();

        /* Options for dual (outer) loop */
        maxDit = opts.getMaxDualIters();
        
        mulPen = opts.getPenaltyMulCoef();
        halfRhoIni = 0.5*opts.getInitialPenalty();

        kktTolIni = opts.getKktTolIni();
        ecTolIni = opts.getEcTolIni();

        /* Absolute tolerances on inner KKT & 2-norm of equality constraints */
        kktTolAbs = opts.getKktTolAbs();
        necTolAbs = opts.getNecTolAbs();

        /* Tolerance on infinity-norm of momentum of primal sequence */
        diffTol = opts.getDiffTol();

        /* Tolerance on magnitude of objective difference */
        objTol = opts.getObjTol();
    }
    
    /* Set initial halved penalty parameter */
    inline void setHalvedRhoIni(double halfRhoIni){
        this->halfRhoIni = halfRhoIni;
    }

    /* Set halved penalty parameter */
    inline void setHalvedRho(double halfRho){
        this->halfRho = halfRho;
    }

    /* Set maximum # primal iterations */
    inline void setMaxPrimalIter(int maxPit){
        this->maxPit = maxPit;
    }

    /* Set maximum # cumulative PCG iters in one primal run */
    inline void setMaxCumCG(const int maxCumCG){
        this->maxCumCG = maxCumCG;
        precSubRefiner.setMaxCumCG(maxCumCG);
    }

    /* Set trust-region radius */
    inline void setTrustRadius(double trRad){
        this->trRad = trRad;
    }

    inline void setBandWidth(int bw){
        precSubRefiner.setBandWidth(bw);
    }

    /* Set bounds on terminal state both equal to xest (from MuShoot) */
    void setPeriodicTerminal();

    /**
    ** Get methods
    **/
    /* Outputs halved penalty parameter */
    inline double getHalvedRho(){
        return halfRho;
    }

    /* Outputs penalty paramter */
    inline double getRho(){
        return 2.*halfRho;
    }

    /* Returns primal optimizer */
    inline double* getPrimal(){
    	return z;
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

    /* Returns global total # CG iterations */
    inline int getGlobCG(){
        return globCG;
    }

    /* Returns cumulative # CG iters over one primal loop */
    inline int getCumCG(){
        return cumCG;
    }

    /* Returns global total # CG restarts */
    inline int getGlobRestart(){
        return globRestart;
    }

    /* Returns cumulative # restarts over one primal loop */
    inline int getCumRestart(){
        return cumRestart;
    }

    /**
    ** Numerical processing methods 
    **/
    /* Initialization method (set primal bounds, performs some precomputations,...)*/
    void init();
    
    /* Solver real-time routine (fixed # primal iterations & 1 dual update) */
    void solveRealTime(const double);
    
    /* Solver main routine (classic augmented Lagrangian loop) */
    void solveFull();

    /* LANCELOT loop */
    void solveLancelot();
    
    /* Solver real-time primal loop (fixed # primal iterations) */  
    void primalRealTime(const double);
    
    /* Solver full primal loop */
    void primalFull(const double);
    
    /* Compute squared KKT on augmented Lagrangian (primal) problem */
    void evalKKTprim();
    
    /* Compute squared 2-norm of equality constraints */
    void evalEconNorm();

    /* Shift primal & dual variables (in an NMPC fashion) */
    void shiftPrimalDual();

    /* Cold-start primal-dual optimizer */
    void coldStartPrimalDual();

    /* Set initial hessian approximation */
    void initHessianApprox();

    /* Compute AL gradient at current primal-dual iterate & penalty parameter */
    void evalALgradient();

    /* Assemble reduced AL hessian */
    void assembleReducedHessian();

    /* Initialize shooting variables for bioreactor */
    void initBioReactor();

    /** 
    ** For debugging 
    **/
    /* Read data from MAT files and store them in z & mu */
    void readMATfiles();

    /* Display members */
    void displayMembers(const int);

    /* Display iterate z */
    void displayIterate();

    /* Display candidate zS */
    void displayCandidate();

    /* Display dual */
    void displayDual();

    /* Display shooting part of AL hessian */
    void displayHalShoot();

    /* Display NLP bounds */
    void displayBounds();

    /* Display initial state xest */
    void displayXest();
#if SCAL_VAR
    /* Display state transfo */
    void displayStateTransfo();

    /* Display input transfo */
    void displayInputTransfo();
#endif

private:
    /**
    ** Processing objects
    **/
    /* Multiple-shooting object (interface for evaluating augmented Lagrangian and its gradient) */
    shooterT shooter;
    
    /* Activity detector */
    ActivityDetector activDetector;

    /* Preconditioned subspace refiner */
    PreconRefine precSubRefiner;

    /* Preconditioned projected refiner (Moré-Toraldo projections, see TRON)) */
    ProjectedRefiner precProjRefiner;

    /* Quadratic model */
    QuadraticModel qmod;
    
    /********************************************************************/
    /************************* Algorithm options ************************/
    /********************************************************************/
    
    /** Options on (primal) inner loop
    */
    /* Initial step-size for Cauchy search */
    double sini;

    /* Test parameters for trust ratio */
    double eta0;
    double eta1;
    double eta2;
    
    /* Coefficients for trust-region radius update */
    double sig1;
    double sig2;
    double sig3;
    
    /* Maximum # primal iterations */
    int maxPit;

    /* Max. # cumulative PCG iter (in one primal run) */
    int maxCumCG;

    /* Options for SR1 update */
    double sr_shoot_ini;
    double sr_path_ini;
    double sr_skip;
    double sr_skip2;

    /* Reg. coef. for making preconditionner positive definite */
    double pr_reg_coef;

    /* Pivot tolerance for LDL' factorization on band preconditioner */
    double pivTol;

    /* Tolerance on trust-region radius */
    double trRadTol;
    
    /** Options on (dual) outer loop
    */
    /* Maximum # dual iterations */
    int maxDit;
    
    /* Multiplicative coefficient to be applied on penalty prm. */
    double mulPen;
    
    /**
    ** Algorithm control parameters
    **/
    /* Initial halved penalty parameter */
    double halfRhoIni;

    /* Halved penalty parameter */
    double halfRho;
    
    /* Trust-region radius */
    double trRad;
    
    /* Primal KKT satisfaction */
    double kkt;

    /* Squared primal KKT satisfaction */
    double kkt2;
    
    /* Absolute KKT tolerance for augmented Lagrangian subproblem */
    double kktTolAbs;

    /* Initial tolerance on KKT conditions of inner problem (LANCELOT) */
    double kktTolIni;

    /* 2-norm of equality constraints */
    double nec;

    /* Squared 2-norm of equality constraints */
    double nec2;

    /* Absolute tolerance on 2-norm of equality constraints */
    double necTolAbs;

    /* Initial tolerance on 2-norm of equality constrains (LANCELOT) */
    double ecTolIni;

    /* Squared CG tolerance (updated at every primal iteration) */
    double cgTol2;

    /* Tolerance on infinity-norm of momentum of primal sequence */
    double diffTol;

    /* Tolerance on magnitude of objective difference */
    double objTol;

    /* Infinity norm of z-zS (stored in negstp)*/
    double norIstp;

    /* Euclidean norm of z-zS */
    double nor2stp;
    
    /* Total number of CG iterations per primal loop */
    int cumCG;

    /* Global total number of CG iterations */
    int globCG;

    /* Total number of restarts per primal loop */
    int cumRestart;

    /* Global total number of restarts */
    int globRestart;

    /* Total number of primal iterations for one dual iter */
    int numPrimIt;

    /* Indicator of error in primal loop */
    int errPrim;

    /**
    ** Solver dimensions
    **/
    /* # shooting intervals */
    int Ns_;

    /* # shooting nodes */
    int Ns;

    /* Shooting primal optimizer dimension */
    int dZshoot;
    int dZshootBL;
    
    /* Slack primal optimizer dimension */
    int dZpath;
    int dZpathBL;

    /* Primal optimizer dimension */
    int dZ;
    int dZBL;
    int dZNT;
    
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
    int dxbl;

	/* State + input dimension */
	int dxu;
  
  	/* State + input + state dimesion */
    int dxux;
    
    /* Path-constraint dimension */
    int dpc;
    
    /* Dimension of separated shooting AL gradient */
    int dGshoot;
    int dGshootBL;
    int dGshootNT;
    
    /* Dimension of slacks gradient */
    int dGpath;
    int dGpathBL;
    int dGpathNT;
    
    /* Dimension of shooting part of AL hessian and path-constraints part of AL hessian */
    int dhS; // Dimension of shooting block
    int dhP; // Dimension of path-constraint block
    int dHshoot; // Shooting AL hessian dimension
    int dHpath; // Path-constraint AL hessian dimension
        
    /*******************************************************************************************/
    /******************************** Solver ALLOCATED variables *******************************/
    /*******************************************************************************************/
    /* Primal optimizer z = (s_0',q_0',s_1',q_1',...,s_{N-1}',q_{N-1}',s_N',r_0',...,r_{N-1}')'
    * 	- s_n: shooting state
    *	- q_n: shooting input
    *	- r_n: path slack
    */
    double* z;
    
    /* Primal optimizer at previous iteration */
    double* zprev;
    
    /* Negative difference z-zS between current iterate and previous iterate (for SR1 update) */
    double* negstp;
    
    /* Candidate point */
    double* zS;
    
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
    
#ifdef BIOREC
    /* Auxiliary vector of size dx+du (primal) */
    double* uxAux;

    /* Auxiliary vector of size dx (dual) */
    double* xAux;

    /* Auxiliary vector of size dpc (primal) */
    double* pcAux;

    /* Auxiliary vector of size dpc (dual) */
    double* pcAuxx;
#endif

    /* Bounds on primal variable (fixed) */ 
    double* zmin;
    double* zmax;
    
    /* Full gradient of augmented Lagrangian at iterate */
    double* gAL;
    
    /* Separated gradient of shooting part of AL at iterate */
    double* gALshoot;

    /* Separated gradient of shooting part of AL at candidate */
    double* gALshootS;
    
    /* Difference between separated shooting AL gradient at candidate and separated shooting
	   AL gradient at iterate */
    double* dffGshoot;
    
    /* Separated slack part of full AL gradient at iterate */
    double* gALpath;

    /* Separated slack part of full AL gradient at candidate */
    double* gALpathS;

    /* Difference between slack AL gradient at candidate and slack AL gradient at iterate */
    double* dffGpath;
    
    /* Separated shooting part of hessian of augmented Lagrangian */
    double* hALshoot;
    
    /* Separated shooting part of reduced AL hessian approx. */
    double* hALshootRed;

    /* Slack part of hessian of augmented Lagrangian (stored in a separable fashion) */
    double* hALpath;
    
    /* Separated slack part of reduced AL hessian approx. */
    double* hALpathRed;

    /* Auxiliary vector for sr1 update on shooting AL hessian */
    double* secVecShoot;

    /* Auxiliary vector for sr1 update on path-constraint AL hessian */
    double* secVecPath;

    /* Product B*(z-x) */
    double* prdHss;

    /* Indices of free variables */
    int* indFreeVars;

    /* Vector containing information about free/active variables over the full primal variable */
    int* freeVars;

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

    /**********************************************************************************
     ********************** File names for primal-dual initial guess ******************
     **********************************************************************************/

    /* MAT file name for primal initial guess */
    char* zFile;
    
    /* MAT file for dual initial guess */
    char* muFile;
    
};
#include "NonlinOCPsolver.tpp"


#endif
