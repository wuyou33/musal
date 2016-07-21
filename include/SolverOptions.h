//
//  SolverOptions.h
//  MusAL
//
//  Created by Jean on 4/23/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_SolverOptions_h
#define MusAL_SolverOptions_h

#include <assert.h>

class SolverOptions {

public:
    SolverOptions();

    ~SolverOptions();
    
    /**
    * Set methods
    **/
    /* Set coefficient for interpolation in Cauchy phase */
    inline void setBetaInter(double betaInter){ 
        assert(betaInter>0);
        assert(betaInter<1);
        this->betaInter = betaInter;
    }

    /* Set coefficient for extrapolation in Cauchy phase */
    inline void setBetaExtra(double betaExtra){
        assert(betaExtra>1);
        this->betaExtra = betaExtra;
    }

    /* Set max # interpolation steps */
    inline void setMaxInterIter(int maxInterIter){
        assert(maxInterIter>=1);
        this->maxInterIter = maxInterIter;
    }

    /* Set max # extrapolation steps */
    inline void setMaxExtraIter(int maxExtraIter){
        assert(maxExtraIter>=1);
        this->maxExtraIter = maxExtraIter;
    }

    /* Set max # refinement steps (PCG+projected search)*/
    inline void setMaxRefIter(int maxRefIter){
        assert(maxRefIter>=1);
        this->maxRefIter = maxRefIter;
    }

    /* Set max # backtracking iter for projected search */
    inline void setMaxItSrch(int maxItSrch){
        assert(maxItSrch>=1);
        this->maxItSrch = maxItSrch;
    }

    /* Set bck coefficient for projected search */
    inline void setBetaSrch(double betaSrch){
        assert(betaSrch>0);
        assert(betaSrch<1);
        this->betaSrch = betaSrch;
    }

    /* Set constant for sufficient decrease in projected search */
    inline void setCstSrch(double cstSrch){
        assert(cstSrch>0);
        assert(cstSrch<1);
        this->cstSrch = cstSrch;
    }

    /* Set initial trust region radius */
    inline void setInitialTrustRadius(double trRadIni){
        assert(trRadIni>0);
        this->trRadIni = trRadIni;
    }

    /* Set tolerance on trust region radius */
    inline void setToleranceTrustRadius(double trRadTol){
        assert(trRadTol>0);
        this->trRadTol = trRadTol;
    }
    
    /* Set decision coef for trust region management */
    inline void setTrustRatioParameters(double eta0,double eta1,double eta2){
        assert(eta0>0);
        assert(eta0<1);
        assert(eta1>0);
        assert(eta1<1);
        assert(eta2>0);
        assert(eta2<1);
        assert(eta1<eta2);   
        this->eta0 = eta0;
        this->eta1 = eta1;
        this->eta2 = eta2;
    }
    
    /* Set update coefs for trust region management */
    inline void setTrustUpdateCoefs(double sig1,double sig2,double sig3){
        assert(sig1<sig2);
        assert(sig1>0);
        assert(sig2>0);
        assert(sig1<1);
        assert(sig2<1);
        assert(sig3>1);
        this->sig1 = sig1;
        this->sig2 = sig2;
        this->sig3 = sig3;
    }
    
    /* Set max # primal iterations */
    inline void setMaxPrimalIters(int maxPit){
        assert(maxPit>=1);
        this->maxPit = maxPit;
    }
    
    /* Set max # dual iterations */
    inline void setMaxDualIters(int maxDit){  
        assert(maxDit>=1);
        this->maxDit = maxDit;
    }
    
    /* Set multiplicative coef on penalty */
    inline void setPenaltyMulCoef(double mulPen){
        assert(mulPen>1);   
        this->mulPen = mulPen;
    }
    
    /* Set initial penalty */
    inline void setInitialPenalty(double rhoIni){
        assert(rhoIni>0);
        this->rhoIni = rhoIni;
    }
    
    /* Set initial step-size for activity detection */
    inline void setInitStepSize(double sini){  
        assert(sini>0);
        this->sini = sini;
    }
    
    /* Set backtracking coefficient for activity detection */
    inline void setBacktrackCoef(double beta){
        assert(beta<1);
        assert(beta>0);
        betaBck = beta;
    }
    
    /* Set backtracking constant (sufficient decrease) for activity detection */
    inline void setBacktrackConst(double c){
        assert(c<0.5);
        assert(c>0);
        cstBck = c;
    }
    
    /* Set maximum # backtracking iterations */
    inline void setBacktrackMaxIter(int maxIt){
        assert(maxIt>=1);
        maxItBck = maxIt;
    }
    
    /* Set tolerance on activity */
    inline void setActivityTolerance(double acTol){
        assert(acTol>0);
        assert(acTol<=1);
        activTol = acTol;
    }

    /* Set absolute KKT tolerance */
    inline void setKktTolAbs(double kktTolAbs){
        assert(kktTolAbs>0);
        this->kktTolAbs = kktTolAbs;
    }

    /* Set absolute tolerance on 2-norm of equality constraints */
    inline void setNecTolAbs(double necTolAbs){
        assert(necTolAbs>0);
        this->necTolAbs = necTolAbs;
    }

    /* Set tolerance on infinity-norm of momentum of primal sequence */
    inline void setDiffTol(double diffTol){
        assert(diffTol>0);
        this->diffTol = diffTol;
    }

    /* Set tolerance on magnitude of objective difference */
    inline void setObjTol(double objTol){
        assert(objTol>0);
        this->objTol = objTol;
    }

    /* Set SR1 scaling coefs */
    inline void setSrShootIni(double sr_shoot_ini){
        this->sr_shoot_ini = sr_shoot_ini;
    }

    inline void setSrPathIni(double sr_path_ini){
        this->sr_path_ini = sr_path_ini;
    }

    inline void setSrSkip(double sr_skip){
        this->sr_skip = sr_skip;
    }

    /* Set regularization coef for preconditionner */
    inline void setRegCoefPrec(double pr_reg_coef){
        assert(pr_reg_coef>0);
        this->pr_reg_coef = pr_reg_coef;
    }

    /* Set bandwidth of band preconditioner */
    inline void setBandWidth(int bwd){
        this->bwd = bwd;
    }

    /* Set pivot tolerance for LDL' factorization on band preconditioner */
    inline void setPivotTolerance(double pivTol){ 
        assert(pivTol>0);
        this->pivTol = pivTol;
    }

    /* Set primal optimality tolerance */
    inline void setKktTolIni(double kktTolIni){
        assert(kktTolIni>0);
        this->kktTolIni = kktTolIni;
    }

    /* Set equality constraints tolerance */
    inline void setEcTolIni(double ecTolIni){ 
        assert(ecTolIni>0);
        this->ecTolIni = ecTolIni;
    }

    /**
    * Get methods
    */    
    inline double getInitialTrustRadius(){
        return trRadIni;
    }

    inline double getToleranceTrustRadius(){
        return trRadTol;
    }
    
    inline double getEtaZero(){
        return eta0;
    }
    
    inline double getEtaOne(){
        return eta1;
    }
    
    inline double getEtaTwo(){
        return eta2;
    }
    
    inline double getSigmaOne(){
        return sig1;
    }
    
    inline double getSigmaTwo(){
        return sig2;
    }
    
    inline double getSigmaThree(){
        return sig3;
    }
    
    inline int getMaxPrimalIters(){
        return maxPit;
    }
    
    inline int getMaxDualIters(){
        return maxDit;
    }
    
    inline double getPenaltyMulCoef(){
        return mulPen;
    }
    
    inline double getInitialPenalty(){
        return rhoIni;
    }

    inline double getInitStepSize(){
        return sini;
    }
    
    inline double getBacktrackCoef(){
        return betaBck;
    }

    inline double getBacktrackConst(){
        return cstBck;
    }
    
    inline int getBacktrackMaxIter(){
        return maxItBck;
    }

    inline double getActivityTolerance(){
        return activTol;
    }

    inline double getSrShootIni(){
        return sr_shoot_ini;
    }

    inline double getSrPathIni(){
        return sr_path_ini;
    }

    inline double getSrSkip(){
        return sr_skip;
    }

    inline double getRegCoefPrec(){
        return pr_reg_coef;
    }

    inline int getBandWidth(){
        return bwd;
    }

    inline double getPivotTolerance(){
        return pivTol;
    }

    inline double getKktTolAbs(){
        return kktTolAbs;
    }

    inline double getNecTolAbs(){
        return necTolAbs;
    }

    inline double getDiffTol(){
        return diffTol;
    }

    inline double getObjTol(){
        return objTol;
    }

    inline double getSpgGam(){
        return gamSpg;
    }

    inline int getSpgHistLen(){
        return histLenSpg;
    }

    inline double getSpgSigOne(){
        return sigOneSpg;
    }

    inline double getSpgSigTwo(){
        return sigTwoSpg;
    }

    inline double getSpgAmin(){
        return aminSpg;
    }

    inline double getSpgAmax(){
        return amaxSpg;
    }

    inline double getSpgAini(){
        return ainiSpg;
    }

    inline int getSpgMaxBckIt(){
        return maxBckItSpg;
    }

    inline double getEcTolIni(){
        return ecTolIni;
    }

    inline double getKktTolIni(){
        return kktTolIni;
    }

    inline int getMaxRefIter(){
        return maxRefIter;
    }   

    inline int getMaxItSrch(){
        return maxItSrch;
    }

    inline double getBetaSrch(){
        return betaSrch;
    }

    inline double getCstSrch(){
        return cstSrch;
    }

    inline double getBetaInter(){
        return betaInter;
    }

    inline double getBetaExtra(){
        return betaExtra;
    }

    inline int getMaxInterIter(){
        return maxInterIter;
    }

    inline int getMaxExtraIter(){
        return maxExtraIter;
    }

private:
    /** 
    * Options for refinement search
    */
    /* Maximum # refinement iterations */
    int maxRefIter;

    /* Maximum # bck iter in projected search */
    int maxItSrch;

    /* Backtracking coefficient for projected search */
    double betaSrch;

    /* Sufficient decrease constant in projected search */
    double cstSrch;

    /**
    * Options for Cauchy search
    */
    /* Initial step-size */
    double sini;
    
    /* Backtracking multiplicative coefficient */
    double betaBck;
    
    /* Backtracking constant */
    double cstBck;
    
    /* Maximum # backtracking iterations */
    int maxItBck;

    /**
    * Options for interpolation-extrapolation
    */
    /* Interpolation coefficient */
    double betaInter;

    /* Extrapolation coefficient */
    double betaExtra;

    /* Max # interpolations */
    int maxInterIter;

    /* Max # extrapolations */
    int maxExtraIter; 

    /**
    * Options for trust-region iterations
    */ 
    /* Initial trust-region radius */
    double trRadIni;
    
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

    /**
    * Options for SR1 update
    */
    /* Initial diagonal scaling coef for shooting hessian approximations */
    double sr_shoot_ini;

    /* Initial diagonal scaling coef for path hessian approximations */
    double sr_path_ini;

    /* Constant for skipping condition |r'*s|>=cst*||r||*||s|| */
    double sr_skip;

    /** 
    * Options for preconditioner
    */
    /* Regularization coef to make Jacobi preconditioner positive definite */
    double pr_reg_coef;

    /* Bandwidth for band preconditioner */
    int bwd;

    /* Pivot tolerance for LDL' on band preconditioner */
    double pivTol;

    /**
    * Options on (dual) outer loop
    */
    /* Maximum # dual iterations */
    int maxDit;
    
    /* Multiplicative coefficient to be applied on penalty prm. */
    double mulPen;
    
    /* Initial penalty parameter */
    double rhoIni;

    /* Initial tolerance on primal optimality */
    double kktTolIni;

    /* Initial tolerance on equality constraints */
    double ecTolIni;

    /**
    * Default tolerances
    */
    /* Tolerance on activity */
    double activTol;

    /* Tolerance on trust region radius, if below abort */
    double trRadTol;

    /* Absolute tolerance on inner KKT */
    double kktTolAbs;

    /* Absolute tolerance on equality constraints */
    double necTolAbs;

    /* Tolerance on infinity-norm of momentum of primal sequence */
    double diffTol;

    /* Tolerance on magnitude of objective difference */
    double objTol;

    /**
    * Options for primal SPG loop
    */  
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

    /* Max. number of SPG backtracking iterations */    
    int maxBckItSpg;
};
    

#endif
