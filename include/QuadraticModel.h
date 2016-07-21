//
//  QuadraticModel.h
//  MusAL
//
//  Created by Jean on 4/9/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#ifndef MusAL_QuadraticModel_h
#define MusAL_QuadraticModel_h

#include <iostream>

#include "macros.h"

/*******
* Implements quadratic model of objective reduction in trust region loop:
*   q(y;x) := gAL(x)'*(y-x)+0.5*(y-x)'*B(x)*(y-x)
************************/
class QuadraticModel {
    
public:
    /* Default constructor */
    QuadraticModel();
    
    /* Copy constructor */
    QuadraticModel(const QuadraticModel&);
    
    /* Assigment operator */
    QuadraticModel& operator= (const QuadraticModel&);
    
    /* Destructor */
    ~QuadraticModel();
    
    /**
     * Get methods
     */
    /* Returns model dimension */
    inline int getDimension(){
    	return dZ;
    }
    
    /* Returns pointer to full AL gradient */
    inline double* getGradient(){
        return grad;
    }
    
    /* Returns pointer to shooting hessian */
    inline double* getHessianShoot(){
        return hessShoot;
    }
    
    /* Returns pointer to path-constraint hessian */
    inline double* getHessianPath(){
        return hessPath;
    }
     
    /**
    * Set methods
    */
    /* Set dimension */
    inline void setDimensions(int dZ,int dZshoot,int dxdu,int dxdudx,int dp){
        this->dZ = dZ;
        dZBL = dZ*SZDBL;
        this->dZshoot = dZshoot;
        this->dxdu = dxdu;
        this->dxdudx = dxdudx;
        ddShoot = dxdudx*dxdudx;
        ddPath = dp*dp;
        dpc = dp;
    }
    
    /* Set # shooting intervals */
    inline void setNumberShoots(int Ns_){
        this->Ns_ = Ns_;
    }
    
    /* Set iterate */
    inline void setIterate(double* z){
        this->z = z;
    }
    
    /* Set Cauchy point */
    inline void setCauchy(double* zC){
        this->zC = zC;
    }
    
    /* Set candidate */
    inline void setCandidate(double* zS){
        this->zS = zS;
    }
    
    /* Set pointer to full AL gradient */
    inline void setGradient(double* grad){
        this->grad = grad;
    }
    
    /* Set hessian of shooting part */
    inline void setHessianShoot(double* hessShoot){
        this->hessShoot = hessShoot;
    }
    
    /* Set hessian of path-constraint part */
    inline void setHessianPath(double* hessPath){
        this->hessPath = hessPath;
    }

    /* Set pointer to product B*(z-x) */
    inline void setProductHessian(double* prdHss){
        this->prdHss = prdHss;
    }

 	/**
    * Processing methods 
    */
    /* Evaluate model */
    double evalModel(double*,double*);

    /* Evaluate model at Cauchy point */
    double evalModelCauchy();
    
    /* Evaluate model at candidate point */
    double evalModelCandidate();
    
    /* Evaluate quadratic term of model 0.5*(zC-z)'*B(z)*(zC-z)
        at Cauchy points */
    double evalQuadTermCauchy();

    /* Evaluate quadratic term of model 0.5*(y-z)'*B(z)*(y-z) */
    double evalQuadTerm(double*);

    /* (v-w)'*H*(v-w) */
    double multHshootDD(double* v,double* w);

    /**
    ** For debugging
    **/
    void displayHessian();

    void displayGradient();

private:
    /* Variable dimension */
    int dZ;
    int dZBL;
    
    /* Shooting variable dimension */
    int dZshoot;
    
    /* Other useful dimensions */
    int dxdu;
    int dxdudx;
    int dpc;
    int ddShoot;
    int ddPath;
    int dHshoot;
    
    /* Number of shooting intervals */
    int Ns_;
    
    /* Pointer to iterate */
    double* z;
    
    /* Pointer to Cauchy point */
    double* zC;
    
    /* Pointer to candidate */
    double* zS;
    
    /* Pointer to full AL gradient */
    double* grad;
    
    /* Pointer to shooting part of hessian */
    double* hessShoot;

    /* Pointer to path-constraint part of hessian */
    double* hessPath;

    /* Pointer to product B*(z-x) */
    double* prdHss;
};

#endif
