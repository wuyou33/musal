//
//  MuShoot.cpp
//  MusAL
//
//  Created by Jean on 3/12/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//
#include <iostream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "macros.h"
#include "userData.h"
#include "MuShoot.h"

/**
** Default constructor 
**/
template <class dynT,class costT,class pconT,class mayerT> MuShoot<dynT,costT,pconT,mayerT>::MuShoot(){

    Tpred = -1.;

    dx = ode.getDx();
    dxbl = dx*SZDBL;
    if (dx<=0){
        std::cerr<<"ERROR<MuShoot>: state-dimension nonpositive."<<std::endl;
        return;
    }

    du = ode.getDu();
    dubl = du*SZDBL;
    if (du<=0){
        std::cerr<<"ERROR<MuShoot>: input-dimension nonpositive."<<std::endl;
        return;
    }
    
    dy = cost.getDy();
    dybl = dy*SZDBL;
    if (dy<=0){
        std::cerr<<"ERROR<MuShoot>: output-dimension nonpositive."<<std::endl;
        return;
    }

    dpc = pcon.getDpc();
    dpcbl = dpc*SZDBL;
    dxu = dx+du;
    dxux = dxu+dx;
    dw = dx+1;
    da = dw+du;
    dper = DPER;

    augStates = NULL;
    augAdjoints = NULL;
    yRef = NULL;
    uRef = NULL;
    econPath = new double[dpc];
    adfin = NULL;
    umin = new double[du];
    umax = new double[du];
    xmin = new double[dx];
    xmax = new double[dx];
    xminT = new double[dx];
    xmaxT = new double[dx];
    augStateIni = new double[dw];
    *(augStateIni+dx) = 0.; // VERY VERY IMPORTANT
    usc = new double[du];
    integAugStates = NULL;
    interpAugStates = NULL;

    xest = NULL;
    z = NULL;
    zS = NULL;
    mu = NULL;
    econ = NULL;
    gAL = NULL;
    gALshoot = NULL;
    gALshootS = NULL;
    dffGpath = NULL;
    gALpath = NULL;
    gALpathS = NULL;

    /* State scaling transformation */
    Ds = new double[dx];
    iDs = new double[dx];
    cs = new double[dx];

    /* Input scaling transformation */
    Dq = new double[du];
    iDq = new double[du];
    cq = new double[du];
}

/**
** Copy constructor 
**/
template <class dynT,class costT,class pconT,class mayerT> MuShoot<dynT,costT,pconT,mayerT>::MuShoot(const MuShoot<dynT,costT,pconT,mayerT>& obj) : erkInteg(obj.erkInteg),ode(obj.ode),cost(obj.cost),pcon(obj.pcon),mayer(obj.mayer){
    
    int s,nsteps;

    Tpred = obj.Tpred;
    
    Ns = obj.Ns;
    Ns_ = obj.Ns_;

    dx = obj.dx;
    dxbl = obj.dxbl;
    
    du = obj.du;
    dubl = obj.dubl;

	dy = obj.dy;
	dybl = obj.dybl;

    dxu = obj.dxu;
    dxux = dxu+dx;
    dw = obj.dw;
    da = obj.da;
    dper = obj.dper;

    dpc = obj.dpc;
    dpcbl = obj.dpcbl;

    dZ = obj.dZ;
    dZshoot = obj.dZshoot;
    dZpath = obj.dZpath;

    izTerm = obj.izTerm;
    igTerm = obj.igTerm;
    igPath = obj.igPath;

    dL = obj.dL;
    dLpath = obj.dLpath;
    dLshoot = obj.dLshoot;

    dGshoot = obj.dGshoot;
    dGshootDBL = obj.dGshootDBL;
    dGpathDBL = obj.dGpathDBL;
    
    dAugS = obj.dAugS;
    dAugA = obj.dAugA;

    dwnsteps = obj.dwnsteps;
    sdwnsteps = obj.sdwnsteps;

    augStateIni = new double[dw];
    std::copy(obj.augStateIni,obj.augStateIni+obj.dw,augStateIni);

    usc = new double[du];
    std::copy(obj.usc,obj.usc+obj.du,usc);
    
    econPath = new double[dpc];
    std::copy(obj.econPath,obj.econPath+obj.dpc,econPath);
    
    umin = new double[du];
    std::copy(obj.umin,obj.umin+obj.du,umin);
    
    umax = new double[du];
    std::copy(obj.umax,obj.umax+obj.du,umax);
    
    xmin = new double[dx];
    std::copy(obj.xmin,obj.xmin+obj.dx,xmin);
    
    xmax = new double[dx];
    std::copy(obj.xmax,obj.xmax+obj.dx,xmax);

    xminT = new double[dx];
    std::copy(obj.xminT,obj.xminT+obj.dx,xminT);

    xmaxT = new double[dx];
    std::copy(obj.xmaxT,obj.xmaxT+obj.dx,xmaxT);

    adfin = new double[Ns_*da];
    std::copy(obj.adfin,obj.adfin+obj.Ns_*obj.da,adfin);

    yRef = new double[Ns*dy];
    std::copy(obj.yRef,obj.yRef+obj.Ns*obj.dy,yRef);

    uRef = new double[Ns*du];
    std::copy(obj.uRef,obj.uRef+obj.Ns*obj.du,uRef);

    augStates = new double[Ns_*dw];
    std::copy(obj.augStates,obj.augStates+obj.Ns_*obj.dw,augStates);

    augAdjoints = new double[da];
    std::copy(obj.augAdjoints,obj.augAdjoints+obj.da,augAdjoints);

    integAugStates = new double[Ns_*dwnsteps];
    std::copy(obj.integAugStates,obj.integAugStates+obj.Ns_*obj.dwnsteps,integAugStates);

    interpAugStates = new double[Ns_*sdwnsteps];
    std::copy(obj.interpAugStates,obj.interpAugStates+obj.Ns_*obj.sdwnsteps,interpAugStates);

    Ds = new double[dx];
    std::copy(obj.Ds,obj.Ds+obj.dx,Ds);

    iDs = new double[dx];
    std::copy(obj.iDs,obj.iDs+obj.dx,iDs);

    Dq = new double[du];
    std::copy(obj.Dq,obj.Dq+obj.du,Dq);

    iDq = new double[du];
    std::copy(obj.iDq,obj.iDq+obj.du,iDq);

    cs = new double[dx];
    std::copy(obj.cs,obj.cs+obj.dx,cs);

    cq = new double[du];
    std::copy(obj.cq,obj.cq+obj.du,cq);

    xest = obj.xest;
    z = obj.z;
    zS = obj.zS;
    mu = obj.mu;
    econ = obj.econ;
    gAL = obj.gAL;
    gALshoot = obj.gALshoot;
    gALshootS = obj.gALshootS;
    dffGpath = obj.dffGpath;
    gALpath = obj.gALpath;
    gALpathS = obj.gALpathS;
}

/**
** Assignment operator 
**/
template <class dynT,class costT,class pconT,class mayerT> MuShoot<dynT,costT,pconT,mayerT>& MuShoot<dynT,costT,pconT,mayerT>::operator=(const MuShoot<dynT,costT,pconT,mayerT>& obj){

    int s,nsteps;
    double *adfinCpy,*uminCpy,*umaxCpy,*xminCpy,*xmaxCpy;
    double *augStateIniCpy;
    double *yRef_,*uRef_;
    double *augStates_,*augAdjoints_;
    double *econPath_;
    double *integAugStates_,*interpAugStates_;
    double *Ds_,*iDs_,*cs_;
    double *Dq_,*iDq_,*cq_;
    double *usc_;
    double *xminT_,*xmaxT_;

    if (this != &obj){
        erkInteg = obj.erkInteg;
        ode = obj.ode;
        cost = obj.cost;
        pcon = obj.pcon;
        mayer = obj.mayer;
        
        Tpred = obj.Tpred;
        
        Ns = obj.Ns;
        Ns_ = obj.Ns_;
        
        /* Primal dimensions */
        dZ = obj.dZ;
        dZshoot = obj.dZshoot;
        dZpath = obj.dZpath;

        /* Dual dimensions */
        dL = obj.dL;
        dLshoot = obj.dLshoot;
        dLpath = obj.dLpath;
        
        /* Gradient dimensions */
        dGshoot = obj.dGshoot;
        dGshootDBL = obj.dGshootDBL;
        dGpathDBL = obj.dGpathDBL;

        /* State dimension */
        dx = obj.dx;
        dxbl = obj.dxbl;
        
        /* Input dimension */
        du = obj.du;
        dubl = obj.dubl;
        
        /* Output dimension */
        dy = obj.dy;
        dybl = obj.dybl;

        /* Periodicity constraints dimension */
        dper = obj.dper;
        
        /* Shooting node dimension */
        dxu = obj.dxu;
        
        /* Separable gradient dimension */
        dxux = obj.dxux;

        /* Augmented-state dimension */
        dw = obj.dw;
        
        /* Adjoint dimension */
        da = obj.da;
        
        /* Path constraint dimension */
        dpc = obj.dpc;
        dpcbl = obj.dpcbl;
        
        /* Augmented states & adjoints dimensions */
        dAugS = obj.dAugS;
        dAugA = obj.dAugA;

        /* Dummy indices */
        izTerm = obj.izTerm;
        igTerm = obj.igTerm;
        igPath = obj.igPath;

        /* Dimensions of local integration & interpolation data */
        dwnsteps = obj.dwnsteps;
        sdwnsteps = obj.sdwnsteps;

        /* Sequence of augmented adjoints */
        if (adfin != NULL)
            delete [] adfin;
        if (obj.adfin != NULL){
            adfinCpy = new double[Ns_*da];
            std::copy(obj.adfin,obj.adfin+Ns_*da,adfinCpy);
            adfin = adfinCpy;
        }
        else {
            adfin = NULL;
        }

        /* Bounds on inputs and states */
        if (umin != NULL)
            delete [] umin;
        if (obj.umin != NULL){
            uminCpy = new double[du];
            std::copy(obj.umin,obj.umin+du,uminCpy);
            umin = uminCpy;
        } 
        else {
            umin = NULL;
        }
        
        /* */
        if (umax != NULL)
            delete [] umax;
        if (obj.umax != NULL){
            umaxCpy = new double[du];
            std::copy(obj.umax,obj.umax+du,umaxCpy);        
            umax = umaxCpy;
        } 
        else {
            umax = NULL;
        }

        /* */
        if (xmin != NULL)
            delete [] xmin;
        if (obj.xmin != NULL){
            xminCpy = new double[dx];
            std::copy(obj.xmin,obj.xmin+dx,xminCpy);
            xmin = xminCpy;
        } 
        else {
            xmin = NULL;
        }

        /* */
        if (xmax != NULL)
            delete [] xmax;
        if (obj.xmax != NULL){
            xmaxCpy = new double[dx];
            std::copy(obj.xmax,obj.xmax+dx,xmaxCpy);
            xmax = xmaxCpy;
        }
        else {
            xmax = NULL;
        }
        
        /* */
        if (xminT != NULL)
            delete [] xminT;
        if (obj.xminT != NULL){
            xminT_ = new double[dx];
            std::copy(obj.xminT,obj.xminT+dx,xminT_);
            xminT = xminT_;
        }
        else {
            xminT = NULL;
        }

        /* */
        if (xmaxT != NULL)
            delete [] xmaxT;
        if (obj.xmaxT != NULL){
            xmaxT_ = new double[dx];
            std::copy(obj.xmaxT,obj.xmaxT+dx,xmaxT_);
            xmaxT = xmaxT_;
        }
        else {
            xmaxT = NULL;
        }

        /* Allocated space for initial augmented state (state+cost)*/
        if (augStateIni != NULL)
            delete [] augStateIni;
        if (obj.augStateIni != NULL){
            augStateIniCpy = new double[dw];
            std::copy(obj.augStateIni,obj.augStateIni+dw,augStateIniCpy);
            augStateIni = augStateIniCpy;
        }
        else {
            augStateIni = NULL;
        }

        /* Scaled input */
        if (usc != NULL)
            delete [] usc;
        if (obj.usc != NULL){
            usc_ = new double[du];
            std::copy(obj.usc,obj.usc+du,usc_);
            usc = usc_;
        }
        else {
            usc = NULL;
        }
        
        /* Output reference */
        if (yRef != NULL)
            delete [] yRef;
        if (obj.yRef != NULL){
            yRef_ = new double[Ns_*dy];
            std::copy(obj.yRef,obj.yRef+Ns_*dy,yRef_);
            yRef = yRef_;
        }
        else {
            yRef = NULL;
        }

        /* Input reference */
        if (uRef != NULL)
            delete [] uRef;
        if (obj.uRef != NULL){
            uRef_ = new double[Ns_*du];
            std::copy(obj.uRef,obj.uRef+Ns_*du,uRef_);
            uRef = uRef_;
        } 
        else {
            uRef = NULL;
        }

        /* Augmented states */
        if (augStates != NULL)
            delete [] augStates;
        if (obj.augStates != NULL){
            augStates_ = new double[Ns_*dw];
            std::copy(obj.augStates,obj.augStates+Ns_*dw,augStates_);
            augStates = augStates_;
        }
        else {
            augStates = NULL;
        }

        /* Augmented adjoints */
        if (augAdjoints != NULL)
            delete [] augAdjoints;
        if (obj.augAdjoints != NULL){
            augAdjoints_ = new double[da];
            std::copy(obj.augAdjoints,obj.augAdjoints+da,augAdjoints_);
            augAdjoints = augAdjoints_;
        }
        else {
            augAdjoints = NULL;
        }

        /* Dummy vector for path-constraints */
        if (econPath != NULL)
            delete [] econPath;
        if (obj.econPath != NULL){
            econPath_ = new double[dpc];
            std::copy(obj.econPath,obj.econPath+dpc,econPath_);
            econPath = econPath_;
        }
        else {
            econPath = NULL;
        }

        /* Array of intermediate integration augmented states */
        if (integAugStates != NULL)
            delete [] integAugStates;
        if (obj.integAugStates != NULL){
            integAugStates_ = new double[Ns_*dwnsteps];
            std::copy(obj.integAugStates,obj.integAugStates+obj.Ns_*obj.dwnsteps,integAugStates_);
            integAugStates = integAugStates_;
        }
        else {
            integAugStates = NULL;
        }    

        /* Array of interpolated augmented states */
        if (interpAugStates != NULL)
            delete [] interpAugStates;
        if (obj.interpAugStates != NULL){
            interpAugStates_ = new double[Ns_*sdwnsteps];
            std::copy(obj.interpAugStates,obj.interpAugStates+obj.Ns_*obj.sdwnsteps,interpAugStates_);
            interpAugStates = interpAugStates_;
        }
        else {
            interpAugStates = NULL;
        }

        /* State scaling transformation */
        if (Ds != NULL)
            delete [] Ds;
        if (obj.Ds != NULL){
            Ds_ = new double[dx];
            std::copy(obj.Ds,obj.Ds+obj.dx,Ds_);
            Ds = Ds_;
        }
        else {
            Ds = NULL;
        }

        if (iDs != NULL)
            delete [] iDs;
        if (obj.iDs != NULL){
            iDs_ = new double[dx];
            std::copy(obj.iDs,obj.iDs+obj.dx,iDs_);
            iDs = iDs_;
        }
        else {
            iDs = NULL;
        }

        if (cs != NULL)
            delete [] cs;
        if (obj.cs != NULL){
            cs_ = new double[dx];
            std::copy(obj.cs,obj.cs+obj.dx,cs_);
            cs = cs_;
        }
        else {
            cs = NULL;
        }

        /* Input scaling transformation */
        if (Dq != NULL)
            delete [] Dq;
        if (obj.Dq != NULL){
            Dq_ = new double[du];
            std::copy(obj.Dq,obj.Dq+obj.du,Dq_);
            Dq = Dq_;
        }
        else {
            Dq = NULL;
        }

        if (iDq != NULL)
            delete [] iDq;
        if (obj.iDq != NULL){
            iDq_ = new double[du];
            std::copy(obj.iDq,obj.iDq+obj.du,iDq_);
            iDq = iDq_;
        }
        else {
            iDq = NULL;
        }

        if (cq != NULL)
            delete [] cq;
        if (obj.cq != NULL){
            cq_ = new double[du];
            std::copy(obj.cq,obj.cq+obj.du,cq_);
            cq = cq_;
        }
        else {
            cq = NULL;
        }

        /* Pointers */
        xest = obj.xest;
        z = obj.z;
        zS = obj.zS;
        mu = obj.mu;
        econ = obj.econ;
        gAL = obj.gAL;
        gALshoot = obj.gALshoot;
        gALshootS = obj.gALshootS;
        dffGpath = obj.dffGpath;
        gALpath = obj.gALpath;
        gALpathS = obj.gALpathS;
    }
    
    return *this;
}

/**
** Destructor 
**/
template <class dynT,class costT,class pconT,class mayerT> MuShoot<dynT,costT,pconT,mayerT>::~MuShoot(){

    if (augStateIni != NULL)
        delete [] augStateIni;

    if (usc != NULL)
        delete [] usc;

    if (umin != NULL)
        delete [] umin;
        
    if (umax != NULL)
        delete [] umax;
        
    if (xmin != NULL)
        delete [] xmin;
        
    if (xmax != NULL)
        delete [] xmax;

    if (xminT != NULL)
        delete [] xminT;

    if (xmaxT != NULL)
        delete [] xmaxT;

    if (econPath != NULL)
        delete [] econPath;

    if (adfin != NULL)
        delete [] adfin;
    
    if (yRef != NULL)
        delete [] yRef;
        
    if (uRef != NULL)
        delete [] uRef;
        
    if (augStates != NULL)
        delete [] augStates;

    if (augAdjoints != NULL)
        delete [] augAdjoints;

    if (integAugStates != NULL)
        delete [] integAugStates;

    if (interpAugStates != NULL)
        delete [] interpAugStates;
 
    if (Ds != NULL)
        delete [] Ds;

    if (iDs != NULL)
        delete [] iDs;

    if (Dq != NULL)
        delete [] Dq;

    if (iDq != NULL)
        delete [] iDq;

    if (cq != NULL)
        delete [] cq;

    if (cs != NULL)
        delete [] cs;
}

/**
** Set bounds on input
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::setUmin(double* umin){
    assert(du>0);
    std::copy(umin,umin+du,this->umin);
}

template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::setUmax(double* umax){
    assert(du>0);
    std::copy(umax,umax+du,this->umax);
}

/**
** Set bounds on state
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::setXmin(double* xmin){
    assert(dx>0);
    std::copy(xmin,xmin+dx,this->xmin);
}

template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::setXmax(double* xmax){
    assert(dx>0);
    std::copy(xmax,xmax+dx,this->xmax);
}

/**
** Set bounds on terminal state 
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::setXminT(double* xminT){
    assert(dx>0);
    std::copy(xminT,xminT+dx,this->xminT);
}

template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::setXmaxT(double* xmaxT){
    assert(dx>0);   
    std::copy(xmaxT,xmaxT+dx,this->xmaxT);
}

/** 
** Set output & input references
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::setReference(double* yRef_,double* uRef_){
    assert(dy>0);
    assert(du>0);
    assert(Ns>0);

    /* Copy output reference */
    std::copy(yRef_,yRef_+dy*Ns,yRef);
    /* Copy input reference */
    std::copy(uRef_,uRef_+du*Ns,uRef);
}

/**
** Initialise members
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::init(){
    assert(Tpred>0.);
    assert(Ns_>=1);
    assert(dx>0);
	assert(du>0);
    assert(da>0);
	assert(dy>0);
    assert(dxu>0);
    assert(dxux>0);

    erkInteg.setTimeInterval(Tpred/Ns_);

    /* Primal dimensions and indices */
    izTerm = Ns_*dxu;
    igTerm = (Ns_-1)*dxux+dxu;
    dZshoot = izTerm+dx; // (s_0',q_0',s_1',q_1',...,s_{N-1}',q_{N-1}',s_N')'
    dZpath = Ns_*dpc; // (r_0',...,r_{N-1}')'
    dZ = dZshoot+dZpath;

    /* Dual dimensions */
    dLshoot = Ns_*dx;
    dLpath = Ns_*dpc;
#ifndef BIOREC 
    dL = dx+dLshoot+dLpath;
#else 
    /* Bioreactor example */
#ifndef PEROCP
    /* For real-time enmpc on bioreactor */
    dL = dx+dLshoot+dper+dLpath;
#else
    /* For bioreactor traj generation: need to set x4(0) & x5(0) to 0 */
    dL = 2+dLshoot+dper+dLpath;
#endif
#endif

    /* Dimension of shooting part of AL gradient */
    dGshoot = dZshoot;
    dGshootDBL = dGshoot*SZDBL;
    dGpathDBL = dZpath*SZDBL;

	/* Dimension of vector of integrated augmented states */
    dAugS = Ns_*dw;
    /* Dimension of vector of integrated augmented adjoints */
    dAugA = Ns_*da;

	/* Vector of integrated augmented states */
    if (augStates != NULL)
        delete [] augStates;
    augStates = new double[dAugS];

    /* Vector of integrated augmented adjoints */
    if (augAdjoints != NULL)
        delete [] augAdjoints;
    augAdjoints = new double[da];

    /* Vector of final augmented adjoints */
    if (adfin != NULL)
        delete [] adfin;
    adfin = new double[Ns_*da];

    /* State reference vector */
    if (yRef != NULL)
        delete [] yRef;
    yRef = new double[Ns*dy];

    /* Input reference vector */
    if (uRef != NULL)
        delete [] uRef;
    uRef = new double[Ns*du];

    /**
    ** Initialize integrator 
    **/   
    /* Set inverse input scaling */
    erkInteg.setInputScaling(Dq);
    /* Set pointer to initial augmented-state */
    erkInteg.setStateIni(augStateIni);
    if (erkInteg.getIntegrationOrder()<=0){ 
    	std::cerr<<"ERROR: integrator data not properly set."<<std::endl;
        return;
    }
    /* Initialize and allocate erkInteg members */   
    erkInteg.init();
    if (erkInteg.getDx() != dw){
        std::cerr<<"ERROR<MuShoot>: integrated state dimension different from augmented state dimension."<<std::endl;
        return;
    }
    if (erkInteg.getNsteps() <= 0){
        std::cerr<<"ERROR<MuShoot>: # integration steps nonpositive."<<std::endl;
        return;
    }
    if (erkInteg.getIntegrationOrder() <= 0){
        std::cerr<<"ERROR<MuShoot>: integration order nonpositive."<<std::endl;
        return;
    }

    /* Set state scaling for Mayer term */
    mayer.setStateScaling(Ds);
    mayer.setStateCentering(cs);

    /**
    ** Allocate arrays for intermediate & interpolated augmented states
    **/
    dwnsteps = (erkInteg.getNsteps()+1)*dw;
    sdwnsteps = erkInteg.getIntegrationOrder()*erkInteg.getNsteps()*dw;

    /* Intermediate augmented states */
    if (integAugStates != NULL)
        delete [] integAugStates;
    integAugStates = new double[Ns_*dwnsteps];
    memset(integAugStates,0,Ns_*dwnsteps*SZDBL);

    /* Interpolated states */
    if (interpAugStates != NULL)
        delete [] interpAugStates;
    interpAugStates = new double[Ns_*sdwnsteps];
    memset(interpAugStates,0,Ns_*sdwnsteps*SZDBL);
}

/**
** Compute state scaling (transformation to bring state variables in [0,1]) and adapt problem bounds
**  s = cs + Ds*\tilde{s}
** with \tilde{s} in [0,1]
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::computeStateScaling(){
    
    double xxmm,xxMM;
    double *xxm,*xxM;
    double *ccx,*ddx,*iddx;

    xxm = xmin;
    xxM = xmax;
    ccx = cs;
    ddx = Ds;
    iddx = iDs;
    while (xxm != xmin+dx){
        xxmm = *xxm;
        xxMM = *xxM;
        if ((xxmm<=NEGINF)||(xxMM>=POSINF)){
            /* No scaling if one of the bounds is infinity */
            *ddx++ = 1.;
            *iddx++ = 1.;
            *ccx++ = 0.;
        }
        else {
            /* s = xmin + ds*\tilde{s} */
            /* Scaling ds */
            *ddx++ = (xxmm==xxMM) ? 1. : (xxMM-xxmm);
            /* Inverse scaling 1/ds */
            *iddx++ = (xxmm==xxMM) ? 1. : 1./(xxMM-xxmm);
            /* Centering cs */
            *ccx++ = xxmm;
            /* Update bounds to [0,1] */
            *xxm = 0.;
            *xxM = 1.;
        }
        ++xxm;
        ++xxM;
    }
}

/**
** Compute input scaling (cq,Dq) such that
** q = cq + Dq*\tilde{q} 
** with \tilde{q} in [0,1]
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::computeInputScaling(){
    
    double uumm,uuMM;
    double *uum,*uuM;
    double *ccu,*ddu,*iddu;

    uum = umin;
    uuM = umax;
    ccu = cq;
    ddu = Dq;
    iddu = iDq;
    while (uum != umin+du){
        uumm = *uum;
        uuMM = *uuM;
        if ((uumm<=NEGINF)||(uuMM>=POSINF)){
            /* No scaling if one of the bounds is infinity */
            *ddu++ = 1.;
            *iddu++ = 1.;
            *ccu++ = 0.;
        }
        else {
            /* q = umin + dq*\tilde{q} */
            /* Scaling dq */
            *ddu++ = (uumm==uuMM) ? 1. : (uuMM-uumm);
            /* Inverse scaling 1/dq */
            *iddu++ = (uumm==uuMM) ? 1. : 1./(uuMM-uumm);
            /* Centering cu */ 
            *ccu++ = uumm;
            /* Update bounds */
            *uum = 0.;
            *uuM = 1.;
        }
        ++uum;
        ++uuM;
    }
}

/**
** Apply inverse state scaling from [0,1] to [xmin,xmax]: s = xmin+(xmax-xmin)*\tilde{s}
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::applyInvStateScaling(){
    
    double *xx,*ccs,*dds;

    xx = augStateIni;
    ccs = cs;
    dds = Ds;
    while (xx != augStateIni+dx){
        *xx *= *dds++;
        *xx++ += *ccs++;
    }
}

/**
** Apply inverse input scaling from [0,1] to [umin,umax]: q = umin+(umax-umin)*\tilde{q}
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::applyInvInputScaling(double* un){

    double *us,*uu,*ccq,*ddq;
        
    uu = un;
    us = usc;
    ccq = cq;
    ddq = Dq;
    while (uu != un+du)
        *us++ = (*ccq++)+*ddq++*(*uu++);
}

/********************************************************************************************************/
/************************** Shooting from candidate zS **************************************************/
/********************************************************************************************************/
/**
* Multiple-shooting loop to get augmented states & adjoints from candidate zS
* zS: candidate (s_0',q_0',s_1',q_1',...,s_{N-1}',q_{N-1}',s_N',r_0',...r_{N-1}')'
* s_n: state node
* q_n: input node
* r_n: slack associated with path-constraint at node n 
* mu: dual optimizer (initial-value embedding + shooting constraints + path-constraints)
* halfRho: halved penalty parameter (AL formulation)
* vcan: AL value at candidate zS
* gALshootS, gALpathS
*/
template <class dynT,class costT,class pconT,class mayerT> double MuShoot<dynT,costT,pconT,mayerT>::shootFromCandidate(double halfRho){

    int i,n;
    double *xn,*un;
    double *xf,*xff;
    double *ee,*es;
    double *ms;
    double *yR,*uR;
    double *adjs,*adjsw,*adfinS;
    double *gALsS,*gALsS_,*gALsTS;
    double *gMay;
    double *egaus,*erpaus;
    double e,r2e,mr2e;
    double may,vcan; 
#ifdef BIOREC
    double *dds;
    double ddds;
    double *sini,*sfin;
    double *epr;
    double *mpr;
#endif
#if DPC
    double *rp,*ep,*mp;
    double *gALpS;
#endif
    /**
    ** Initialize pointers
    **/
    /* Initial shooting state s_n */
    xn = zS;
    /* Shooting input q_n */
    un = xn+dx;
    /* Final shooting state w(t_{n+1};s_n,q_n) (integrator output) */
    xf = augStates;
    /* Final shooting state s_{n+1} (node, shooting constraint xff-xf=0) */ 
    xff = zS+dxu;
    /* Intermediate augmented states */
    egaus = integAugStates;
    /* Interpolated augmented states */
    erpaus = interpAugStates;

    /* Output reference */
    yR = yRef; 
    /* Input reference */
    uR = uRef; 
    
    /* Pointer to initial-value embedding constraint */
    ee = econ;
    /* Pointer to shooting constraints */
    es = ee+dx;
#ifdef BIOREC
    /* Pointer to  constraint */
    epr = es+dLshoot;
#endif

    /* Pointer to shooting duals */
    ms = mu+dx; 
#ifdef BIOREC
    /* Pointer to periodicity dual */
    mpr = ms+dLshoot;
#endif 

    /* Pointer to initial shooting augmented adjoint  */
    adjs = augAdjoints;
    /* Current final augmented adjoint (shooting part) */
    adfinS = adfin; 
    /* Pointer to begin shooting node in separated AL gradient */
    gALsS = gALshootS;
    /* Pointer to end shooting node in separated AL gradient */
    gALsS_ = gALshootS+dxu; 
    /* Pointer to terminal-state part of separated AL gradient (shooting) */ 
    gALsTS = gALshootS+igTerm;
#if DPC
    /* Path-constraints slack r_n */
    rp = zS+dZshoot;
    /* Pointer to slacked path-constraints */
    ep = es+dLshoot;
    /* Pointer to path duals */
    mp = ms+dLshoot;
    /* Path-constraint part of AL gradient */
    gALpS = gALpathS; 
#endif
    /**
    ** Loop over shooting intervals [t_n,t_{n+1}], n=0,...,Ns_
    **/
    vcan = 0.;
    for (n=0;n<Ns_;++n){
        /* Set parameters for integrator (input,reference,initial state) */
        erkInteg.setOutputRef(yR); // Set pointer to output-reference
        erkInteg.setInputRef(uR); // Set pointer to input-reference
#if UNSCAL_VAR
        erkInteg.setInput(un); // Set pointer to input
#endif
#if SCAL_VAR
        erkInteg.setInput(usc); // Set pointer to scaled input
        applyInvInputScaling(un); // Scale input from [0,1] to [umin,umax]
#endif
        memcpy(augStateIni,xn,dxbl); // Copy initial state into initial augmented-state
#if SCAL_VAR
        applyInvStateScaling();
#endif  
        erkInteg.setStateFin(xf); // Set pointer to final augmented-state
        erkInteg.setStates(egaus); // Set pointer to intermediate augmented states 
        erkInteg.setStatesInterp(erpaus); // Set pointer to interpolated augmented states

        /* Augmented-state forward integration */
        erkInteg.integrateStateFixed(); // Forward integration (augmented state): get final shooting state and cost
        
        /* Add stage-cost */
        vcan += *(xf+dx);

/*        std::cout<<"System input n= "<<n<<std::endl;
        for (i=0;i<du;++i)
            std::cout<<un[i]<<std::endl;
        std::cout<<"Integrator input-output n= "<<n<<std::endl;
        for (i=0;i<dw;++i)
            std::cout<<augStateIni[i]<<", "<<xf[i]<<std::endl; */

        /* Compute shooting equality constraint, final adjoint and partially separable AL gradient wrt s_{n+1} */
        for (i=0;i<dx;++i){
#if UNSCAL_VAR
            e = *(xff+i)-*(xf+i);
#endif
#if SCAL_VAR
            e = *(xff+i)+*(iDs+i)*(*(cs+i)-*(xf+i));
#endif
            r2e = halfRho*e;
            mr2e = *(ms+i)+r2e;
            vcan += mr2e*e;
            mr2e += r2e;
            /* Store in shooting constraint */
            *es++ = e;
            /* Initialize final augmented-adjoint (state part) */
#if UNSCAL_VAR
            *(adfinS+i) = -mr2e;
#endif
#if SCAL_VAR
            *(adfinS+i) = -*(iDs+i)*mr2e;
#endif
            /* AL gradient wrt s_{i+1} */
            *(gALsS_+i) = mr2e;
        }
        /* Set final augmented-adjoint */
        *(adfinS+dx) = 1.; // Final cost-adjoint
        memset(adfinS+dw,0,dubl); // Final input-adjoint
        
        /* Backwards state interpolation over shooting interval: get intermediate states for backwards integration */
        erkInteg.interpState();
        
        /* Augmented-adjoint backward integration */
        erkInteg.setAdjointFin(adfinS); // Set final adjoint
        erkInteg.setAdjointIni(adjs); // Set initial adjoint (output of backwards integration)
        erkInteg.integrateAdjointFixed(); // Backwards-in-time integration (augmented adjoint): get initial shooting adjoint
        
        /* AL gradient wrt s_n */
#if UNSCAL_VAR
        memcpy(gALsS,adjs,dxbl);
#endif
#if SCAL_VAR
        for (i=0;i<dx;++i){
            *(gALsS+i) = *(Ds+i)**(adjs+i);
        }
#endif
        /* AL gradient wrt q_n */
        memcpy(gALsS+dx,adjs+dw,dubl);
#if DPC
        /* Evaluate path-constraints */
        pcon.eval(econPath,xn,un,rp);
        /* Compute AL relaxation of path-constraint */
        for (i=0;i<dpc;++i){
            e = *(econPath+i);
            r2e = halfRho*e;
            mr2e = *(mp+i)+r2e;
            vcan += mr2e*e;
            mr2e += r2e;
            /* Store in path-constraints */
            *ep++ = e;
            /* Store (mu+rho*e) */
            *(gALpS+i) = mr2e;
        }
        /* Evaluate path-constraints jacobians transpose */
        pcon.evalJacXt(xn,un); // wrt state
        pcon.evalJacUt(xn,un); // wrt input
        /* AL gradient wrt s_n, q_n */
        pcon.prdJacXt(gALsS,gALpS); // Product against state-jacobian transpose
        pcon.prdJacUt(gALsS+dx,gALpS); // Product against input jacobian transpose
        /* Update path pointers */
        rp += dpc;
        mp += dpc;
        gALpS += dpc;
#endif
        /* Update shooting and reference pointers */
        xn += dxu;
        un += dxu;
        xf += dw;
        xff += dxu;
        egaus += dwnsteps;
        erpaus += sdwnsteps;
        yR += dy;
        uR += du;
        ms += dx;
        adfinS += da;
        gALsS += dxux;
        gALsS_ += dxux;
    }
    /**
    * Constraint on initial shooting node, s_0-xest=0 
    */
    for (i=0;i<dx;++i){
#if UNSCAL_VAR
        e = *(zS+i)-*(xest+i);
#endif
#if SCAL_VAR
        e = *(cs+i)+*(Ds+i)**(zS+i)-*(xest+i);
#endif
        r2e = halfRho*e;
        mr2e = *(mu+i)+r2e;
        vcan += mr2e*e;
        mr2e += r2e;
        /* Store in initial-value embedding */
        *ee++ = e;
        /* AL gradient wrt s_0 */
#if UNSCAL_VAR
        *(gALshootS+i) += mr2e;
#endif
#if SCAL_VAR
        *(gALshootS+i) += *(Ds+i)*mr2e;
#endif
    }   
#ifdef BIOREC
    /**
    ** Periodicity constraint 
    **/
    /* Pointer to begin shooting node in separated AL gradient */
    gALsS = gALshootS; 
    /* Pointer to terminal-state part of separated AL gradient (shooting) */ 
    gALsTS = gALshootS+igTerm;
    /* Pointer to initial shooting node s_0 */
    sini = zS;
    /* Pointer to final shooting node s_N */
    sfin = sini+izTerm;
    /* Pointer to state scaling */
    dds = Ds;
    while (epr != econ+dL){
#if SCAL_VAR
        ddds = *dds++;
        e = ddds*(*sini++-*sfin++); 
#endif
#if UNSCAL_VAR
        e = *sini++-*sfin++;
#endif
        r2e = halfRho*e;
        mr2e = (*mpr++)+r2e;
        vcan += mr2e*e;
        mr2e += r2e;
        /* Store periodicity constraint */
        *epr++ = e;
#if SCAL_VAR
        /* AL gradient wrt s_0 */
        *gALsS++ += ddds*mr2e;
        /* AL gradient wrt s_N */
        *gALsTS++ -= ddds*mr2e;
#endif
#if UNSCAL_VAR
        /* AL gradient wrt s_0 */
        *gALsS++ += mr2e;
        /* AL gradient wrt s_N */
        *gALsTS++ -= mr2e;
#endif
    }
#endif
#if MAY
    /**
    ** Mayer term
    **/
    /* Mayer objective */
    may = mayer.eval(zS+izTerm);
    vcan += may; 
    /* Mayer gradient */
    mayer.grad(zS+izTerm);   
    gMay = mayer.getGrad(); 
    for (i=0;i<dx;++i){
#if UNSCAL_VAR
        *(gALsTS+i) += *(gMay+i);
#endif
#if SCAL_VAR
        *(gALsTS+i) += *(Ds+i)**(gMay+i);
#endif
    }
#endif

    return vcan;
}

/**
** Multiple-shooting loop to get augmented states and adjoints from candidate zS with periodicity constraint s_0-s_N=0
** and without initial constraint s_0-xest=0 (except on 2 auxiliary states)
**/
template <class dynT,class costT,class pconT,class mayerT> double MuShoot<dynT,costT,pconT,mayerT>::shootFromCandidatePeriodic(double halfRho){
    
//    std::cout<<"Shoot from candidate periodic"<<std::endl;

    int i,n;
    double *xn,*un;
    double *xf,*xff;
    double *es;
    double *ms;
    double *yR,*uR;
    double *adjs,*adjsw,*adfinS;
    double *gALsS,*gALsS_,*gALsTS;
    double *gMay;
    double *egaus,*erpaus;
    double e,r2e,mr2e;
    double may,vcan;
    /* Pointers for periodicity and initial constraint */ 
    double *sini,*sfin;
    double *epr;
    double *mpr;
    double *ein;
    double *min;
#if DPC
    double *rp,*ep,*mp;
    double *gALpS;
#endif
    /**
    ** Initialize pointers
    **/
    /* Initial shooting state s_n */
    xn = zS;
    /* Shooting input q_n */
    un = xn+dx;
    /* Final shooting state w(t_{n+1};s_n,q_n) (integrator output) */
    xf = augStates;
    /* Final shooting state s_{n+1} (node, shooting constraint xff-xf=0) */ 
    xff = zS+dxu;
    /* Intermediate augmented states */
    egaus = integAugStates;
    /* Interpolated augmented states */
    erpaus = interpAugStates;

    /* Output reference */
    yR = yRef; 
    /* Input reference */
    uR = uRef; 
    
    ein = econ;
    /* Pointer to shooting constraints */
    es = ein+2;
    /* Pointer to periodicity constraint */
    epr = es+dLshoot;

    min = mu;
    /* Pointer to shooting duals */
    ms = min+2; 
    /* Pointer to periodicity dual */
    mpr = ms+dLshoot;

    /* Pointer to initial shooting augmented adjoint  */
    adjs = augAdjoints;
    /* Current final augmented adjoint (shooting part) */
    adfinS = adfin; 
    /* Pointer to begin shooting node in separated AL gradient */
    gALsS = gALshootS;
    /* Pointer to end shooting node in separated AL gradient */
    gALsS_ = gALshootS+dxu; 
    /* Pointer to terminal-state part of separated AL gradient (shooting) */ 
    gALsTS = gALshootS+igTerm;

#if DPC
    /* Path-constraints slack r_n */
    rp = zS+dZshoot;
    /* Pointer to slacked path-constraints */
    ep = epr+dper;
    /* Pointer to path duals */
    mp = mpr+dper;
    /* Path-constraint part of AL gradient */
    gALpS = gALpathS; 
#endif
    /**
    ** Loop over shooting intervals [t_n,t_{n+1}], n=0,...,Ns_
    **/
    vcan = 0.;
    for (n=0;n<Ns_;++n){
        /* Set parameters for integrator (input,reference,initial state) */
        erkInteg.setOutputRef(yR); // Set pointer to output-reference
        erkInteg.setInputRef(uR); // Set pointer to input-reference
#if UNSCAL_VAR
        erkInteg.setInput(un); // Set pointer to input
#endif
#if SCAL_VAR
        erkInteg.setInput(usc); // Set pointer to scaled input
        applyInvInputScaling(un); // Scale input from [0,1] to [umin,umax]
#endif
        memcpy(augStateIni,xn,dxbl); // Copy initial state into initial augmented-state
#if SCAL_VAR
        applyInvStateScaling();
#endif  
        erkInteg.setStateFin(xf); // Set pointer to final augmented-state
        erkInteg.setStates(egaus); // Set pointer to intermediate augmented states 
        erkInteg.setStatesInterp(erpaus); // Set pointer to interpolated augmented states

        /* Augmented-state forward integration */
        erkInteg.integrateStateFixed(); // Forward integration (augmented state): get final shooting state and cost
        
        /* Add stage-cost */
        vcan += *(xf+dx);

/*        std::cout<<"System input n= "<<n<<std::endl;
        for (i=0;i<du;++i)
            std::cout<<un[i]<<std::endl;
        std::cout<<"Integrator input-output n= "<<n<<std::endl;
        for (i=0;i<dw;++i)
            std::cout<<augStateIni[i]<<", "<<xf[i]<<std::endl; */

        /* Compute shooting equality constraint, final adjoint and partially separable AL gradient wrt s_{n+1} */
        for (i=0;i<dx;++i){
#if UNSCAL_VAR
            e = *(xff+i)-*(xf+i);
#endif
#if SCAL_VAR
            e = *(xff+i)+*(iDs+i)*(*(cs+i)-*(xf+i));
#endif
            r2e = halfRho*e;
            mr2e = *(ms+i)+r2e;
            vcan += mr2e*e;
            mr2e += r2e;
            /* Store in shooting constraint */
            *es++ = e;
            /* Initialize final augmented-adjoint (state part) */
#if UNSCAL_VAR
            *(adfinS+i) = -mr2e;
#endif
#if SCAL_VAR
            *(adfinS+i) = -*(iDs+i)*mr2e;
#endif
            /* AL gradient wrt s_{i+1} */
            *(gALsS_+i) = mr2e;
        }
        /* Set final augmented-adjoint */
        *(adfinS+dx) = 1.; // Final cost-adjoint
        memset(adfinS+dw,0,dubl); // Final input-adjoint
        
        /* Backwards state interpolation over shooting interval: get intermediate states for backwards integration */
        erkInteg.interpState();
        
        /* Augmented-adjoint backward integration */
        erkInteg.setAdjointFin(adfinS); // Set final adjoint
        erkInteg.setAdjointIni(adjs); // Set initial adjoint (output of backwards integration)
        erkInteg.integrateAdjointFixed(); // Backwards-in-time integration (augmented adjoint): get initial shooting adjoint
        
        /* AL gradient wrt s_n */
#if UNSCAL_VAR
        memcpy(gALsS,adjs,dxbl);
#endif
#if SCAL_VAR
        for (i=0;i<dx;++i){
            *(gALsS+i) = *(Ds+i)**(adjs+i);
        }
#endif
        /* AL gradient wrt q_n */
        memcpy(gALsS+dx,adjs+dw,dubl);
#if DPC
        /* Evaluate path-constraints */
        pcon.eval(econPath,xn,un,rp);
        /* Compute AL relaxation of path-constraint */
        for (i=0;i<dpc;++i){
            e = *(econPath+i);
            r2e = halfRho*e;
            mr2e = *(mp+i)+r2e;
            vcan += mr2e*e;
            mr2e += r2e;
            /* Store in path-constraints */
            *ep++ = e;
            /* Store (mu+rho*e) */
            *(gALpS+i) = mr2e;
        }
        /* Evaluate path-constraints jacobians transpose */
        pcon.evalJacXt(xn,un); // wrt state
        pcon.evalJacUt(xn,un); // wrt input
        /* AL gradient wrt s_n, q_n */
        pcon.prdJacXt(gALsS,gALpS); // Product against state-jacobian transpose
        pcon.prdJacUt(gALsS+dx,gALpS); // Product against input jacobian transpose
        /* Update path pointers */
        rp += dpc;
        mp += dpc;
        gALpS += dpc;
#endif
        /* Update shooting and reference pointers */
        xn += dxu;
        un += dxu;
        xf += dw;
        xff += dxu;
        egaus += dwnsteps;
        erpaus += sdwnsteps;
        yR += dy;
        uR += du;
        ms += dx;
        adfinS += da;
        gALsS += dxux;
        gALsS_ += dxux;
    }
    /**
    ** Initial condition (on last 2 states of bioreactor)
    ** x(4) = 0
    ** x(5) = 0
    **/
    sini = zS+3;
    gALsS = gALshootS+3;
    while (ein != econ+2){
        e = *sini++;
        r2e = halfRho*e;
        mr2e = (*min++)+r2e;
        vcan += mr2e*e;
        mr2e += r2e; 
        /* Store initial constraint */
        *ein++ = e;
        /* Gradient wrt s_0 */
        *gALsS++ += mr2e;
    }
    /**
    ** Periodicity constraint on first 3 states
    **/
    /* Pointer to begin shooting node in separated AL gradient */
    gALsS = gALshootS; 
    /* Pointer to terminal-state part of separated AL gradient (shooting) */ 
    gALsTS = gALshootS+igTerm;
    /* Pointer to initial shooting node s_0 */
    sini = zS;
    /* Pointer to final shooting node s_N */
    sfin = sini+izTerm;
    while (epr != econ+dL){
        e = *sini-*sfin;
        ++sini;
        ++sfin;
        r2e = halfRho*e;
        mr2e = (*mpr++)+r2e;
        vcan += mr2e*e;
        mr2e += r2e;
        /* Store periodicity constraint */
        *epr++ = e;
        /* AL gradient wrt s_0 */
        *gALsS++ += mr2e;
        /* AL gradient wrt s_N */
        *gALsTS++ -= mr2e;
    }
#if MAY
    /**
    * Mayer term
    */
    /* Mayer objective */
    may = mayer.eval(zS+izTerm);
    vcan += may; 
    /* Mayer gradient */
    mayer.grad(zS+izTerm);   
    gMay = mayer.getGrad(); 
    for (i=0;i<dx;++i){
#if UNSCAL_VAR
        *(gALsTS+i) += *(gMay+i);
#endif
#if SCAL_VAR
        *(gALsTS+i) += *(Ds+i)**(gMay+i);
#endif
    }
#endif

    return vcan;
}

/**
* Multiple-shooting loop to get augmented states, augmented Lagrangian and part of the gradient of the augmented Lagrangian
* zS: candidate (s_0',q_0',s_1',q_1',...,s_{N-1}',q_{N-1}',s_N',r_0',...r_{N-1}')'
* s_n: state node
* q_n: input node
* r_n: slack associated with path-constraint at node n 
* mu: dual optimizer (initial-value embedding + shooting constraints + path-constraints)
* halfRho: halved penalty parameter (AL formulation)
* vcan: AL value at candidate
* gALshootS, gALpathS: fill only the s_{n+1} parts (available from shooting constraints evaluation) and s_0
*/
template <class dynT,class costT,class pconT,class mayerT> double MuShoot<dynT,costT,pconT,mayerT>::shootStatesFromCandidate(double halfRho){

    int i,n;
    double *xn,*un;
    double *xf,*xff;
    double *ee,*es;
    double *ms;
    double *yR,*uR;
    double *adfinS;
    double *gALsS;
    double *egaus,*erpaus;
    double e,r2e,mr2e;
    double may,vcan;
#if DPC
    double *rp,*ep,*mp;
    double *gALpS;
#endif

//    std::cout<<"Shoot sta from candidate"<<std::endl;

    /**
    * Initialize pointers
    */
    /* Initial shooting state s_n */
    xn = zS;
    /* Shooting input q_n */
    un = xn+dx;
    /* Final shooting state w(t_{n+1};s_n,q_n) (integrator output) */
    xf = augStates;
    /* Final shooting state s_{n+1} (node, shooting constraint xff-xf=0) */ 
    xff = zS+dxu;
    /* Intermediate augmented states */
    egaus = integAugStates;
    /* Interpolated augmented states */
    erpaus = interpAugStates;
	/* Output reference */
    yR = yRef; 
    /* Input reference */
    uR = uRef; 
    /* Pointer to initial-value embedding constraint */
    ee = econ;
    /* Pointer to shooting constraints */
    es = ee+dx;
    /* Pointer to shooting duals */
    ms = mu+dx;  
    /* Current final augmented adjoint (shooting part) */
    adfinS = adfin; 
    /* Current AL gradient copy (shooting part) */
    gALsS = gALshootS+dxu;

#if DPC 
    /* Path-constraints slack r_n */
    rp = zS+dZshoot;   
    /* Pointer to slacked path-constraints */
    ep = es+dLshoot;
    /* Pointer to path duals */
    mp = ms+dLshoot;
    /* Path-constraint part of AL gradient */
    gALpS = gALpathS;
#endif
    /**
    * Iterate over shooting intervals, s_{n+1}-x(t_{n+1};s_n,q_n) = 0 with augmented-state integration on [t_n,t_{n+1}]
    */
    vcan = 0.;
    for (n=0;n<Ns_;++n){
        /* Set parameters for integrator (input,reference,initial state) */
        erkInteg.setOutputRef(yR); // Set pointer to output-reference
        erkInteg.setInputRef(uR); // Set pointer to input-reference
        erkInteg.setInput(un); // Set pointer to input
        memcpy(augStateIni,xn,dxbl); // Copy initial state into initial augmented-state
        erkInteg.setStateFin(xf); // Set pointer to final augmented-state
        erkInteg.setStates(egaus); // Set pointer to intermediate augmented states 
        erkInteg.setStatesInterp(erpaus); // Set pointer to interpolated augmented states
#if SCAL_VAR
        /* Apply inverse state scaling*/
        applyInvStateScaling();
        /* Apply inverse input scaling */
        erkInteg.applyInvInputScaling();
#endif
        /* Augmented-state forward integration */
        erkInteg.integrateStateFixed(); // Forward integration (augmented state): get final shooting state and cost

        /* Add stage-cost */
        vcan += *(xf+dx);

        /* Compute shooting equality constraint, final adjoint and partially separable AL gradient wrt s_{n+1} */
        for (i=0;i<dx;++i){
#if UNSCAL_VAR
            e = *(xff+i)-*(xf+i);
#endif
#if SCAL_VAR
            e = *(xff+i)-*(cs+i)-*(Ds+i)**(xf+i);
#endif
            r2e = halfRho*e;
            mr2e = *(ms+i)+r2e;
            vcan += mr2e*e;
            mr2e += r2e;
            /* Store in shooting constraint */
            *es++ = e;
            /* Initialize final augmented-adjoint (state part) */
            *(adfinS+i) = -mr2e;
            /* AL gradient wrt s_{i+1} */
#if UNSCAL_VAR
            *(gALsS+i) = mr2e;
#endif
#if SCAL_VAR
            *(gALsS+i) = *(iDs+i)*mr2e;
#endif
        }

#if DPC
        /* Evaluate path-constraints */
        pcon.eval(econPath,xn,un,rp);
        
        for (i=0;i<dpc;++i){ 
            e = *(econPath+i);
            r2e = halfRho*e;
            mr2e = *(mp+i)+r2e;
            vcan += mr2e*e;
            mr2e += r2e;
            /* Store in path-constraints */
            *ep++ = e;
            /* Store (mu+rho*e) */
            *(gALpS+i) = mr2e;
        }
        /* Update path pointers */
        rp += dpc;
        mp += dpc;
        gALpS += dpc;
#endif
        /* Update shooting and reference pointers */
        xn += dxu;
        un += dxu;
        xf += dw;
        xff += dxu;
        egaus += dwnsteps;
        erpaus += sdwnsteps;
        yR += dy;
        uR += du;
        ms += dx;
        adfinS += da;
        gALsS += dxux;
    }
    /**
    * Constraint on initial shooting node, s_0-xest = 0
    */
    for (i=0;i<dx;++i){
        e = *(zS+i)-*(xest+i);
        r2e = halfRho*e;
        mr2e = *(mu+i)+r2e;
        vcan += mr2e*e;
        mr2e += r2e;
        /* Store in initial-value embedding */
        *ee++ = e;
        /* AL gradient wrt s_0 */
        *(gALshootS+i) = mr2e;
    }
#if MAY
    /**
    * Mayer term 
    */
    may = mayer.eval(zS+izTerm);
    vcan += may; 
#endif
    return vcan;
}

/**
* Multiple-shooting loop to get augmented adjoints and rest of gradient of augmented Lagrangian
* zS: candidate for primal optimizer (s_0',q_0',s_1',q_1',...,s_{N-1}',q_{N-1}',s_N',r_0',...r_{N-1}')'
* s_n: state node
* q_n: input node
* r_n: slack associated with path-constraint at node n 
* mu: dual optimizer (initial-value embedding + shooting constraints + path-constraints)
* gALshootS
*/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::shootAdjointsFromCandidate(){

 
}

/*******************************************************************************************************/
/****************************************** Shooting from iterate z *************************************/
/*******************************************************************************************************/
/**
* Multiple-shooting loop to get augmented states & adjoints from iterate z
*/
template <class dynT,class costT,class pconT,class mayerT> double MuShoot<dynT,costT,pconT,mayerT>::shootFromIterate(double halfRho){ 
        
    int i,n;
    double *xn,*un;
    double *xf,*xff;
    double *ee,*es;
    double *ms;
    double *yR,*uR;
    double *adjs,*adfinS;
    double *gALs,*gALs_,*gALsT,*gALt;
    double *g,*gMay;
    double *egaus,*erpaus;
    double e,r2e,mr2e;
    double may,vAL,gm;
#if BIOREC
    double *dds;
    double ddds;
    double *sini,*sfin;
    double *epr;
    double *mpr;
#endif 
#if DPC
    double *rp,*ep,*mp;
    double *gALp;
#endif
    /**
    ** Initialize pointers
    **/
    /* Initial shooting state s_n */
    xn = z;
    /* Shooting input q_n */
    un = xn+dx;
    /* Final shooting state w(t_{n+1};s_n,q_n) (integrator output) */
    xf = augStates;
    /* Final shooting state s_{n+1} (node, shooting constraint xff-xf=0) */ 
    xff = z+dxu;
    /* Intermediate augmented states */
    egaus = integAugStates;
    /* Interpolated augmented states */
    erpaus = interpAugStates;
    /* Output reference */
    yR = yRef; 
    /* Input reference */
    uR = uRef; 
    
    /* Pointer to initial-value embedding constraint */
    ee = econ;
    /* Pointer to shooting constraints */
    es = ee+dx;
#if BIOREC
    /* Pointer to periodicity constraint */
    epr = es+dLshoot;
#endif
    
    /* Pointer to shooting duals */
    ms = mu+dx; 
#if BIOREC
    /* Pointer to periodicity dual */
    mpr = mu+dLshoot;
#endif

    /* Pointer to initial shooting augmented adjoint  */
    adjs = augAdjoints;
    /* Current final augmented adjoint (shooting part) */
    adfinS = adfin; 

    /* Pointer to full AL gradient */
    g = gAL;
    /* Pointer to begin shooting node in separated AL gradient */
    gALs = gALshoot;
    /* Pointer to end shooting node in separated AL gradient */
    gALs_ = gALshoot+dxu;
    /* Pointer to terminal-state part of separated AL gradient (shooting) */ 
    gALsT = gALs+igTerm; 
    /* Pointer to terminal-state part of full AL gradient */ 
    gALt = gAL+izTerm;
#if DPC
    /* Path-constraints slack r_n */
    rp = z+dZshoot;
    /* Pointer to slacked path-constraints */
    ep = es+dLshoot;
    /* Pointer to path duals */
    mp = ms+dLshoot; 
    /* Path-constraint part of AL gradient */
    gALp = gALpath; 
#endif

    /* Set shooting part of full AL gradient to 0 */
    memset(gAL,0,dGshootDBL);

    /** Loop over shooting intervals [t_n,t_{n+1}], n=0,...,Ns_ **/
    vAL = 0.;
    for (n=0;n<Ns_;++n){
        /* Set parameters for integrator (input,reference,initial state) */
        erkInteg.setOutputRef(yR); // Set pointer to output-reference
        erkInteg.setInputRef(uR); // Set pointer to input-reference
#if UNSCAL_VAR
        erkInteg.setInput(un); // Set pointer to input
#endif
#if SCAL_VAR
        applyInvInputScaling(un);
        erkInteg.setInput(usc); // Set pointer to scaled input 
#endif
        memcpy(augStateIni,xn,dxbl); // Copy initial state into initial augmented-state
#if SCAL_VAR
        applyInvStateScaling();
#endif
        erkInteg.setStateFin(xf); // Set pointer to final augmented-state
        erkInteg.setStates(egaus); // Set pointer to intermediate augmented states 
        erkInteg.setStatesInterp(erpaus); // Set pointer to interpolated augmented states

        /* Augmented-state forward integration */
//        std::cout<<"Start integration"<<std::endl;
        erkInteg.integrateStateFixed(); // Forward integration (augmented state): get final shooting state and cost
//        std::cout<<"...integration done."<<std::endl;

        /* Add stage-cost */
        vAL += *(xf+dx);

/*        std::cout<<"System input n= "<<n<<std::endl;
        for (i=0;i<du;++i)
            std::cout<<un[i]<<std::endl;
        std::cout<<"Integrator input-output n= "<<n<<std::endl;
        for (i=0;i<dw;++i)
            std::cout<<augStateIni[i]<<", "<<xf[i]<<std::endl; */

        /* Compute shooting equality constraint, final adjoint and partially separable AL gradient wrt s_{n+1} */
        for (i=0;i<dx;++i){
#if UNSCAL_VAR
            e = *(xff+i)-*(xf+i);
#endif
#if SCAL_VAR
            e = *(xff+i)+*(iDs+i)*(*(cs+i)-*(xf+i));
#endif
            r2e = halfRho*e;
            mr2e = *(ms+i)+r2e;
            vAL += mr2e*e;
            mr2e += r2e;
            /* Store in shooting constraint */
            *es++ = e;
            /* Initialize final augmented-adjoint (state part) */
#if UNSCAL_VAR
            *(adfinS+i) = -mr2e;
#endif
#if SCAL_VAR
            *(adfinS+i) = -*(iDs+i)*mr2e;
#endif
            /* AL gradient wrt s_{i+1} */
            *(gALs_+i) = mr2e;
        }

        /* Set final augmented-adjoint */
        *(adfinS+dx) = 1.; // Final cost-adjoint
        memset(adfinS+dw,0,dubl); // Final input-adjoint
        
        /* Backwards state interpolation over shooting interval: get intermediate states for backwards integration */
        erkInteg.interpState();
        
        /* Augmented-adjoint backward integration */
        erkInteg.setAdjointFin(adfinS); // Set final adjoint
        erkInteg.setAdjointIni(adjs); // Set initial adjoint (output of backwards integration)
        erkInteg.integrateAdjointFixed(); // Backwards-in-time integration (augmented adjoint): get initial shooting adjoint
        
        /* AL gradient wrt s_n */
#if UNSCAL_VAR
        memcpy(gALs,adjs,dxbl);
#endif
#if SCAL_VAR
        for (i=0;i<dx;++i)
            *(gALs+i) = *(Ds+i)**(adjs+i);
#endif
        /* AL gradient wrt q_n */
        memcpy(gALs+dx,adjs+dw,dubl);

#if DPC
        /* Evaluate path-constraints */
        pcon.eval(econPath,xn,un,rp);
        
        for (i=0;i<dpc;++i){
            e = *(econPath+i);
            r2e = halfRho*e;
            mr2e = *(mp+i)+r2e;
            vAL += mr2e*e;
            mr2e += r2e;
            /* Store in path-constraints */
            *ep++ = e;
            /* Store (mu+rho*e) */
            *(gALp+i) = mr2e;
        }
        /* Evaluate path-constraints jacobians transpose */
        pcon.evalJacXt(xn,un); // wrt state
        pcon.evalJacUt(xn,un); // wrt input

        /* AL gradient wrt s_n, q_n */
        pcon.prdJacXt(gALs,gALp); // Product against state-jacobian transpose
        pcon.prdJacUt(gALs+dx,gALp); // Product against input jacobian transpose
        
        /* Update path pointers */
        rp += dpc;
        mp += dpc;
        gALp += dpc;
#endif
        /* Build full gradient */
        for (i=0;i<dxux;++i)
            *(g+i) += *gALs++;

        /* Update shooting and reference pointers */
        xn += dxu;
        un += dxu;
        xf += dw;
        xff += dxu;
        egaus += dwnsteps;
        erpaus += sdwnsteps;
        yR += dy;
        uR += du;
        ms += dx;
        adfinS += da;
        g += dxu;
        gALs_ += dxux;
    }
    /**
    * Constraint on initial shooting node, s_0-xest=0 
    */
    for (i=0;i<dx;++i){
#if UNSCAL_VAR
        e = *(z+i)-*(xest+i);
#endif
#if SCAL_VAR
        e = *(Ds+i)**(z+i)+*(cs+i)-*(xest+i);
#endif
        r2e = halfRho*e;
        mr2e = *(mu+i)+r2e;
        vAL += mr2e*e;
        mr2e += r2e;
        /* Store in initial-value embedding */
        *ee++ = e;
        /* AL gradient wrt s_0 */
#if UNSCAL_VAR
        *(gALshoot+i) += mr2e;
        *(gAL+i) += mr2e;
#endif
#if SCAL_VAR
        *(gALshoot+i) += *(Ds+i)*mr2e;
        *(gAL+i) += *(Ds+i)*mr2e;
#endif
    }
#if BIOREC
    /**
    * Periodicity constraint 
    */
    /* Pointer to full AL gradient */
    g = gAL;
    /* Pointer to begin shooting node in separated AL gradient */
    gALs = gALshoot;
    /* Pointer to terminal-state part of full AL gradient */ 
    gALt = gAL+izTerm;
    /* Pointer to terminal-state part of separated AL gradient (shooting) */ 
    gALsT = gALs+igTerm; 
    /* Pointer to intial and final shooting nodes */
    sini = z;
    sfin = sini+izTerm;
    /* Pointer to state scaling */
    dds = Ds;
    while (epr != econ+dL){
#if SCAL_VAR
        ddds = *dds++;
        e = ddds*(*sini++-*sfin++);
#endif
#if UNSCAL_VAR
        e = *sini++-*sfin++;
#endif
        r2e = halfRho*e;
        mr2e = (*mpr++)+r2e;
        vAL += mr2e*e;
        mr2e += r2e;
        /* Store periodicity constraint */
        *epr++ = e;
#if SCAL_VAR
        /* AL gradient wrt s_0 */
        *g++ += ddds*mr2e;
        *gALs++ += ddds*mr2e;
        /* AL gradient wrt s_N */
        *gALt++ -= ddds*mr2e;
        *gALsT++ -= ddds*mr2e; 
#endif
#if UNSCAL_VAR
        /* AL gradient wrt s_0 */
        *g++ += mr2e;
        *gALs++ += mr2e;
        /* AL gradient wrt s_N */
        *gALt++ -= mr2e;
        *gALsT++ -= mr2e; 
#endif
    }
#endif
#if MAY    
    /**
    * Mayer term
    */
    /* Mayer objective */
    may = mayer.eval(z+izTerm);
    vAL += may; 
    /* Mayer gradient */
    mayer.grad(z+izTerm);   
    gMay = mayer.getGrad(); 
    for (i=0;i<dx;++i){
#if UNSCAL_VAR
        gm = *(gMay+i);
#endif
#if SCAL_VAR
        gm = *(Ds+i)**(gMay+i);
#endif
        *(gALsT+i) += gm;
        *(gALt+i) += gm;
    }
#endif
    return vAL;
}

/**
** Multiple-shooting loop to get augmented states and adjoints from iterate z with periodicity constraints,
** build full AL gradient 
**/
template <class dynT,class costT,class pconT,class mayerT> double MuShoot<dynT,costT,pconT,mayerT>::shootFromIteratePeriodic(double halfRho){
    
//    std::cout<<"Shoot from iterate periodic"<<std::endl;

    int i,n;
    double *xn,*un;
    double *xf,*xff;
    double *ee,*es;
    double *ms;
    double *yR,*uR;
    double *adjs,*adfinS;
    double *gALs,*gALs_,*gALsT,*gALt;
    double *g,*gMay;
    double *egaus,*erpaus;
    double e,r2e,mr2e;
    double may,vAL,gm;
    /* Pointers for periodicity */
    double *sini,*sfin;
    double *ein;
    double *min;
    double *epr;
    double *mpr;
#if DPC
    double *rp,*ep,*mp;
    double *gALp;
#endif

    /**
    ** Initialize pointers
    **/
    /* Initial shooting state s_n */
    xn = z;
    /* Shooting input q_n */
    un = xn+dx;
    /* Final shooting state w(t_{n+1};s_n,q_n) (integrator output) */
    xf = augStates;
    /* Final shooting state s_{n+1} (node, shooting constraint xff-xf=0) */ 
    xff = z+dxu;
    /* Intermediate augmented states */
    egaus = integAugStates;
    /* Interpolated augmented states */
    erpaus = interpAugStates;
    /* Output reference */
    yR = yRef; 
    /* Input reference */
    uR = uRef; 
    
    /* Pointer to constraint on initial state (2 aux. states)*/
    ein = econ;
    /* Pointer to shooting constraints */
    es = ein+2;
    /* Pointer to periodicity constraint */
    epr = es+dLshoot;
    
    /* Pointer to dual var associated with initial state (2 aux. states) */
    min = mu;
    /* Pointer to shooting duals */
    ms = min+2; 
    /* Pointer to periodicity dual */
    mpr = ms+dLshoot;

    /* Pointer to initial shooting augmented adjoint  */
    adjs = augAdjoints;
    /* Current final augmented adjoint (shooting part) */
    adfinS = adfin; 

    /* Pointer to full AL gradient */
    g = gAL;
    /* Pointer to begin shooting node in separated AL gradient */
    gALs = gALshoot;
    /* Pointer to end shooting node in separated AL gradient */
    gALs_ = gALshoot+dxu;
    /* Pointer to terminal-state part of separated AL gradient (shooting) */ 
    gALsT = gALs+igTerm; 
    /* Pointer to terminal-state part of full AL gradient */ 
    gALt = gAL+izTerm;
#if DPC
    /* Path-constraints slack r_n */
    rp = z+dZshoot;
    /* Pointer to slacked path-constraints */
    ep = epr+dper;
    /* Pointer to path duals */
    mp = mpr+dper; 
    /* Path-constraint part of AL gradient */
    gALp = gALpath; 
#endif

    /* Set shooting part of full AL gradient to 0 */
    memset(gAL,0,dGshootDBL);

    /** Loop over shooting intervals [t_n,t_{n+1}], n=0,...,Ns_ **/
    vAL = 0.;
    for (n=0;n<Ns_;++n){
        /* Set parameters for integrator (input,reference,initial state) */
        erkInteg.setOutputRef(yR); // Set pointer to output-reference
        erkInteg.setInputRef(uR); // Set pointer to input-reference
#if UNSCAL_VAR
        erkInteg.setInput(un); // Set pointer to input
#endif
#if SCAL_VAR
        applyInvInputScaling(un);
        erkInteg.setInput(usc); // Set pointer to scaled input 
#endif
        memcpy(augStateIni,xn,dxbl); // Copy initial state into initial augmented-state
#if SCAL_VAR
        applyInvStateScaling();
#endif
        erkInteg.setStateFin(xf); // Set pointer to final augmented-state
        erkInteg.setStates(egaus); // Set pointer to intermediate augmented states 
        erkInteg.setStatesInterp(erpaus); // Set pointer to interpolated augmented states

        /* Augmented-state forward integration */
//        std::cout<<"Start integration"<<std::endl;
        erkInteg.integrateStateFixed(); // Forward integration (augmented state): get final shooting state and cost
//        std::cout<<"...integration done."<<std::endl;

        /* Add stage-cost */
        vAL += *(xf+dx);

/*        std::cout<<"System input n= "<<n<<std::endl;
        for (i=0;i<du;++i)
            std::cout<<un[i]<<std::endl;
        std::cout<<"Integrator input-output n= "<<n<<std::endl;
        for (i=0;i<dw;++i)
            std::cout<<augStateIni[i]<<", "<<xf[i]<<std::endl; */

        /* Compute shooting equality constraint, final adjoint and partially separable AL gradient wrt s_{n+1} */
        for (i=0;i<dx;++i){
#if UNSCAL_VAR
            e = *(xff+i)-*(xf+i);
#endif
#if SCAL_VAR
            e = *(xff+i)-*(cs+i)-*(Ds+i)**(xf+i);
#endif
            r2e = halfRho*e;
            mr2e = *(ms+i)+r2e;
            vAL += mr2e*e;
            mr2e += r2e;
            /* Store in shooting constraint */
            *es++ = e;
            /* Initialize final augmented-adjoint (state part) */
#if UNSCAL_VAR
            *(adfinS+i) = -mr2e;
#endif
#if SCAL_VAR
            *(adfinS+i) = -*(Ds+i)*mr2e;
#endif
            /* AL gradient wrt s_{i+1} */
            *(gALs_+i) = mr2e;
        }
        /* Set final augmented-adjoint */
        *(adfinS+dx) = 1.; // Final cost-adjoint
        memset(adfinS+dw,0,dubl); // Final input-adjoint
        
        /* Backwards state interpolation over shooting interval: get intermediate states for backwards integration */
        erkInteg.interpState();
        
        /* Augmented-adjoint backward integration */
        erkInteg.setAdjointFin(adfinS); // Set final adjoint
        erkInteg.setAdjointIni(adjs); // Set initial adjoint (output of backwards integration)
        erkInteg.integrateAdjointFixed(); // Backwards-in-time integration (augmented adjoint): get initial shooting adjoint
        
        /* AL gradient wrt s_n */
#if UNSCAL_VAR
        memcpy(gALs,adjs,dxbl);
#endif
#if SCAL_VAR
        for (i=0;i<dx;++i)
            *(gALs+i) = *(iDs+i)**(adjs+i);
#endif
        /* AL gradient wrt q_n */
        memcpy(gALs+dx,adjs+dw,dubl);

#if DPC
        /* Evaluate path-constraints */
        pcon.eval(econPath,xn,un,rp);
        
        for (i=0;i<dpc;++i){
            e = *(econPath+i);
            r2e = halfRho*e;
            mr2e = *(mp+i)+r2e;
            vAL += mr2e*e;
            mr2e += r2e;
            /* Store in path-constraints */
            *ep++ = e;
            /* Store (mu+rho*e) */
            *(gALp+i) = mr2e;
        }
        /* Evaluate path-constraints jacobians transpose */
        pcon.evalJacXt(xn,un); // wrt state
        pcon.evalJacUt(xn,un); // wrt input

        /* AL gradient wrt s_n, q_n */
        pcon.prdJacXt(gALs,gALp); // Product against state-jacobian transpose
        pcon.prdJacUt(gALs+dx,gALp); // Product against input jacobian transpose
        
        /* Update path pointers */
        rp += dpc;
        mp += dpc;
        gALp += dpc;
#endif
        /* Build full gradient */
        for (i=0;i<dxux;++i)
            *(g+i) += *gALs++;

        /* Update shooting and reference pointers */
        xn += dxu;
        un += dxu;
        xf += dw;
        xff += dxu;
        egaus += dwnsteps;
        erpaus += sdwnsteps;
        yR += dy;
        uR += du;
        ms += dx;
        adfinS += da;
        g += dxu;
        gALs_ += dxux;
    }
    /**
    * Constraint on initial states 4 & 5
    */
    g = gAL+3;
    gALs = gALshoot+3;
    sini = z+3;
    while (ein != econ+2){
        e = *sini++;
        r2e = halfRho*e;
        mr2e = (*min++)+r2e;
        vAL += mr2e*e;
        mr2e += r2e; 
        /* Store initial constraint */
        *ein++ = e;
        /* Gradient wrt s_0 */
        *gALs++ += mr2e;
        *g++ += mr2e;
    }
    /**
    * Periodicity constraint on first 3 states
    */
    /* Pointer to full AL gradient */
    g = gAL;
    /* Pointer to begin shooting node in separated AL gradient */
    gALs = gALshoot;
    /* Pointer to terminal-state part of separated AL gradient (shooting) */ 
    gALsT = gALs+igTerm; 
    /* Pointer to terminal-state part of full AL gradient */ 
    gALt = gAL+izTerm;
    /* Pointer to intial and final shooting nodes */
    sini = z;
    sfin = sini+izTerm;
    while (epr != econ+dL){
        e = *sini-*sfin;
        ++sini;
        ++sfin;
        r2e = halfRho*e;
        mr2e = (*mpr++)+r2e;
        vAL += mr2e*e;
        mr2e += r2e;
        /* Store periodicity constraint */
        *epr++ = e;
        /* AL gradient wrt s_0 */
        *g++ += mr2e;
        *gALs++ += mr2e;
        /* AL gradient wrt s_N */
        *gALsT++ -= mr2e;
        *gALt++ -= mr2e; 
    }
#if MAY    
    /**
    * Mayer term
    */
    /* Mayer objective */
    may = mayer.eval(z+izTerm);
    vAL += may; 
    /* Mayer gradient */
    mayer.grad(z+izTerm);   
    gMay = mayer.getGrad(); 
    for (i=0;i<dx;++i){
#if UNSCAL_VAR
        gm = *(gMay+i);
#endif
#if SCAL_VAR
        gm = *(iDs+i)**(gMay+i);
#endif
        *(gALsT+i) += gm;
        *(gALt+i) += gm;
    }
#endif
    return vAL;
}

/**
* Multiple-shooting loop to get augmented states, augmented Lagrangian and part of the gradient of the augmented Lagrangian
* z: primal optimizer (s_0',q_0',s_1',q_1',...,s_{N-1}',q_{N-1}',s_N',r_0',...r_{N-1}')'
* s_n: state node
* q_n: input node
* r_n: slack associated with path-constraint at node n 
* mu: dual optimizer (initial-value embedding + shooting constraints + path-constraints)
* halfRho: halved penalty parameter (AL formulation)
* vAL: AL value at iterate
* gALshoot, dffGpath: copies of AL gradient parts (fill only the s_{n+1} parts (available from shooting constraints evaluation) and s_0)
*/
template <class dynT,class costT,class pconT,class mayerT> double MuShoot<dynT,costT,pconT,mayerT>::shootStatesFromIterate(double halfRho){

    int i,n;
    double *xn,*un;
    double *xf,*xff;
    double *es,*ee;
    double *ms;
    double *yR,*uR;
    double *adfinS;
    double *gALs;
    double *egaus,*erpaus;
    double e,r2e,mr2e;
    double may,vAL;
#if DPC    
    double *rp,*ep,*mp;
    double *gALp;
#endif   
//    std::cout<<"shoot sta from iterate: "<<std::endl;

    /** 
    * Initialize pointers
    */
    /* Initial shooting state s_n */
    xn = z;
    /* Shooting input q_n */
    un = xn+dx;
    /* Final shooting state w(t_{n+1};s_n,q_n) (integrator output) */
    xf = augStates;
    /* Final shooting state s_{n+1} (node, shooting constraint xff-xf=0) */ 
    xff = z+dxu;
    /* Intermediate augmented states */
    egaus = integAugStates;
    /* Interpolated augmented states */
    erpaus = interpAugStates;

    /* Output reference */
    yR = yRef; 
    /* Input reference */
    uR = uRef; 
    
    /* Pointer to shooting duals, skip initial-value embedding */
    ms = mu+dx;  
    
    /* Pointer to initial-value embedding constraint */
    ee = econ;
    /* Pointer to shooting constraints */
    es = ee+dx;
    
    /* Current final augmented adjoint (shooting part) */
    adfinS = adfin; 
    /* Current AL gradient copy (shooting part) */
    gALs = gALshoot+dxu;

#if DPC
    /* Path-constraints slack r_n */
    rp = z+dZshoot; 
    /* Pointer to path duals */
    mp = ms+dLshoot;
    /* Pointer to slacked path-constraints */
    ep = es+dLshoot;
    /* Path-constraint part of AL gradient */
    gALp = gALpath;
#endif

    /** 
    * Iterate over shooting intervals, s_{n+1}-x(t_{n+1};s_n,q_n) = 0 with augmented-state integration on [t_n,t_{n+1}]
    */
    vAL = 0.;
    for (n=0;n<Ns_;++n){
        /* Set parameters for integrator (input,reference,initial state) */
        erkInteg.setOutputRef(yR); // Set pointer to output-reference
        erkInteg.setInputRef(uR); // Set pointer to input-reference
        erkInteg.setInput(un); // Set pointer to input
        memcpy(augStateIni,xn,dxbl); // Copy initial state into initial augmented-state
        erkInteg.setStateFin(xf); // Set pointer to final augmented-state
        erkInteg.setStates(egaus); // Set pointer to intermediate augmented states 
        erkInteg.setStatesInterp(erpaus); // Set pointer to interpolated augmented states

        /* Augmented-state forward integration */
        erkInteg.integrateStateFixed(); // Forward integration (augmented state): get final shooting state and cost

        /* Add stage-cost */
        vAL += *(xf+dx);

        /* Compute shooting equality constraint, final adjoint and partially separable AL gradient wrt s_{n+1} */
        for (i=0;i<dx;++i){
            e = *(xff+i)-*(xf+i);
            r2e = halfRho*e;
            mr2e = *(ms+i)+r2e;
            vAL += mr2e*e;
            mr2e += r2e;
            /* Store in shooting constraint */
            *es++ = e;
            /* Initialize final augmented-adjoint (state part) */
            *(adfinS+i) = -mr2e;
            /* AL gradient wrt s_{i+1} */
            *(gALs+i) = mr2e;
        }

#if DPC
        /* Evaluate path-constraints */
        pcon.eval(econPath,xn,un,rp);
        
        for (i=0;i<dpc;++i){
            e = *(econPath+i);
            r2e = halfRho*e;
            mr2e = *(mp+i)+r2e;
            vAL += mr2e*e;
            mr2e += r2e;
            /* Store in slacked path-constraint */
            *ep++ = e;
            /* Store (mu+rho*e) */
            *(gALp+i) = mr2e;
        }
        /* Update path pointers */
        rp += dpc;
        mp += dpc;
        gALp += dpc;
#endif
        /* Update shooting & reference pointers */
        xn += dxu;
        un += dxu;
        xf += dw;
        xff += dxu;
        egaus += dwnsteps;
        erpaus += sdwnsteps;
        yR += dy;
        uR += du;
        ms += dx;
        adfinS += da;
        gALs += dxux;
    }

    /**
    * Initial value embedding, s_0-xest = 0
    */
    for (i=0;i<dx;++i){  
        e = *(z+i)-*(xest+i);
        r2e = halfRho*e;
        mr2e = *(mu+i)+r2e;
        vAL += mr2e*e;
        mr2e += r2e;
        /* Store in initial-value embedding constraint */
        *ee++ = e;
        /* AL gradient wrt s_0 */
        *(gALshoot+i) = mr2e;
    }
#if MAY
    /**
    * Mayer term 
    */
    may = mayer.eval(z+izTerm);
    vAL += may; 
#endif

    return vAL;
}

/**
* Multiple-shooting loop to get augmented adjoints and rest of gradient of augmented Lagrangian
* z: primal optimizer (s_0',q_0',s_1',q_1',...,s_{N-1}',q_{N-1}',s_N',r_0',...r_{N-1}')'
* s_n: state node
* q_n: input node
* r_n: slack associated with path-constraint at node n 
* mu: dual optimizer (initial-value embedding + shooting constraints + path-constraints)
* gALshoot
*/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::shootAdjointsFromIterate(){

    int i,n;
    double gm;
    double *xn,*un,*xf;
    double *yR,*uR;
    double *gALs,*gALp,*gALt,*gALsT;
    double *adjs,*adfinS;
    double *gMay;
    double *g;
    double *egaus,*erpaus;
    
//    std::cout<<"Shoot adj from iterate: "<<std::endl;

    /**
    * Initialize pointers 
    */
    /* Pointer to initial shooting augmented adjoint  */
    adjs = augAdjoints;
    /* Pointer to final shooting augmented adjoint */ 
    adfinS = adfin; 
    /* Separated gradient (s_n,q_n,s_{n+1}) */
    gALs = gALshoot;
    /* Pointer to gradient of path-constraint part of AL */
    gALp = gALpath;
    /* Pointer to terminal-state part of full AL gradient */
    gALt = gAL+izTerm;
    /* Pointer to terminal-state part of separated AL gradient */
    gALsT = gALs+igTerm;
    /* Pointer to full gradient */
    g = gAL;
    /* Pointer to output reference */
    yR = yRef; 
    /* Pointer to input reference */
    uR = uRef; 
    /* Pointer to initial shooting state */
    xn = z; 
    /* Pointer to shooting input */ 
    un = xn+dx;
    /* Pointer to final shooting augmented-state */
    xf = augStates; 
    /* Pointer to intermediate augmented states */
    egaus = integAugStates;
    /* Pointer to interpolated augmented states */
    erpaus = interpAugStates;
    
    /* Set shooting part of full AL gradient to 0 */
    memset(gAL,0,dGshootDBL);
    /* Copy initial value embedding AL term into full gradient gAL */
    memcpy(gAL,gALshoot,dxbl);

    /**
    * Loop over shooting intervals with augmented-adjoint integration on [t_n,t_{n+1}]
    */
    for (n=0;n<Ns_;++n){    
        /* Set integration parameters */
        erkInteg.setOutputRef(yR); // Set pointer to output-reference
        erkInteg.setInputRef(uR); // Set pointer to input-reference
        erkInteg.setInput(un); // Set pointer to input (parameter)
        memcpy(augStateIni,xn,dxbl); // Copy initial state into initial augmented-state
        erkInteg.setStateFin(xf); // Set pointer to final augmented-state           
        erkInteg.setStates(egaus); // Set pointer to intermediate augmented states 
        erkInteg.setStatesInterp(erpaus); // Set pointer to interpolated augmented states
        /* Set final augmented-adjoint */
        *(adfinS+dx) = 1.; // Final cost-adjoint
        memset(adfinS+dw,0,dubl); // Final input-adjoint
        /* Backwards state interpolation over shooting interval: get intermediate states for backwards integration */
        erkInteg.interpState();
        /* Augmented-adjoint backward integration */
        erkInteg.setAdjointFin(adfinS); // Set final adjoint
        erkInteg.setAdjointIni(adjs); // Set initial adjoint (output of backwards integration)
        erkInteg.integrateAdjointFixed(); // Backwards-in-time integration (augmented adjoint): get initial shooting adjoint
        /* AL gradient wrt s_n */
        memcpy(gALs,adjs,dxbl);
        /* AL gradient wrt q_n */
        memcpy(gALs+dx,adjs+dw,dubl);
#if DPC           
        /* Evaluate path-constraints jacobians transpose */
        pcon.evalJacXt(xn,un); // wrt state
        pcon.evalJacUt(xn,un); // wrt input
        /* AL gradient wrt s_n, q_n */
        pcon.prdJacXt(gALs,gALp); // Product against state-jacobian transpose
        pcon.prdJacUt(gALs+dx,gALp); // Product against input jacobian transpose
        /* Update path gradient pointer */
        gALp += dpc;
#endif
        /* Build full gradient */
        for (i=0;i<dxux;++i)
            *(g+i) += *gALs++;
        /* Update shooting and reference pointers */
        xn += dxu;
        un += dxu;
        xf += dw;
        egaus += dwnsteps;
        erpaus += sdwnsteps;
        yR += dy;
        uR += du;
//        adjs += da;
        adfinS += da;
        g += dxu;
    }  
    /* Initial-value embedding in separated AL gradient */
    memcpy(gALshoot,gAL,dxbl);
#if MAY
    /* Mayer term */
    mayer.grad(z+izTerm);   
    gMay = mayer.getGrad(); 
    for (i=0;i<dx;++i){
        gm = *(gMay+i);
        *(gALt+i) += gm;
        *(gALsT+i) += gm;
    }
#endif
}

/************************
******* For debugging ***
************************/
/**
** Display state transformation
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::displayStateTransfo(){
    
    int i;

    std::cout<<"State centering cs:"<<std::endl;
    for (i=0;i<dx;++i)
        std::cout<<"cs["<<i<<"]="<<cs[i]<<std::endl;

    std::cout<<"State scaling Ds:"<<std::endl;
    for (i=0;i<dx;++i)
        std::cout<<"Ds["<<i<<"]="<<Ds[i]<<std::endl;

    std::cout<<"Inverse state scaling iDs:"<<std::endl;
    for (i=0;i<dx;++i)
        std::cout<<"iDs["<<i<<"]="<<iDs[i]<<std::endl;
}

/**
** Display input transformation 
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::displayInputTransfo(){
    
    int i;

    std::cout<<"Input centering cq:"<<std::endl;
    for (i=0;i<du;++i)
        std::cout<<"cq["<<i<<"]="<<cq[i]<<std::endl;

    std::cout<<"Input scaling Dq:"<<std::endl;
    for (i=0;i<du;++i)
        std::cout<<"Dq["<<i<<"]="<<Dq[i]<<std::endl;
    
    std::cout<<"Inverse state scaling iDq:"<<std::endl;
    for (i=0;i<du;++i)
        std::cout<<"iDq["<<i<<"]="<<iDq[i]<<std::endl;

}

/* Display multiple-shooting parameters */
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::displayShootingParameters(){

    std::cout<<"===========> Shooting parameters <=============="<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"State dimension dx: "<<dx<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Input dimension du: "<<du<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Adjoint dimension da: "<<da<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Augmented-state dimension dw: "<<dw<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Output dimension dy: "<<dy<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Path-constraint dimension dpc: "<<dpc<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Prediction time Tpred: "<<Tpred<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Number of shooting nodes Ns: "<<Ns<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Number of shooting intervals Ns_: "<<Ns_<<std::endl;
    std::cout<<"============================================================="<<std::endl;
    std::cout<<""<<std::endl;

    erkInteg.displayParameters();
    erkInteg.displayButcher();
}

/**
** Display infinity norm of augmented state 
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::displayAsInfNorm(){
    
    double nIas;
    double *as;

    nIas = -1.;
    as = augStates;
    while (as != augStates+dAugS){
        nIas = MAX(nIas,ABS(*as));
        ++as;
    }   
    std::cout<<"<<<MuShoot>>> Infinity norm of augmented states: "<<nIas<<std::endl;
}

/**
** Display infinity norm of augmented adjoints 
**/
template <class dynT,class costT,class pconT,class mayerT> void MuShoot<dynT,costT,pconT,mayerT>::displayAaInfNorm(){
    
    double nIaa;
    double *aa;
    
    nIaa = -1.;
    aa = augAdjoints;
    while (aa != augAdjoints+da){
        nIaa = MAX(nIaa,ABS(*aa));
        ++aa;
    }
    std::cout<<"<<<MuShoot>>> Infinity norm of augmented adjoints: "<<nIaa<<std::endl;
}


