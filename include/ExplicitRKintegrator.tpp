//
//  ExplicitRKintegrator.cpp
//  MusAL
//
//  Created by Jean on 3/12/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//
#include <iostream>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "macros.h"
#include "ExplicitRKintegrator.h"

/** 
* Default constructor 
*/
template <class dynT> ExplicitRKintegrator<dynT>::ExplicitRKintegrator(){

    dx = ode.getDx();
    dxbl = dx*SZDBL;
    if (dx<=0){
        std::cerr<<"ERROR<ExplicitRKintegrator>: state-dimension nonpositive."<<std::endl;
        return;
    }

    du = ode.getDu();
	dubl = du*SZDBL;
    if (du<=0){
        std::cerr<<"ERROR<ExplicitRKintegrator>: input-dimension nonpositive."<<std::endl;
        return;
    }

    da = ode.getDa();
	dabl = da*SZDBL;

    A = NULL;
    b = NULL;
    c = NULL;

	states = NULL;
	adjoints = NULL;
	statesInterp = NULL;
    stateNodes = NULL;
    adjointNodes = NULL; 

    /* Allocate auxiliary state & adjoint */
    auxState = new double[dx];
    auxAdjoint = new double[da];

    hsta = NULL;
    hadj = NULL;
    xini = NULL;
    xfin = NULL;
    adfin = NULL;
    adini = NULL;
    u = NULL;
    
    /* Set default-value for max. # integration steps */
    nstepsMax = NSMAX;
    
    /* Alocate step-size history for state & adjoint integration */
    hsta = new double[nstepsMax];
    hadj = new double[nstepsMax];
    
    /* Set integration interval and step-sizes to default values */
    timeInterval = -1.;
    hfix = -1.;
    hini = -1.;
    hmin = -1.;
    hmax = -1.;
}

/**
* Constructor from data
*/
/*template <class dynT> ExplicitRKintegrator<dynT>::ExplicitRKintegrator(int s,double* A,double* b,double* c,int dx,int du,int nsMax,double hmin,double hmax){
    
    assert(s>0);
    assert(dx>0);
    assert(du>0);
    assert(nsMax>0);
    assert(hmin>0);
    assert(hmax>0);
    
    int ss = s*s;
    
    this->s = s;
    
    this->A = new double[ss];
    std::copy(A,A+ss,this->A);
    
    this->b = new double[s];
    std::copy(b,b+s,this->b);
    
    this->c = new double[s];
    std::copy(c,c+s,this->c);
    
    this->dx = dx;
    this->du = du;
    this->da = dx+du;
    
    this->nstepsMax = nsMax;
    this->hmin = hmin;
    this->hmax = hmax;
    
    this->states = new double[dx*(nsMax+1)];
	this->adjoints = new double[da*(nsMax+1)];
	this->statesInterp = new double[s*dx*nsMax];
    
    this->hsta = new double[nsMax];
    this->hadj = new double[nsMax];
    
    this->stateNodes = new double[s*dx];
    this->adjointNodes = new double[s*da];
    
    this->auxState = new double[dx];
    this->auxAdjoint = new double[da];
    
    xini = NULL;
    xfin = NULL;
    adfin = NULL;
    adini = NULL;
    u = NULL;
}
*/

/**
* Copy constructor
*/
template <class dynT> ExplicitRKintegrator<dynT>::ExplicitRKintegrator(const ExplicitRKintegrator<dynT>& obj){
    
    int ss;
    
    s = obj.s;
    ss = s*s;
    
    A = new double[ss];
    std::copy(obj.A,obj.A+ss,A);
    
    b = new double[s];
    std::copy(obj.b,obj.b+obj.s,b);
    
    c = new double[s];
    std::copy(obj.c,obj.c+obj.s,c);
    
    stateNodes = new double[s*dx];
    std::copy(obj.stateNodes,obj.stateNodes+obj.s*obj.dx,stateNodes);
    
    adjointNodes = new double[s*da];
    std::copy(obj.adjointNodes,obj.adjointNodes+obj.s*obj.da,adjointNodes);
    
    auxState = new double[dx];
    std::copy(obj.auxState,obj.auxState+obj.dx,auxState);
    
    auxAdjoint = new double[da];
    std::copy(obj.auxAdjoint,obj.auxAdjoint+obj.da,auxAdjoint);
  
    statesInterp = NULL;  
//    statesInterp = new double[obj.s*obj.dx*obj.nsteps];
//    std::copy(obj.statesInterp,obj.statesInterp+obj.s*obj.dx*obj.nsteps,statesInterp);
    
 	states = NULL;
//  states = new double[obj.dx*(obj.nsteps+1)];
// 	std::copy(obj.states,obj.states+obj.dx*(obj.nsteps+1),states);
 	
 	adjoints = new double[obj.da*(obj.nsteps+1)];
 	std::copy(obj.adjoints,obj.adjoints+obj.da*(obj.nsteps+1),adjoints);

    ode = obj.ode;
    
    timeInterval = obj.timeInterval;

    u = obj.u;
    adini = obj.adini;
    adfin = obj.adfin;
    xfin = obj.xfin;
    xini = obj.xini;
    
    nstepsMax = obj.nstepsMax;
    nsteps = obj.nsteps;
    sdx = obj.sdx;
    sdxnsteps = obj.sdxnsteps;
    dxnsteps = obj.dxnsteps;
    
    hadj = new double[nstepsMax];
    std::copy(obj.hadj,obj.hadj+obj.nstepsMax,hadj);
    
    hsta = new double[nstepsMax];
    std::copy(obj.hsta,obj.hsta+obj.nstepsMax,hsta);
    
    hfix = obj.hfix;
    hmin = obj.hmin;
    hmax = obj.hmax;
    hini = obj.hini;
    da = obj.da;
    dabl = obj.dabl;
    du = obj.du;
    dubl = obj.dubl;
    dx = obj.dx;
    dxbl = obj.dxbl;
}

/**
** Destructor
**/
template <class dynT> ExplicitRKintegrator<dynT>::~ExplicitRKintegrator(){
    
//    std::cout<<"~ExplicitRKintegrator()"<<std::endl;
    if (A != NULL)
        delete [] A;

    if (b != NULL)
        delete [] b;

    if (c != NULL)
        delete [] c;
    
//    if (states != NULL)
//    	delete [] states;
    	
  	if (adjoints != NULL)	
    	delete [] adjoints;
    
//    if (statesInterp != NULL)
//    	delete [] statesInterp;
    
    if (stateNodes != NULL)
        delete [] stateNodes;
    
    if (adjointNodes != NULL)
        delete [] adjointNodes;
    
    if (auxState != NULL)
        delete [] auxState;
    
    if (auxAdjoint != NULL)
        delete [] auxAdjoint;
    
    if (hsta != NULL)
        delete [] hsta;
    
    if (hadj != NULL)
        delete [] hadj;
//    std::cout<<"~ExplicitRKintegrator()"<<std::endl;
}

/**
** Assignment operator
**/
template <class dynT> ExplicitRKintegrator<dynT>& ExplicitRKintegrator<dynT>::operator= (const ExplicitRKintegrator<dynT>& obj){
    
    int ss;
    double *cA,*cb,*cc;
    double *stateNodes_,*adjointNodes_;
    double *auxState_,*auxAdjoint_;
    double *c_hsta,*c_hadj;
//	double *states_,*statesInterp_;
    double *adjoints_;

    if (this != &obj){

        dx = obj.dx;
        dxbl = obj.dxbl;
        du = obj.du;
        dubl = obj.dubl;
        da = obj.da;
        dabl = obj.dabl;

        hini = obj.hini;
        hmax = obj.hmax;
        hmin = obj.hmin;
        hfix = obj.hfix;
        
        sdx = obj.sdx;
        nsteps = obj.nsteps;
        nstepsMax = obj.nstepsMax;
        sdxnsteps = obj.sdxnsteps;
        dxnsteps = obj.dxnsteps;

        xini = obj.xini;
        xfin = obj.xfin;
        adfin = obj.adfin;
        adini = obj.adini;
        u = obj.u;

        states = obj.states;
        statesInterp = obj.statesInterp;

        timeInterval = obj.timeInterval;
        
        ode = obj.ode;
        
        s = obj.s;
        ss = s*s;

        /* hsta */
        if (hsta != NULL)
            delete [] hsta;
        if (obj.hsta != NULL){
            c_hsta = new double[nstepsMax];
            std::copy(obj.hsta,obj.hsta+nstepsMax,c_hsta);
            hsta = c_hsta;
        }
        else {
            hsta = NULL;
        }

        /* hadj */
        if (hadj != NULL)
            delete [] hadj;
        if (obj.hadj != NULL){
            c_hadj = new double[nstepsMax];
            std::copy(obj.hadj,obj.hadj+nstepsMax,c_hadj);
            hadj = c_hadj;
        }
        else {
            hadj = NULL;
        }
        
        /* A */
        if (A != NULL)
            delete [] A;
        if (obj.A != NULL){
            cA = new double[ss];
            std::copy(obj.A,obj.A+ss,cA);
            A = cA;
        }
        else {
            A = NULL;
        }
        
        /* b */
        if (b != NULL)
            delete [] b;
        if (obj.b != NULL){
            cb = new double[s];
            std::copy(obj.b,obj.b+s,cb);
            b = cb;
        }
        else {
            b = NULL;
        }   
        
        /* c */
        if (c != NULL)
            delete [] c;
        if (obj.c != NULL){
            cc = new double[s];
            std::copy(obj.c,obj.c+s,cc);
            c = cc;
        }
        else {
            c = NULL;
        }

        /* state Nodes */
        if (stateNodes != NULL)
            delete [] stateNodes;
        if (obj.stateNodes != NULL){
            stateNodes_ = new double[s*dx];
            std::copy(obj.stateNodes,obj.stateNodes+s*dx,stateNodes_);
            stateNodes = stateNodes_;
        }
        else {
            stateNodes = NULL;
        }

        /* adjointNodes  */
        if (adjointNodes != NULL)
            delete [] adjointNodes;
        if (obj.adjointNodes != NULL){
            adjointNodes_ = new double[s*da];
            std::copy(obj.adjointNodes,obj.adjointNodes+s*da,adjointNodes_);
            adjointNodes = adjointNodes_;
        } 
        else {
            adjointNodes = NULL;
        }

        /* auxState */
        if (auxState != NULL)
            delete [] auxState;
        if (obj.auxState != NULL){
            auxState_ = new double[dx];
            std::copy(obj.auxState,obj.auxState+dx,auxState_);
            auxState = auxState_;
        } 
        else {
            auxState = NULL;
        }

        /* auxAdjoint */
        if (auxAdjoint != NULL)
            delete [] auxAdjoint;
        if (obj.auxAdjoint != NULL){
            auxAdjoint_ = new double[da];
            std::copy(obj.auxAdjoint,obj.auxAdjoint+da,auxAdjoint_);
            auxAdjoint = auxAdjoint_;
        } 
        else {
            auxAdjoint = NULL;
        }

        /* states */
    /*    if (states != NULL)
        	delete [] states;
        
        if (obj.states != NULL){
            states_ = new double[(nsteps+1)*dx];
            std::copy(obj.states,obj.states+dx*(nsteps+1),states_);
            states = states_;
        } 
        else {
            states = NULL;
        } */

        /* adjoints */
        if (adjoints != NULL)
        	delete [] adjoints;
        if (obj.adjoints != NULL){
            adjoints_ = new double[(nsteps+1)*da];
            std::copy(obj.adjoints,obj.adjoints+da*(nsteps+1),adjoints_);
            adjoints = adjoints_;
        } 
        else {
            adjoints = NULL;
        }

        /* statesInterp */
    /*    if (statesInterp != NULL)
        	delete [] statesInterp;
        	
        if (obj.statesInterp != NULL){
       	    statesInterp_ = new double[nsteps*s*dx];
            std::copy(obj.statesInterp,obj.statesInterp+nsteps*s*dx,statesInterp_);
            statesInterp = statesInterp_;
        } 
        else {
            statesInterp = NULL;
        } */
    }
    return *this;
}

/**
** Set Butcher tableau and modify data accordingly 
**/
template <class dynT> void ExplicitRKintegrator<dynT>::setButcher(int s,double* Anew,double* bnew,double* cnew){

	assert(nstepsMax>0);
	assert(dx>0);

    int ss;

    ss = s*s;
    this->s = s;
    
    this->sdx = s*dx;

    /* Reallocate A */
    if (A != NULL)
        delete [] A;    
    double* cA = new double[ss];
    std::copy(Anew,Anew+ss,cA);
    A = cA;
    
    /* Reallocate b */
    if (b != NULL)
        delete [] b;
    double* cb = new double[s];
    std::copy(bnew,bnew+s,cb);
    b = cb;
    
    /* Reallocate c */
    if (c != NULL)
        delete [] c;
    double* cc = new double[s];
    std::copy(cnew,cnew+s,cc);
    c = cc;
    
    /* Reallocate stateNodes */
    if (stateNodes != NULL)
        delete [] stateNodes;
    double* cstateNodes = new double[s*dx];
    stateNodes = cstateNodes; 

    /* Reallocate adjointNodes */
    if (adjointNodes != NULL)
        delete [] adjointNodes;
    double* cadjointNodes = new double[s*da];
    adjointNodes = cadjointNodes; 
}

/** 
** Set integration parameters 
**/
template <class dynT> void ExplicitRKintegrator<dynT>::setIntegSteps(int nsteps,int nstepsMax,double hmin,double hmax,double hini,double hfix){

	assert(sdx>0);
	assert(nsteps<nstepsMax);
	assert(da>0);
	assert(dx>0);
	assert(hmin>0);
	assert(hmax>0);
	assert(hfix>0);
	assert(hini>0);

    int nsMAXmin;

    nsMAXmin = MIN(this->nstepsMax,nstepsMax);

    this->nsteps = nsteps;
    this->nstepsMax = nstepsMax;
    this->sdxnsteps = sdx*nsteps;
    
    this->hmin = hmin;
    this->hmax = hmax;
    this->hini = hini;
    this->hfix = hfix;
    
    /* Reallocate hsta */
    double* chsta = new double[nstepsMax];
    if (hsta != NULL){
        std::copy(hsta,hsta+nsMAXmin,chsta);
        delete [] hsta;
    }
    hsta = chsta;

    /* Reallocate hadj */
    double* chadj = new double[nstepsMax];
    if (hadj != NULL){
        std::copy(hadj,hadj+nsMAXmin,chadj);
        delete [] hadj;
    }
    hadj = chadj;
    
    /* Reallocate states */
/*    if (states != NULL)
    	delete [] states;
    
    double* states_ = new double[dx*(nstepsMax+1)];
	states = states_; */

	/* Reallocate adjoints */
	if (adjoints != NULL)
		delete [] adjoints;
	double* adjoints_ = new double[da*(nstepsMax+1)];
	adjoints = adjoints_;
	
	/* Reallocate statesInterp */
/*	if (statesInterp != NULL)
		delete [] statesInterp;
	
	double* statesInterp_ = new double[nstepsMax*sdx];
	statesInterp = statesInterp_; */	
}

/**
** Initialization method (called once Butcher data & step-size set)
**/
template <class dynT> void ExplicitRKintegrator<dynT>::init(){

	assert(dx>=1);
	assert(s>=1);
	assert(sdx>=1);
	assert(da>=1);
	assert(timeInterval>0.);
    assert(nsteps>=1);

    if (nsteps>nstepsMax){
        std::cerr<<"ERROR<ExplicitRKintegrator>: computed number of steps larger than maximal allowed number."<<std::endl;
        return;
    }

	hfix = timeInterval/nsteps;

    if (hfix<=0.){
        std::cerr<<"ERROR<ExplicitRKintegrator>: integration step-size nonpositive."<<std::endl;
        return;
    }

	sdxnsteps = sdx*nsteps;
    dxnsteps = dx*(nsteps+1);
	
    /* states */
/*    if (states != NULL)
        delete [] states;
    states = new double[dxnsteps]; */
	
    /* adjoints */
    if (adjoints != NULL)
        delete [] adjoints;
    adjoints = new double[da*(nsteps+1)];

    /* statesInterp */
/*    if (statesInterp != NULL)
        delete [] statesInterp;
	statesInterp = new double[sdxnsteps]; */
}   

/** 
** Display Butcher-tableau 
**/
template <class dynT> void ExplicitRKintegrator<dynT>::displayButcher(){

    int i,ss,s_;
    
    ss = s*s;
    s_ = s-1;

    std::cout<<"===========> BUTCHER TABLEAU <==========="<<std::endl;
    std::cout<<"Integration order: "<<s<<std::endl;
    std::cout<<" "<<std::endl;
    std::cout<<"A matrix:"<<std::endl;
    for (i=0;i<ss;++i){
        if (i%s == s_){
            std::cout<<A[i]<<std::endl;
        }
        else {
            std::cout<<A[i]<<" ";
        }
    }
    std::cout<<" "<<std::endl;
    std::cout<<"b vector:"<<std::endl;
    for (i=0;i<s;++i)
        std::cout<<b[i]<<std::endl;
    std::cout<<" "<<std::endl;
    std::cout<<"c vector:"<<std::endl;
    for (i=0;i<s;++i)
        std::cout<<c[i]<<std::endl;
    std::cout<<" "<<std::endl;
    std::cout<<"========================================"<<std::endl;
    std::cout<<" "<<std::endl;

}

/**
** Display integration parameters
**/
template <class dynT> void ExplicitRKintegrator<dynT>::displayParameters(){

    std::cout<<"===========> ExplicitRKintegrator parameters <=============="<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"State dimension dx: "<<dx<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Input dimension du: "<<du<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Adjoint dimension da: "<<da<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Fixed integration step-size hfix: "<<hfix<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Number of integration steps: "<<nsteps<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Max. number of integration steps: "<<nstepsMax<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"Time interval: "<<timeInterval<<std::endl;
    std::cout<<"============================================================="<<std::endl;
    std::cout<<""<<std::endl;
}

/**  
** Backwards-in-time Linear state interpolation between xini and xfin based on Butcher-vector c 
**/
template <class dynT> void ExplicitRKintegrator<dynT>::interpState(){
    
    int i,j;
    double cc;
    double *xb,*xe,*xint;

    /* Loop over integration steps */   
//    for (n=0;n<nsteps;++n){
    	/* Linear interpolation on nth integration interval */
//    	for (i=0;i<s;++i){
//    		cc = *(c+i);	
//    		for (j=0;j<dx;++j)
//    			*xint++ = *(xb+j)*(1-cc)+*(xe+j)*cc;
//    	}  	
//    	xb += dx;
//    	xe += dx;
//    }

    /* Backwards loop over integration steps */
    xe = states+dxnsteps-dx;
    xb = xe-dx;

 //   std::cout<<"xe="<<std::endl;
 //   for (i=0;i<dx;++i)
 //       std::cout<<xe[i]<<std::endl;
 //   std::cout<<"xb="<<std::endl;
 //   for (i=0;i<dx;++i)
 //       std::cout<<xb[i]<<std::endl;
    xint = statesInterp+sdxnsteps-dx;
    while (1){
        /* Linear interpolation on nth integration interval */
        for (i=0;i<s;++i){
            cc = *(c+i);
            for (j=0;j<dx;++j)
                *(xint+j) = *(xe+j)*(1-cc)+*(xb+j)*cc;
            xint -= dx;
        }
        xe -= dx;
        xb -= dx;
        if (xe==states)
            break;
    }

 /*   xint = statesInterp;
    for (n=0;n<nsteps;++n){

        std::cout<<"integ step "<<n<<std::endl;

        for (i=0;i<s;++i){
            std::cout<<"stage "<<i<<std::endl;
            for (j=0;j<dx;++j)
                std::cout<<*xint++<<std::endl;
        }

    } */
   /* 	for (j=0;j<dx;++j){
        	xbj = *(xini+j);
        	xej = *(xfin+j);
        	xij = xinter+j;
        	ix = 0;
        	for (i=0;i<s;++i){
            	cc = *(c+i);
            	*(xij+ix) = (1-cc)*xbj+cc*xej;
            	ix += dx;
        	}
    	} */
    
}

/**
** State integration loop (fixed step-size)
**/
template <class dynT> void ExplicitRKintegrator<dynT>::integrateStateFixed(){
    
    assert(nsteps>=1);

    int i,j,k,ii,is,kx,n;
    double aa,hh;
    double *Ais,*noj,*sb,*se;

	/* states <- xini */
    memcpy(states,xini,dxbl);

	/* Forward-in-time integration loop */
	sb = states;
	se = states+dx;
    n = 0;
    while (1){
//        std::cout<<"ERK forward integration step "<<n+1<<std::endl;
        /* Loop over integration stages */
        ode.rhs(stateNodes,sb,u);
        ii = dx;
        is = s;
        for (i=1;i<s;++i){
            Ais = A+is;
            for (j=0;j<dx;++j){
                noj = stateNodes+j;
                aa = 0.;
                kx = 0;
                for (k=0;k<i;++k){
                    aa += *(Ais+k)**(noj+kx);
                    kx += dx;
                }
                *(auxState+j) = *(sb+j)+hfix*aa;
            }
            ode.rhs(stateNodes+ii,auxState,u);
            ii += dx;
            is += s;
        }

//        std::cout<<"State nodes:"<<std::endl;
//        for (i=0;i<sdx;++i)
//            std::cout<<stateNodes[i]<<std::endl;
    
        /* Update final state */
        for (j=0;j<dx;++j){
            noj = stateNodes+j;
            aa = 0.;
            ii = 0;
            for (i=0;i<s;++i){
                aa += *(b+i)**(noj+ii);
                ii += dx;
            }
            *se++ = hfix*aa+*sb++;
        }
		/* Go to next integration step */
        ++n;
        if (n>=nsteps)
            break;
    }
    /* xfin <- se */
    memcpy(xfin,se-dx,dxbl);

//    std::cout<<"states in sta="<<std::endl;
//    for (i=0;i<dxnsteps;++i)
//        std::cout<<states[i]<<std::endl; 
}

/**
** Adjoint integration loop (BACKWARDS in time with fixed step-size):
** - assumes interpolated states to be time-sorted and last interpolate to be the final state of the integration step
**/
template <class dynT> void ExplicitRKintegrator<dynT>::integrateAdjointFixed(){
    
    assert(nsteps>=1);

    int n,i,ii,is,j,jj,k,kx;
    double aa,hh;
    double *noj,*Ais,*ab,*ae,*si;

	/* adjoints <- adfin */
	memcpy(adjoints,adfin,dabl);

	/* Backward-in-time integration loop */
	ae = adjoints;
	ab = adjoints+da;
	/* Move backwards on interpolated states */
	si = statesInterp+sdxnsteps-dx;

//    displayInterpolates();

//    std::cout<<"states in adj="<<std::endl;
//    for (i=0;i<dxnsteps;++i)
//        std::cout<<states[i]<<std::endl;

	n = 0;
	while (1){
//        std::cout<<"ERK backward integration step "<<n+1<<std::endl;
		/* Loop over stages */   
        ode.adjRhs(adjointNodes,si,u,ae);
//        std::cout<<"====================="<<std::endl;
//        std::cout<<"adjointNode "<<1<<": "<<std::endl;
//        for (k=0;k<da;++k)
//            std::cout<<adjointNodes[k]<<std::endl;
//        std::cout<<"Interpolated state "<<1<<": "<<std::endl;
//        for (k=0;k<dx;++k)
//            std::cout<<*(si+k)<<std::endl;
//        std::cout<<"Final state:"<<std::endl;
//        for (k=0;k<dx;++k)
//            std::cout<<*(xfin+k)<<std::endl;
        si -= dx;
    	jj = da;
    	is = s;
    	for (i=1;i<s;++i){
        	Ais = A+is;
//            std::cout<<"Aux adjoint "<<i+1<<": "<<std::endl;
        	for (j=0;j<da;++j){
            	noj = adjointNodes+j;
            	aa = 0.;
            	kx = 0;
            	for (k=0;k<i;++k){
                	aa += *(Ais+k)**(noj+kx);
                	kx += da;
            	}
            	*(auxAdjoint+j) = *(ae+j)-hfix*aa;
//                std::cout<<*(auxAdjoint+j)<<std::endl;
        	}
//            std::cout<<"State interp "<<std::endl;
//            for (j=0;j<dx;++j)
//                std::cout<<si[j]<<std::endl;
        	ode.adjRhs(adjointNodes+jj,si,u,auxAdjoint);
//            std::cout<<"adjointNode "<<i+1<<": "<<std::endl;
//            for (k=0;k<da;++k)
//                std::cout<<*(adjointNodes+jj+k)<<std::endl;
//            std::cout<<"Interpolated states "<<i+1<<": "<<std::endl;
//            for (k=0;k<dx;++k)
//                std::cout<<*(si+k)<<std::endl;
            si -= dx;
            jj += da;
        	is += s;
    	}
        
//        std::cout<<"Initial state:"<<std::endl;
//        for (k=0;k<dx;++k)
//            std::cout<<xini[k]<<std::endl;

//        std::cout<<"Adjoint nodes:"<<std::endl;
//        for (i=0;i<s*da;++i)
//            std::cout<<adjointNodes[i]<<std::endl;
//        std::cout<<"=================="<<std::endl;

//        std::cout<<"Initial interm state:"<<std::endl;
//        for (j=0;j<dx;++j)
//            std::cout<<states[dx*(nsteps+1)-(i+1)*dx-j]<<std::endl;

    	/* Update initial adjoint */
    	for (j=0;j<da;++j){
        	noj = adjointNodes+j;
        	aa = 0.;
        	ii = 0;
        	for (i=0;i<s;++i){
            	aa += *(b+i)**(noj+ii);
            	ii += da;
        	}
        	*ab++ = (*ae++)-hfix*aa;
    	}

		++n;

        if (n>=nsteps)
            break;
	}
    /* adini <- ab */
    memcpy(adini,ab-da,dabl);
}

/**
** Display states 
**/
template <class dynT> void ExplicitRKintegrator<dynT>::displayStates(){

    int j;

    std::cout<<"Initial state:"<<std::endl;
    for (j=0;j<dx;++j)
        std::cout<<xini[j]<<std::endl;

    std::cout<<"States: "<<std::endl;
    for (j=0;j<dxnsteps;++j)
        std::cout<<states[j]<<std::endl;

    std::cout<<"Final state:"<<std::endl;
    for (j=0;j<dx;++j)
        std::cout<<xfin[j]<<std::endl;
}

/**
** Display interpolated states
**/
template <class dynT> void ExplicitRKintegrator<dynT>::displayInterpolates(){

    int i,j;
    double* psi,*psiEnd;

    std::cout<<"Initial state:"<<std::endl;
    for (j=0;j<dx;++j)
        std::cout<<xini[j]<<std::endl;

    std::cout<<"Interpolated states:"<<std::endl;
    psi = statesInterp;
    psiEnd = psi+sdxnsteps;
    while (psi != psiEnd){

        for (i=0;i<s;++i){
    
            std::cout<<"States interp "<<i+1<<std::endl;
            for (j=0;j<dx;++j)
                std::cout<<*psi++<<std::endl;
            
        }
    } 
//    for (j=0;j<sdxnsteps;++j)
//        std::cout<<statesInterp[j]<<std::endl;

    std::cout<<"Final state:"<<std::endl;
    for (j=0;j<dx;++j)
        std::cout<<xfin[j]<<std::endl; 
}

/**
** Perform one single forward integration step 
**/
template <class dynT> void ExplicitRKintegrator<dynT>::performStepState(double* xend,double* xini,double* h){
    
    //    std::cout<<"Perform step state"<<std::endl;
    
    int i,j,k,ii,is,kx;
    double aa,hh;
    double *Ais,*noj;
    
    hh = *h;
    
    /* Loop over stages */
    ode.rhs(stateNodes,xini,u);
    ii = dx;
    is = s;
    for (i=1;i<s;++i){
        Ais = A+is;
        for (j=0;j<dx;++j){
            noj = stateNodes+j;
            aa = 0.;
            kx = 0;
            for (k=0;k<i;++k){
                aa += *(Ais+k)**(noj+kx);
                kx += dx;
            }
            *(auxState+j) = *(xini+j)+hh*aa;
        }
        ode.rhs(stateNodes+ii,auxState,u);
        ii += dx;
        is += s;
    }
    
    /* Update final state */
    for (j=0;j<dx;++j){
        noj = stateNodes+j;
        aa = 0.;
        ii = 0;
        for (i=0;i<s;++i){
            aa += *(b+i)**(noj+ii);
            ii += dx;
        }
        *(xend+j) = *(xini+j)+hh*aa;
    }
    
}

/**
** Perform one single adjoint integration step 
**/
template <class dynT> void ExplicitRKintegrator<dynT>::performStepAdjoint(double* adjIni,double* adjFin,double* h){
    
    int i,ii,is,j,jj,k,kx;
    double aa,hh;
    double *noj,*Ais;
    
    hh = *h;
    
    /* Loop over stages */
    ode.adjRhs(adjointNodes,statesInterp,u,adjFin);
    jj = da;
    is = s;
    ii = dx;
    for (i=1;i<s;++i){
        Ais = A+is;
        for (j=0;j<da;++j){
            noj = adjointNodes+j;
            aa = 0.;
            kx = 0;
            for (k=0;k<i;++k){
                aa += *(Ais+k)**(noj+kx);
                kx += dx;
            }
            *(auxAdjoint+j) = *(xini+j)-hh*aa;
        }
        ode.adjRhs(adjointNodes+jj,statesInterp+ii,u,auxAdjoint);
        jj += da;
        ii += dx;
        is += s;
    }
    
    /* Update initial adjoint */
    for (j=0;j<da;++j){
        noj = adjointNodes+j;
        aa = 0.;
        ii = 0;
        for (i=0;i<s;++i){
            aa += *(b+i)**(noj+ii);
            ii += da;
        }
        *(adjIni+j) = *(adjFin+j)-hh*aa;
    }
    
}

/**
** State integration main loop (varying step-size)
**/
template <class dynT> void ExplicitRKintegrator<dynT>::integrateStateAdaptive(){
    

	std::cout<<"ERROR<ExplicitRKintegrator>: State integration with adaptive step-size not implemented."<<std::endl; 

 /*   
 	assert(nsteps>=0);
 
 
 	int t,tx;
    double *elems,*tims;
    
    elems = stateGrid.getElems();
    tims = stateGrid.getTime();
    
        
    *hsta = hini;
    performStepState(elems,xini,hsta);
    *tims = hini; // Remove for speed
    tx = dx;
        
    for (t=1;t<nstepsm;++t){
        *(hsta+t) = hini;
        performStepState(elems+tx,elems+tx-dx,hsta+t);
        *(tims+t) = *(tims+t-1)+hini; // Remove for speed
        tx += dx;
    }
        
    *(hsta+nstepsm) = hini;
    performStepState(xfin,elems+tx-dx,hsta+nstepsm);
    *(tims+nstepsm) = *(tims+nstepsm-1)+hini; // Remove for speed
 */
}

/**
** Adjoint integration main loop (backwards in time with varying step-size)
**/
template <class dynT> void ExplicitRKintegrator<dynT>::integrateAdjointAdaptive(){
   
  	std::cout<<"ERROR<ExplicitRKintegrator>: Adjoint integration with adaptive step-size not implemented."<<std::endl; 
    
 /*   assert(nsteps>=0);
    int t,ta;
    double *elems,*tims;
    elems = adjointGrid.getElems();
    tims = adjointGrid.getTime();
    
    *hadj = hini;
    performStepAdjoint(elems,adfin,hadj);
    *tims = hini; // Remove for speed
    ta = da;   
    for (t=1;t<nstepsm;++t){
        *(hadj+t) = hini;
        performStepAdjoint(elems+ta,elems+ta-da,hadj+t);
        *(tims+t) = *(tims+t-1)+hini; // Remove for speed
        ta += da;
    }
    *(hadj+nstepsm) = hini;
    performStepAdjoint(adini,elems+ta-da,hadj+nstepsm);
    *(tims+nstepsm) = *(tims+nstepsm-1)+hini;
    */      
}
