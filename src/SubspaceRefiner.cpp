//
//  SubspaceRefiner.cpp
//  MusAL
//
//  Created by Jean on 4/20/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//
#include <string.h>
#include <math.h>
#include <time.h>

#include "SubspaceRefiner.h"
#include "defaultOptions.h"

/**
 * Default constructor
 */
SubspaceRefiner::SubspaceRefiner(){
    
    diagHal = NULL;
    diagPrec = NULL;
    
    auxFul = NULL;
    vecFulIn = NULL;
    vecFulOut = NULL;
    
    zRed = NULL;
    rRed = NULL;
    pRed = NULL;
    yRed = NULL;
    xRed = NULL;
    
    zRedmin = NULL;
    zRedmax = NULL;
    zmin = NULL;
    zmax = NULL;
    
    gAL = NULL;
    hALshoot = NULL;
    hALpath = NULL;
    
    freeVars = NULL;
    
    z = NULL;
    zNex = NULL;
    zC = NULL;
    zeval = NULL;
    
    dZ = -1;
    dZBL = -1;
    dZshoot = -1;
    dZshootBL = -1;
    dZpath = -1;
    dZpathBL = -1;
    dHshoot = -1;
    dHpath = -1;
    
    Ns_ = -1;
    dx = -1;
    du = -1;
    dxu = -1;
    dxux = -1;
    dpc = -1;
    
    nfree = -1;
    dimRed = -1;
    dimRedBL = -1;
    
    z = NULL;
    zNex = NULL;
    zC = NULL;
    freeVars = NULL;

    activTol = ACTIV_TOL;
}

/**
** Copy constructor
**/
SubspaceRefiner::SubspaceRefiner(const SubspaceRefiner& obj) : qmod(obj.qmod){
 
    Ns_ = obj.Ns_;
    
    dZshoot = obj.dZshoot;
    dZshootBL = obj.dZshootBL;
    dZpath = obj.dZpath;
    dZpathBL = obj.dZpathBL;
    dZ = obj.dZ;
    dZBL = obj.dZBL;
    
    dHpath = obj.dHpath;
    dHshoot = obj.dHshoot;
    
    dx = obj.dx;
    du = obj.du;
    dxu = obj.dxu;
    dxux = obj.dxux;
    dpc = obj.dpc;
    
    maxIt = obj.maxIt;
    tol2 = obj.tol2;
    activTol = obj.activTol;
    
    nfree = obj.nfree;
    
    dimRed = obj.dimRed;
	dimRedBL = obj.dimRedBL;

    zmin = obj.zmin;
    zmax = obj.zmax;

    diagHal = new double[dZ];
    std::copy(obj.diagHal,obj.diagHal+obj.dZ,diagHal);
    
    diagPrec = new double[dZ];
    std::copy(obj.diagPrec,obj.diagPrec+obj.dZ,diagPrec);
    
    auxFul = new double[dZ];
    std::copy(obj.auxFul,obj.auxFul+obj.dZ,auxFul);
    
    vecFulIn = new double[dZ];
    std::copy(obj.vecFulIn,obj.vecFulIn+obj.dZ,vecFulIn);
    
    vecFulOut = new double[dZ];
    std::copy(obj.vecFulOut,obj.vecFulOut+obj.dZ,vecFulOut);
    
    zRed = new double[dZ];
    std::copy(obj.zRed,obj.zRed+obj.dZ,zRed);
    
    rRed = new double[dZ];
    std::copy(obj.rRed,obj.rRed+obj.dZ,rRed);
    
    pRed = new double[dZ];
    std::copy(obj.pRed,obj.pRed+obj.dZ,pRed);
    
    yRed = new double[dZ];
    std::copy(obj.yRed,obj.yRed+obj.dZ,yRed);
    
    xRed = new double[dZ];
    std::copy(obj.xRed,obj.xRed+obj.dZ,xRed);
    
    zRedmin = new double[dZ];
    std::copy(obj.zRedmin,obj.zRedmin+obj.dZ,zRedmin);
    
    zRedmax = new double[dZ];
    std::copy(obj.zRedmax,obj.zRedmax+obj.dZ,zRedmax);

    zeval = new double[dZ];
    std::copy(obj.zeval,obj.zeval+obj.dZ,zeval);
}

/**
** Assignment operator
**/
SubspaceRefiner& SubspaceRefiner::operator=(const SubspaceRefiner& obj){

	double *diagPrec_,*diagHal_,*auxFul_,*vecFulIn_,*vecFulOut_;
    double *zRed_,*rRed_,*pRed_,*yRed_,*xRed_;
    double *zRedmin_,*zRedmax_;
    double *zeval_;
    
    if (this != &obj){
        
        qmod = obj.qmod;

        dZpath = obj.dZpath;
        dZpathBL = obj.dZpathBL;
        
        dZshoot = obj.dZshoot;
        dZshootBL = obj.dZshootBL;
        
        dZ = obj.dZ;
        dZBL = obj.dZBL;
        
        dHshoot = obj.dHshoot;
        dHpath = obj.dHpath;
        
        Ns_ = obj.Ns_;
        
        dx = obj.dx;
        du = obj.du;
        dxu = obj.dxu;
        dxux = obj.dxux;
        dpc = obj.dpc;
        
        z = obj.z;
        zNex = obj.zNex;
        zC = obj.zC;
        
        freeVars = obj.freeVars;
        nfree = obj.nfree;
        
        dimRed = obj.dimRed;
        dimRedBL = obj.dimRedBL;
        
        gAL = obj.gAL;
        
        hALshoot = obj.hALshoot;
        hALpath = obj.hALpath;
        
        zmin = obj.zmin;
        zmax = obj.zmax;
        
        maxIt = obj.maxIt;
        tol2 = obj.tol2;
        activTol = obj.activTol;
        
        /* diagHal */
        if (diagHal != NULL)
            delete [] diagHal;
        diagHal_ = new double[dZ];
        std::copy(obj.diagHal,obj.diagHal+obj.dZ,diagHal_);
        diagHal = diagHal_;
        
        /* diagPrec */
        if (diagPrec != NULL)
            delete [] diagPrec;
        diagPrec_ = new double[dZ];
        std::copy(obj.diagPrec,obj.diagPrec+obj.dZ,diagPrec_);
        diagPrec = diagPrec_;
        
        /* auxFul */
        if (auxFul != NULL)
        	delete [] auxFul;
        auxFul_ = new double[dZ];
        std::copy(obj.auxFul,obj.auxFul+obj.dZ,auxFul_);
        auxFul = auxFul_;
        
        /* vecFulIn */
        if (vecFulIn != NULL)
        	delete [] vecFulIn;
        vecFulIn_ = new double[dZ];
        std::copy(obj.vecFulIn,obj.vecFulIn+obj.dZ,vecFulIn_);
        vecFulIn = vecFulIn_;
        
        /* vecFulOut */
        if (vecFulOut != NULL)
        	delete [] vecFulOut;
        vecFulOut_ = new double[dZ];	
        std::copy(obj.vecFulOut,obj.vecFulOut+obj.dZ,vecFulOut_);
        vecFulOut = vecFulOut_;
        
        /* zRed */
        if (zRed != NULL)
            delete [] zRed;
        zRed_ = new double[dZ];
        std::copy(obj.zRed,obj.zRed+obj.dZ,zRed_);
        zRed = zRed_;
        
        /* rRed */
        if (rRed != NULL)
            delete [] rRed;
        rRed_ = new double[dZ];
        std::copy(obj.rRed,obj.rRed+obj.dZ,rRed_);
        rRed = rRed_;
        
        /* pRed */
        if (pRed != NULL)
            delete [] pRed;
        pRed_ = new double[dZ];
        std::copy(obj.pRed,obj.pRed+obj.dZ,pRed_);
        pRed = pRed_;
        
        /* yRed */
        if (yRed != NULL)
            delete [] yRed;
        yRed_ = new double[dZ];
        std::copy(obj.yRed,obj.yRed+obj.dZ,yRed_);
        yRed = yRed_;
        
        /* xRed */
        if (xRed != NULL)
        	delete [] xRed;
        xRed_ = new double[dZ];
        std::copy(obj.xRed,obj.xRed+obj.dZ,xRed_);
        xRed = xRed_;
        
        /* zRedmin */
        if (zRedmin != NULL)
            delete [] zRedmin;
        zRedmin_ = new double[dZ];
        std::copy(obj.zRedmin,obj.zRedmin+obj.dZ,zRedmin_);
        zRedmin = zRedmin_;
        
        /* zRedmax */
        if (zRedmax != NULL)
            delete [] zRedmax;
        zRedmax_ = new double[dZ];
        std::copy(obj.zRedmax,obj.zRedmax+obj.dZ,zRedmax_);
        zRedmax = zRedmax_;

        /* zeval */
        if (zeval != NULL)
            delete [] zeval;
        zeval_ = new double[dZ];
        std::copy(obj.zeval,obj.zeval+obj.dZ,zeval_);
        zeval = zeval_;
    }
    
    return *this;
}

/**
** Destructor
**/
SubspaceRefiner::~SubspaceRefiner(){

    if (diagPrec != NULL)
        delete [] diagPrec;
    
    if (diagHal != NULL)
        delete [] diagHal;
    
    if (auxFul != NULL)
    	delete [] auxFul;
    	
    if (vecFulIn != NULL)	
    	delete [] vecFulIn;
    	
   	if (vecFulOut != NULL)
    	delete [] vecFulOut;
    
    if (zRed != NULL)
        delete [] zRed;
    
    if (rRed != NULL)
        delete [] rRed;
    
    if (pRed != NULL)
        delete [] pRed;
    
    if (yRed != NULL)
        delete [] yRed;
    
    if (xRed != NULL)
    	delete [] xRed;
    
    if (zRedmin != NULL)
        delete [] zRedmin;
    
    if (zRedmax != NULL)
        delete [] zRedmax;

    if (zeval != NULL)
        delete [] zeval;
}

/**
** Initialize local data
**/
void SubspaceRefiner::init(){

    /* diagHal */
    if (diagHal != NULL)
        delete [] diagHal;
    diagHal = new double[dZ];

    /* diagPrec */
    if (diagPrec != NULL)
        delete [] diagPrec;
    diagPrec = new double[dZ];

    /* auxFul */
    if (auxFul != NULL)
        delete [] auxFul;
    auxFul = new double[dZ];

    /* zRed */
    if (zRed != NULL)
        delete [] zRed;
    zRed = new double[dZ]; 
    
    /* rRed */
    if (rRed != NULL)
        delete [] rRed;
    rRed = new double[dZ];
    
    /* yRed */
    if (yRed != NULL)
        delete [] yRed;
    yRed = new double[dZ];

    /* pRed */
    if (pRed != NULL)
        delete [] pRed;
    pRed = new double[dZ];

    /* xRed */
    if (xRed != NULL)
        delete [] xRed;
    xRed = new double[dZ];
    
    /* vecFulIn */
    if (vecFulIn != NULL)
        delete [] vecFulIn;
    vecFulIn = new double[dZ];
    
    /* vecFulOut */
    if (vecFulOut != NULL)
        delete [] vecFulOut;
    vecFulOut = new double[dZ];
    
    /* zRedmin */
    if (zRedmin != NULL)
        delete [] zRedmin;
    zRedmin = new double[dZ];

    /* zRedmax */
    if (zRedmax != NULL)
        delete [] zRedmax;
    zRedmax = new double[dZ];

    /* zeval */
    if (zeval != NULL)
        delete [] zeval;
    zeval = new double[dZ];
}

/**
* Display shooting part of AL hessian 
*/
void SubspaceRefiner::displayHalShoot(){

    int n,i;
    double *hal;

    std::cout<<"SUBSPACE-REFINER<hALshoot>: "<<std::endl;    
    hal = hALshoot;
    for (n=0;n<Ns_;++n){
        std::cout<<"Block "<<n<<std::endl;
        for (i=0;i<dxux*dxux;++i){  
            if (i%dxux == dxux-1){
                std::cout<<*hal++<<std::endl;
            }
            else {
                std::cout<<*hal++<<" ";
            }
        }
    }
    std::cout<<" "<<std::endl;
}

/**
** Display diagonal of shooting part of AL hessian
**/
void SubspaceRefiner::displayDiagHalShoot(){

    int i;

    std::cout<<"SUBSPACE-REFINER<diagHalShoot>: "<<std::endl;
    for (i=0;i<dZshoot;++i)
        std::cout<<diagHal[i]<<std::endl;
    std::cout<<" "<<std::endl;
}

/**
** Display preconditionner for AL hessian
**/
void SubspaceRefiner::displayDiagPrec(){
 
    int i;

    std::cout<<"SUBSPACE-REFINER<diagPrec>:"<<std::endl;
    for (i=0;i<dimRed;++i)
        std::cout<<diagPrec[i]<<std::endl;
    std::cout<<" "<<std::endl;
}

/**
** Refinement method (TPCG on null-space of active constraints stored in freeVars)
**/
void SubspaceRefiner::findCandidate(const double trustRadius){

    int n,i,j,it;
    int *fr;
    double a1,a2,aux,auxr,auxy,ddd;
    double rh1,rh2,n2res,beta,px,curv;
    double zzr;
    double *yr,*pr,*rr,*zr,*xr;
    double *axf,*axfs,*axfp,*aaxfs;
	double *gg;
    double *hs,*hhs;
    double *hp,*hhp;
    double *dh,*dp;
    double *dhs,*dhp,*ddhs;
    double *zzs,*zzCs,*zzp,*zzCp;
    double *zzz,*zzzC,*zn;
    double *zm,*zrm,*zM,*zrM;
    double *vfis,*vvfis,*vfip,*vvfip,*vfos,*vvfos,*vfop;
    
//    std::cout<<"SubspaceRefiner findCandidate..."<<std::endl;
    
    /**
    * In the same loop, construct diagonal preconditionner and compute 1st half of reduced residual
     */
    /* Extract AL hessian diagonal and compute auxFul <- -1*H*(zC-z) */
    memset(diagHal,0,dZBL);
    memset(auxFul,0,dZBL);
	memset(vecFulIn,0,dZBL);

//    displayHalShoot();

    dhs = diagHal;
    dhp = dhs+dZshoot;
    hs = hALshoot;
    hhs = hALshoot;
    hp = hALpath;
    hhp = hALpath;
    axfs = auxFul;
    axfp = auxFul+dZshoot;
    zzs = z;
    zzCs = zC;
    zzp = zzs+dZshoot;
    zzCp = zzCs+dZshoot;
    for (n=0;n<Ns_;++n){           
        /* For shooting part */
        ddhs = dhs;
        aaxfs = axfs;
        for (i=0;i<dxux;++i){
        	/* Compute diagonal term */
            *ddhs += *(hs+i);
            hs += dxux;
            ++ddhs;
            /* Compute product against shooting block */
            zzz = zzs;
    		zzzC = zzCs;
            aux = 0.;
    		for (j=0;j<dxux;++j){
				aux += *hhs*(*zzz-*zzzC);
                ++hhs;
				++zzz;
				++zzzC;			
    		}
    		*aaxfs += aux;
    		++aaxfs;
        }
        dhs += dxu;
        axfs += dxu;
    	zzs += dxu;
    	zzCs += dxu;  
        /* For path-constraint part */
        for (i=0;i<dpc;++i){
        	/* Compute diagonal term */
            *dhp += *(hp+i);
            hp += dpc;
            ++dhp;
            /* Compute product against path-constraint block */
            zzz = zzp;
    		zzzC = zzCp;
            aux = 0.;
    		for (j=0;j<dpc;++j){
    			aux += *hhp*(*zzz-*zzzC);
    			++hhp;
    			++zzz;
    			++zzzC;
    		}
    		*axfp += aux;
    		++axfp;
        }		
		zzp += dpc;
    	zzCp += dpc;
    }
//    displayDiagHalShoot();


    /* Find smallest positive element of diagHal */
/*  double mndg;  
    dh = diagHal;
    mndg = POSINF;
    fr = freeVars;
    for (i=0;i<dZ;++i){
        aux = *dh;
        if ((*fr==1)&&(aux>=POSZER))
            mndg = MIN(mndg,aux);
        ++dh;
        ++fr;
    } */
//    std::cout<<"SUBSPACE-REFINER<mndg>="<<mndg<<std::endl;
 

    /** - Compute reduced residual
    *   - Apply prec. to reduced residual and store in reduced update direction 
    *   - Initialize reduced bounds, optimizer and CG direction 
    */
    fr = freeVars;
    axf = auxFul;
    dh = diagHal;
    dp = diagPrec;
    gg = gAL;
    zm = zmin;
    zM = zmax;
    zrm = zRedmin;
    zrM = zRedmax;
    rr = rRed;
    zr = zRed;
    yr = yRed;
    zzzC = zC;
    for (i=0;i<dZ;++i){
        if (*fr==1){
            ddd = *axf-*gg;
            /* Fill-in reduced residual */
            *rr = ddd; 
            //aux = MAX(*dh,mndg);
            aux = MAX(*dh,pr_reg_coef);
            //aux = 1.;
            //aux = *dh;
            /* Apply preconditionner and fill in reduced update direction */
            *yr = ddd/aux;
            /* Store preconditionner */
            *dp = aux;
            /* Init. reduced bounds */
            *zrm = *zm;
            *zrM = *zM;
            /* Init. optimizer on Cauchy point */ 
            *zr = *zzzC;
            ++rr;
            ++yr;
            ++dp;
            ++zrm;
            ++zrM;
            ++zr;
        }
        ++fr;
        ++dh;
        ++gg;
        ++axf;
        ++zm;
        ++zM;
        ++zzzC;
    }    
    //displayDiagPrec();

    /* Init. CG direction */
    memset(pRed,0,dimRedBL);
        
    /**
    * Safeguarded TPCG main loop
     */
    rh1 = 1.;
    /* Compute squared 2-norm of model residual n2res and preconditioned residual rh2 */
    rh2 = 0.;
    n2res = 0.;
    rr = rRed; 
    yr = yRed;
	for (i=0;i<dimRed;++i){
        auxr = *rr++;
		rh2 += *yr++*auxr;
        n2res += auxr*auxr;
    }

    it = 0;
    while (1) {
//        std::cout<<"SUBSPACE-REFINER iteration: "<<it<<std::endl;
//        std::cout<<"Reduced subspace dimension: "<<dimRed<<std::endl;
//        std::cout<<"n2res="<<n2res<<", rh2="<<rh2<<std::endl;
//        std::cout<<"tol2="<<tol2<<std::endl;
//        std::cout<<"Model value ="<<evalModelReducedCG(zRed)<<std::endl;

        /* Check tolerance and # iterations, stop if satisfied */
        if ((n2res<=tol2)||(it>=maxIt)){
            zr = zRed;
            fr = freeVars;
            zn = zNex;
            for (i=0;i<dZ;++i){
                if (*fr==1){
                    *zn = *zr;
                    ++zr;
                }
                ++zn;
                ++fr;
            }
        //    std::cout<<"SUBSPACE-REFINER: tolerance or max. num. iters hit."<<std::endl;
        //    std::cout<<"n2res="<<n2res<<", tol2="<<tol2<<std::endl;
        //    std::cout<<"it="<<it<<", maxIt="<<maxIt<<std::endl;
            break;
        }
        
        /* Step-size computation */
        beta = rh2/rh1;
        
        /* Conjugate direction update pr <- yr + beta*pr */
        pr = pRed;
        yr = yRed;
        for (i=0;i<dimRed;++i){
            aux = *pr;
            *pr = *yr+beta*aux;
            ++pr;
            ++yr;
        }
        
        /**
        * Product with reduced AL hessian xRed <- Hred*pRed
        */
        /* Map reduced conjugate direction pRed into vecFulIn */
 		fr = freeVars;
 		vfis = vecFulIn;
 		pr = pRed;
 		for (i=0;i<dZ;++i){  
 			if (*fr==1){
 				*vfis = *pr;
                ++pr;
            }
 			++fr;
			++vfis;
        }

        /* Product with AL hessian vecFulOut <- H*vecFulIn */
        memset(vecFulOut,0,dZBL);
        
        hs = hALshoot;
        hp = hALpath;
        vfos = vecFulOut;
        vfop = vfos+dZshoot;
        vfis = vecFulIn;
        vfip = vfis+dZshoot;
        for (n=0;n<Ns_;++n){
        	
        	/* Shooting block */
        	vvfos = vfos;
        	for (i=0;i<dxux;++i){
        		vvfis = vfis;
        		aux = 0.;
        		for (j=0;j<dxux;++j){
        			aux += *hs**vvfis;
        			++hs;
        			++vvfis;
        		}
        		*vvfos += aux;
        		++vvfos;
        	}
        	vfos += dxu;
        	vfis += dxu;
        	
        	/* Path-constraint block */
			for (i=0;i<dpc;++i){
				vvfip = vfip;			
				aux = 0.;
				for (j=0;j<dpc;++j){        	
        			aux += *hp**vvfip;
        			++hp;
        			++vvfip;
        		}
        		*vfop = aux;
        		++vfop;
        	}
        	vfip += dpc;
        }
        
        /* Map vecFulOut back into xRed */
        fr = freeVars;
        vfos = vecFulOut;
        xr = xRed;
        for (i=0;i<dZ;++i){
        	if (*fr==1){
        		*xr = *vfos;
                ++xr;
            }
        	++fr;
        	++vfos;
        }
                
        /*  1) Compute curvature
            2) Find largest step-size such that (trust-region or set) bound is hit */
        a1 = POSINF;
        curv = 0.;
        zrm = zRedmin;
        zrM = zRedmax;
        pr = pRed;
        zr = zRed;
        xr = xRed;
        for (i=0;i<dimRed;++i){
            zzr = *zr;
            px = *pr++;
            curv += *xr++*px;

		    if (px>POSZER){
                a1 = MIN((MIN(*zrM,zzr+trustRadius)-zzr)/px,a1);
                ++zr;
                ++zrm;
                ++zrM;
                continue;
            }
            
            if (px<NEGZER){
                a1 = MIN((MAX(*zrm,zzr-trustRadius)-zzr)/px,a1);
                ++zr;
                ++zrm;
                ++zrM;
                continue;
            }
            
        //    a1 = MIN(POSINF,a1);

            ++zr;
            ++zrm;
            ++zrM;
        }   
        
        /**
        * Curvature check and update
        */
        /* If nonpositive curvature encountered, move to boundary and stop */
        if (curv<=POSZER){
            zr = zRed;
            pr = pRed;
            fr = freeVars;
            zn = zNex;
            for (i=0;i<dZ;++i){
                if (*fr==1){
                    *zn = *zr+*pr*a1;
                    ++zr;
                    ++pr;
                }
                ++zn;
                ++fr;
            }
        //    std::cout<<"SUBSPACE-REFINER: nonpositive curvature after "<<it<<" iters, curv="<<curv<<std::endl;
            break;
        }
        
        a2 = rh2/curv;
        
        /* If CG step-size too large, hit constraint and stop */
        if (a2>a1){
            zr = zRed;
            pr = pRed;
            fr = freeVars;
            zn = zNex;
            for (i=0;i<dZ;++i){
                if (*fr==1){
                    *zn = *zr+*pr*a1;
                    ++zr;
                    ++pr;
                }
                ++zn;
                ++fr;
            }
        //    std::cout<<"SUBSPACE-REFINER: constraint hit after "<<it<<" iters."<<std::endl;
            break;
        }
        
        /* PCG update */
        zr = zRed;
        pr = pRed;
        xr = xRed;
        rr = rRed;
        yr = yRed;
        dp = diagPrec;
        rh1 = rh2;
        rh2 = 0.;
        n2res = 0.;
        for (i=0;i<dimRed;++i){
            *zr += *pr*a2;
            auxr = *rr-*xr*a2;
            *rr = auxr;
            auxy = auxr/(*dp);
            *yr = auxy;

            rh2 += auxy*auxr;
            n2res += auxr*auxr;

            ++zr;
            ++xr;
            ++pr;
            ++rr;
            ++yr;
            ++dp;
        }

        ++it;
    }
    
}

/**
** Refinement method (TPCG on null-space of active constraints stored in freeVars) with restart if problem constraint encountered
**/
void SubspaceRefiner::findCandidateRestart(const double trustRadius){

    int n,i,j,it,nwlb,nwub;
    int *fr;
    double a1,a2,aux,auxr,auxy,ddd;
    double rh1,rh2,n2res,beta,px,curv;
    double ax,al,au;
    double zzr;
    double *yr,*pr,*rr,*zr,*xr;
    double *axf,*axfs,*axfp,*aaxfs;
    double *gg;
    double *hs,*hhs;
    double *hp,*hhp;
    double *dh,*dp;
    double *dhs,*dhp,*ddhs;
    double *zzs,*zzCs,*zzp,*zzCp;
    double *zzz,*zzzC,*zzn;
    double *zn,*zns,*znp;
    double *zm,*zrm,*zM,*zrM,*zmm,*zMM;
    double *vfis,*vvfis,*vfip,*vvfip,*vfos,*vvfos,*vfop;

/*    clock_t str,end;
    double multime;
    double maptime;

    multime = 0.;
    maptime = 0.; */
    
//    std::cout<<"SubspaceRefiner findCandidate..."<<std::endl;
    
    /**
    * In the same loop, construct diagonal preconditionner and compute 1st half of reduced residual
    */
    /* Extract AL hessian diagonal and compute auxFul <- -1*H*(zC-z) */
    memset(diagHal,0,dZBL);
    memset(auxFul,0,dZBL);
    memset(vecFulIn,0,dZBL);

    dhs = diagHal;
    dhp = dhs+dZshoot;
    hs = hALshoot;
    hhs = hALshoot;
    hp = hALpath;
    hhp = hALpath;
    axfs = auxFul;
    axfp = auxFul+dZshoot;
    zzs = z;
    zzCs = zC;
    zzp = zzs+dZshoot;
    zzCp = zzCs+dZshoot;
    for (n=0;n<Ns_;++n){    
        /* Shooting part */
        ddhs = dhs;
        aaxfs = axfs;
        for (i=0;i<dxux;++i){
            /* Compute diagonal term */
            *ddhs += *(hs+i);
            hs += dxux;
            ++ddhs;
            /* Compute product against shooting block */
            zzz = zzs;
            zzzC = zzCs;
            aux = 0.;
            for (j=0;j<dxux;++j){
                aux += *hhs*(*zzz-*zzzC);
                ++hhs;
                ++zzz;
                ++zzzC;         
            }
            *aaxfs += aux;
            ++aaxfs;
        }
        dhs += dxu;
        axfs += dxu;
        zzs += dxu;
        zzCs += dxu;
                
        /* Path-constraint part */
        for (i=0;i<dpc;++i){
            /* Compute diagonal term */
            *dhp += *(hp+i);
            hp += dpc;
            ++dhp;
            /* Compute product against path-constraint block */
            zzz = zzp;
            zzzC = zzCp;
            aux = 0.;
            for (j=0;j<dpc;++j){
                aux += *hhp*(*zzz-*zzzC);
                ++hhp;
                ++zzz;
                ++zzzC;
            }
            *axfp += aux;
            ++axfp;
        }       
        zzp += dpc;
        zzCp += dpc;
    }
//    displayDiagHalShoot();

    /* Find smallest positive element of diagHal */
/*  double mndg;  
    dh = diagHal;
    mndg = POSINF;
    fr = freeVars;
    for (i=0;i<dZ;++i){
        aux = *dh;
        if ((*fr==1)&&(aux>=POSZER))
            mndg = MIN(mndg,aux);
        ++dh;
        ++fr;
    } */
//    std::cout<<"SUBSPACE-REFINER<mndg>="<<mndg<<std::endl; 

    /** - Compute reduced residual
    *   - Apply prec. to reduced residual and store in reduced update direction 
    *   - Initialize reduced bounds, optimizer and CG direction 
    */
    fr = freeVars;
    axf = auxFul;
    dh = diagHal;
    dp = diagPrec;
    gg = gAL;
    zm = zmin;
    zM = zmax;
    zrm = zRedmin;
    zrM = zRedmax;
    rr = rRed;
    zr = zRed;
    yr = yRed;
    zzzC = zC;
    for (i=0;i<dZ;++i){
        if (*fr==1){
            ddd = *axf-*gg;
            /* Fill-in reduced residual */
            *rr = ddd; 
            //aux = MAX(*dh,mndg);
            aux = MAX(*dh,pr_reg_coef);
            //aux = 1.;
            //aux = *dh;
            /* Apply preconditionner and fill in reduced update direction */
            *yr = ddd/aux;
            /* Store preconditionner */
            *dp = aux;
            /* Init. reduced bounds */
            *zrm = *zm;
            *zrM = *zM;
            /* Init. optimizer on Cauchy point */ 
            *zr = *zzzC;
            ++rr;
            ++yr;
            ++dp;
            ++zrm;
            ++zrM;
            ++zr;
        }
        ++fr;
        ++dh;
        ++gg;
        ++axf;
        ++zm;
        ++zM;
        ++zzzC;
    }    

    /* Init. CG direction */
    memset(pRed,0,dimRedBL);
            
    /* Compute squared 2-norm of model residual n2res and preconditioned residual rh2 */        
    rh1 = 1.;
    rh2 = 0.;
    n2res = 0.;
    rr = rRed; 
    yr = yRed;
    for (i=0;i<dimRed;++i){
        auxr = *rr++;
        rh2 += *yr++*auxr;
        n2res += auxr*auxr;
    }
    /**
    * Safeguarded TPCG main loop
     */
    it = 0;
    while (1) {
//        std::cout<<"SUBSPACE-REFINER iteration: "<<it<<std::endl;
//        std::cout<<"Reduced subspace dimension: "<<dimRed<<std::endl;
//        std::cout<<"n2res="<<n2res<<", rh2="<<rh2<<std::endl;
//        std::cout<<"tol2="<<tol2<<std::endl;

        /* Check tolerance and # iterations, stop if satisfied */
        //if (it>=100){
        if ((n2res<=tol2)||(it>=maxIt)){
            zr = zRed;
            fr = freeVars;
            zn = zNex;
            for (i=0;i<dZ;++i){
                if (*fr==1){
                    *zn = *zr;
                    ++zr;
                }
                ++zn;
                ++fr;
            }
        //    std::cout<<"SUBSPACE-REFINER: tolerance or max. num. iters hit."<<std::endl;
        //    std::cout<<"n2res="<<n2res<<", tol2="<<tol2<<std::endl;
        //    std::cout<<"it="<<it<<", maxIt="<<maxIt<<std::endl;
            break;
        }
        
        /* Step-size computation */
        beta = rh2/rh1;
        
        /* Conjugate direction update pr <- yr + beta*pr */
        pr = pRed;
        yr = yRed;
        for (i=0;i<dimRed;++i){
            aux = *pr;
            *pr = *yr+beta*aux;
            ++pr;
            ++yr;
        }
        
        /**
        * Product with reduced AL hessian xRed <- Hred*pRed
        */
    //    str = clock();
        /* Map reduced conjugate direction pRed into vecFulIn */
        fr = freeVars;
        vfis = vecFulIn;
        pr = pRed;
        for (i=0;i<dZ;++i){  
            if (*fr==1){
                *vfis = *pr;
                ++pr;
            }
            ++fr;
            ++vfis;
        }
    //    end = clock();
    //    maptime += (double)(end-str)/CLOCKS_PER_SEC;

    //    str = clock();
        /* Product with AL hessian vecFulOut <- H*vecFulIn */
        memset(vecFulOut,0,dZBL);
        
        hs = hALshoot;
        hp = hALpath;
        vfos = vecFulOut;
        vfop = vfos+dZshoot;
        vfis = vecFulIn;
        vfip = vfis+dZshoot;
        for (n=0;n<Ns_;++n){
            /* Shooting block */
            vvfos = vfos;
            for (i=0;i<dxux;++i){
                vvfis = vfis;
                aux = 0.;
                for (j=0;j<dxux;++j){
                    aux += *hs**vvfis;
                    ++hs;
                    ++vvfis;
                }
                *vvfos += aux;
                ++vvfos;
            }
            vfos += dxu;
            vfis += dxu;
            /* Path-constraint block */
            for (i=0;i<dpc;++i){
                vvfip = vfip;           
                aux = 0.;
                for (j=0;j<dpc;++j){            
                    aux += *hp**vvfip;
                    ++hp;
                    ++vvfip;
                }
                *vfop = aux;
                ++vfop;
            }
            vfip += dpc;
        }
    //    end = clock();
    //    multime += (double)(end-str)/CLOCKS_PER_SEC;
        
    //    str = clock();
        /* Map vecFulOut back into xRed */
        fr = freeVars;
        vfos = vecFulOut;
        xr = xRed;
        for (i=0;i<dZ;++i){
            if (*fr==1){
                *xr = *vfos;
                ++xr;
            }
            ++fr;
            ++vfos;
        }
    //    end = clock();
    //    maptime += (double)(end-str)/CLOCKS_PER_SEC;

        /**  1) Compute curvature
             2) Find largest step-size such that (trust-region or NLP) bound is hit */
        a1 = POSINF;
        curv = 0.;
        zrm = zRedmin;
        zrM = zRedmax;
        pr = pRed;
        zr = zRed;
        xr = xRed;
        for (i=0;i<dimRed;++i){
            zzr = *zr;
            px = *pr;
            curv += *xr*px;

            if (px>POSZER)
                aux = (MIN(*zrM,zzr+trustRadius)-zzr)/px;
            
            if (px<NEGZER)
                aux = (MAX(*zrm,zzr-trustRadius)-zzr)/px;
            
            if ((NEGZER<=px)&&(px<=POSZER))
                aux = POSINF;
        
            a1 = MIN(aux,a1);

            ++xr;
            ++pr;
            ++zr;
            ++zrm;
            ++zrM;
        }
        
        /**
        * Curvature check and update
        */
        /* If nonpositive curvature encountered, move to boundary and stop */
        if (curv<=POSZER){
            zr = zRed;
            pr = pRed;
            fr = freeVars;
            zn = zNex;
            for (i=0;i<dZ;++i){
                if (*fr==1){
                    *zn = *zr+*pr*a1;
                    ++zr;
                    ++pr;
                }
                ++zn;
                ++fr;
            }
        //    std::cout<<"SUBSPACE-REFINER: nonpositive curvature after "<<it<<" iters, curv="<<curv<<std::endl;
            break;
        }
    
        a2 = rh2/curv;
        
        /** If CG step-size too large, constraint hit:
            - case 1: trust-region bound => stop
            - case 2: NLP bounds => restart on updated subspace **/
        if (a2>a1){
//            std::cout<<"CG step-size too large."<<std::endl;
//            std::cout<<"a2="<<a2<<", a1="<<a1<<std::endl;
            zr = zRed;
            pr = pRed;
            zn = zNex;
            fr = freeVars;
            zmm = zmin;
            zMM = zmax;
            nwlb = 0;
            nwub = 0;
            for (i=0;i<dZ;++i){
                if (*fr==1){
                    ax = *zr+*pr*a1;
                    al = *zmm;
                    au = *zMM;
                //    std::cout<<"al="<<al<<std::endl;
                //    std::cout<<"au="<<au<<std::endl;

                    /* Lower bound activated */
                    if (ABS(ax-al)<=activTol){
                    //    std::cout<<"Lower bound activated"<<std::endl;
                    //    std::cout<<"al="<<al<<", ax="<<ax<<std::endl;
                        *zn = al;
                        *fr = 0; // Fix variable
                        ++nwlb;

                        ++zr;
                        ++pr;
                        ++zmm;
                        ++zMM;
                        ++zn;
                        ++fr;
                        continue;
                    }

                    /* Upper bound activated */
                    if (ABS(ax-au)<=activTol){
                    //    std::cout<<"Upper bound activated"<<std::endl;
                    //    std::cout<<"au="<<au<<", ax="<<ax<<std::endl;
                        *zn = au;
                        *fr = 0; // Fix variable 
                        ++nwub;

                        ++zr;
                        ++pr;
                        ++zmm;
                        ++zMM;
                        ++zn;
                        ++fr;
                        continue;
                    }    
                    *zn = ax;
                    ++zr;
                    ++pr;
                }
                ++zmm;
                ++zMM;
                ++zn;
                ++fr;
            }

            if ((nwlb==0)&&(nwub==0)){
            //    std::cout<<"SUBSPACE-REFINER: TR constraint hit after "<<it<<" iters."<<std::endl;
                break;
            }

            /** 
             * Prepare for restart with updated freeVars & zNex 
             */
             /* Update dimension of reduced subspace */
            dimRed -= nwlb;
            dimRed -= nwub;
            dimRedBL = dimRed*SZDBL;
            nfree = dimRed;
            maxIt = nfree;

            /* Reset auxiliary vectors to zero */
            memset(auxFul,0,dZBL);
            memset(vecFulIn,0,dZBL);

            /* Recompute hessian product part of residual */
            hs = hALshoot;
            hp = hALpath;
            axfs = auxFul;
            axfp = axfs+dZshoot;
            zzs = z;
            zns = zNex;
            zzp = zzs+dZshoot;
            znp = zns+dZshoot;
            for (n=0;n<Ns_;++n){  
                /* For shooting part */
                aaxfs = axfs;
                for (i=0;i<dxux;++i){
                    /* Compute product against shooting block */
                    zzz = zzs;
                    zzn = zns;
                    aux = 0.;
                    for (j=0;j<dxux;++j){
                        aux += *hs*(*zzz-*zzn);
                        ++hs;
                        ++zzz;
                        ++zzn;         
                    }
                    *aaxfs += aux;
                    ++aaxfs;
                }
                axfs += dxu;
                zzs += dxu;
                zns += dxu;
                
                /* For path-constraint part */
                for (i=0;i<dpc;++i){
                    /* Compute product against path-constraint block */
                    zzz = zzp;
                    zzn = znp;
                    aux = 0.;
                    for (j=0;j<dpc;++j){
                        aux += *hp*(*zzz-*zzn);
                        ++hp;
                        ++zzz;
                        ++zzn;
                    }
                    *axfp += aux;
                    ++axfp;
                }       
                zzp += dpc;
                znp += dpc;
            }

            /** - Compute reduced residual
            *   - Apply prec. to reduced residual and store in reduced update direction 
            *   - Initialize reduced bounds, optimizer and CG direction 
            */
            fr = freeVars;
            axf = auxFul;
            dh = diagHal;
            dp = diagPrec;
            gg = gAL;
            zn = zNex;
            zm = zmin;
            zM = zmax;
            zrm = zRedmin;
            zrM = zRedmax;
            rr = rRed;
            zr = zRed;
            yr = yRed;
            for (i=0;i<dZ;++i){
                if (*fr==1){
                    ddd = *axf-*gg;
                    /* Fill-in reduced residual */
                    *rr = ddd; 
                    //aux = MAX(*dh,mndg);
                    aux = MAX(*dh,pr_reg_coef);
                    //aux = 1.;
                    //aux = *dh;
                    /* Apply preconditionner and fill in reduced update direction */
                    *yr = ddd/aux;
                    /* Store preconditionner */
                    *dp = aux;
                    /* Init. reduced bounds */
                    *zrm = *zm;
                    *zrM = *zM;
                    /* Init. reduced optimizer on free variables of last full optimizer */ 
                    *zr = *zn;
                    ++rr;
                    ++yr;
                    ++dp;
                    ++zrm;
                    ++zrM;
                    ++zr;
                }
                ++fr;
                ++dh;
                ++gg;
                ++axf;
                ++zm;
                ++zM;
                ++zn;
            }    

            /* Init. CG direction */
            memset(pRed,0,dimRedBL);
        
            rh1 = 1.;
            /* Compute squared 2-norm of model residual n2res and preconditioned residual rh2 */
            rh2 = 0.;
            n2res = 0.;
            rr = rRed; 
            yr = yRed;
            for (i=0;i<dimRed;++i){
            //    std::cout<<"rRed="<<rRed[i]<<std::endl;
            //    std::cout<<"yRed="<<yRed[i]<<std::endl;
                auxr = *rr++;
                rh2 += *yr++*auxr;
                n2res += auxr*auxr;
            }

            it = 0;

            /* Restart */
//            std::cout<<"SUBSPACE-REFINER<RESTART>"<<std::endl;
//            std::cout<<"rh2="<<rh2<<", n2res="<<n2res<<std::endl;
//            std::cout<<"nwlb="<<nwlb<<std::endl;
//            std::cout<<"nwub="<<nwub<<std::endl;
            continue;
        }
        /* PCG update */
        zr = zRed;
        pr = pRed;
        xr = xRed;
        rr = rRed;
        yr = yRed;
        dp = diagPrec;
        rh1 = rh2;
        rh2 = 0.;
        n2res = 0.;
        for (i=0;i<dimRed;++i){
            *zr += *pr*a2;
            auxr = *rr-*xr*a2;
            *rr = auxr;
            auxy = auxr/(*dp);
            *yr = auxy;

            rh2 += auxy*auxr;
            n2res += auxr*auxr;

            ++zr;
            ++xr;
            ++pr;
            ++rr;
            ++yr;
            ++dp;
        }
        ++it;
    }
//    std::cout<<"maptime ="<<maptime<<std::endl;
//    std::cout<<"multime = "<<multime<<std::endl;
}

/**
* Compute model value at current TPCG iterate (not usable in a restarted version)
*/
double SubspaceRefiner::evalModelReducedCG(double* zCG){

    int i,j,n,*fr;
    double vmodl,vmodq,vmod,dyz;
    double *zr,*zn;
    double *hsn,*hpn;
    double *ys,*yp,*zs,*zp;

    memcpy(zeval,zC,dZBL);

    zr = zCG; // Pointer to reduced-space vector
    zn = zeval; // Pointer to full-space vector
    fr = freeVars; // Pointer to free variables
    for (i=0;i<dZ;++i){
        if (*fr==1){
            *zn = *zr;
            ++zr;
        }
        ++zn;
        ++fr;
    }

    /* Evaluate gAL'*(zCG-z) */
    vmodl = 0.;
    for (i=0;i<dZ;++i){
        vmodl += gAL[i]*(zeval[i]-z[i]);
    }

    /* Evaluate 0.5*(zCG-z)'*hAL*(zCG-z)*/
    hsn = hALshoot;
    hpn = hALpath;
    
    ys = zeval;
    yp = zeval+dZshoot;
    
    zs = z;
    zp = z+dZshoot;

    vmodq = 0.;
    for (n=0;n<Ns_;++n){
    
        /* Shooting block of approx AL hessian */
        for (i=0;i<dxux;++i){
            
            dyz = *(ys+i)-*(zs+i);
            
            /* Hessian cross product */
            for (j=0;j<dxux;++j){
                vmodq += *hsn*dyz*(*(ys+j)-*(zs+j));
                ++hsn;
            }
        }
        ys += dxu;
        zs += dxu;
        
        /* Path constraint block of approx AL hessian */
        for (i=0;i<dpc;++i){
            
            dyz = *(yp+i)-*(zp+i);
            
            /* Hessian cross product */
            for (j=0;j<dpc;++j){
                vmodq += *hpn*dyz*(*(yp+j)-*(zp+j));
                ++hpn;
            }
        }
        yp += dpc;
        zp += dpc;
    }
    
    vmod = vmodl+0.5*vmodq;

    return vmod;
}


/**
* Real-time refinement method (stopped after fixed # iterations)
*/
void SubspaceRefiner::findCandidateRealTime(const double trustRadius){
    
    








}