#include <iostream>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "PreconRefine.h"
#include "defaultOptions.h"
#include "userData.h"

/**
* Default constructor
*/
PreconRefine::PreconRefine(){

	z = NULL;
    zNex = NULL;
    
    gAL = NULL;
    
    prdHss = NULL;

    hALshoot = NULL;
    hALpath = NULL;
    hALshootRed = NULL;
    hALpathRed = NULL;
    
    zmin = NULL;
    zmax = NULL;

    zRed = NULL; 
    rRed = NULL; 
    yRed = NULL; 
    pRed = NULL; 
   	xRed = NULL; 

    zRedmin = NULL;
    zRedmax = NULL;

    freeVars = NULL;
    indFreeVars = NULL;
    varStatus = NULL;
    ivarFreeState = NULL;
    ivarFreeInput = NULL;
    ivarFreeSlack = NULL;
    ivarFreeGroupShoot = NULL;
    nfreeState = NULL;
    nfreeInput = NULL;
    nfreeSlack = NULL;
    nfreeGroupShoot = NULL;

    activTol = ACTIV_TOL;
    pr_reg_coef = PREC_REG_COEF;
    bwd = BAND_WIDTH;

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
    
    dimRed = -1;
    dimRedBL = -1;
    dimRedShoot = -1;
    dimRedPath = -1;

    cumCG = -1;
    cumRestart = -1;
}

/**
** Destructor 
**/
PreconRefine::~PreconRefine(){

	if (zRed != NULL)
    	delete [] zRed;

    if (rRed != NULL)
    	delete [] rRed;

    if (yRed != NULL)
    	delete [] yRed;
    
    if (pRed != NULL)
    	delete [] pRed;
    
    if (xRed != NULL)
    	delete [] xRed; 

    if (zRedmin != NULL)
    	delete [] zRedmin;
    
    if (zRedmax != NULL)
    	delete [] zRedmax;
}

/**
** Copy constructor
**/
PreconRefine::PreconRefine(const PreconRefine& obj) : qmod(obj.qmod),precon(obj.precon){
    
    dZ = obj.dZ;
    dZBL = obj.dZBL;
    
    dZshoot = obj.dZshoot;
   	dZshootBL = obj.dZshootBL;
    dZpath = obj.dZpath;
    dZpathBL = obj.dZpath;
    
    dHshoot = obj.dHshoot;
    dHpath = obj.dHpath;
    dhS = obj.dhS;
    dhP = obj.dhP;
    
    Ns_ = obj.Ns_;
    Ns_BL = obj.Ns_BL;
    NsBL = obj.NsBL;
    Ns_NT = obj.Ns_NT;
    NsNT = obj.NsNT;
    
    dx = obj.dx;
    du = obj.du;
    dxu = obj.dxu;
    dxux = obj.dxux;
    dpc = obj.dpc;
    
    dimRed = obj.dimRed;
    dimRedBL = obj.dimRed;
    dimRedShoot = obj.dimRedShoot;
    dimRedPath = obj.dimRedPath;

    z = obj.z;
    zNex = obj.zNex;
    
    gAL = obj.gAL;
    
    prdHss = obj.prdHss;

    hALshoot = obj.hALshoot;
    hALpath = obj.hALpath;
    hALshootRed = obj.hALshootRed;
    hALpathRed = obj.hALpathRed;
    
    zmin = obj.zmin;
    zmax = obj.zmax;

    freeVars = obj.freeVars; 
    indFreeVars = obj.indFreeVars;
    varStatus = obj.varStatus;
    ivarFreeState = obj.ivarFreeState;
    ivarFreeInput = obj.ivarFreeInput;
    ivarFreeSlack = obj.ivarFreeSlack;
    ivarFreeGroupShoot = obj.ivarFreeGroupShoot;
    nfreeState = obj.nfreeState;
    nfreeInput = obj.nfreeInput;
    nfreeSlack = obj.nfreeSlack;
    nfreeGroupShoot = obj.nfreeGroupShoot;

    maxIt = obj.maxIt;
    maxCumCG = obj.maxCumCG;
    pr_reg_coef = obj.pr_reg_coef;
   	activTol = obj.activTol;
    bwd = obj.bwd;

    cumCG = obj.cumCG;
    cumRestart = obj.cumRestart;

    zRed = new double[dZ];
    std::copy(obj.zRed,obj.zRed+dZ,zRed);

    rRed = new double[dZ];
    std::copy(obj.rRed,obj.rRed+dZ,rRed);

    yRed = new double[dZ]; 
    std::copy(obj.yRed,obj.yRed+dZ,yRed);

    pRed = new double[dZ];
    std::copy(obj.pRed,obj.pRed+dZ,pRed);

    xRed = new double[dZ]; 
    std::copy(obj.xRed,obj.xRed+dZ,xRed);

    zRedmin = new double[dZ];
    std::copy(obj.zRedmin,obj.zRedmin+dZ,zRedmin);
    
    zRedmax = new double[dZ];
    std::copy(obj.zRedmax,obj.zRedmax+dZ,zRedmax);
}

/**
* Assignment operator
*/
PreconRefine& PreconRefine::operator= (const PreconRefine& obj){

	double *diagHal_,*diagPrec_;
	double *auxFul_;
	double *zRed_,*rRed_,*yRed_,*pRed_,*xRed_;
	double *zRedmin_,*zRedmax_;

	if (this != &obj){

		qmod = obj.qmod;
        precon = obj.precon;
    
    	dZ = obj.dZ;
    	dZBL = obj.dZBL;
    	dZshoot = obj.dZshoot;
   		dZshootBL = obj.dZshootBL;
    	dZpath = obj.dZpath;
    	dZpathBL = obj.dZpath;
    
        prdHss = obj.prdHss;

    	dHshoot = obj.dHshoot;
    	dHpath = obj.dHpath;
        dhS = obj.dhS;
        dhP = obj.dhP;
    
    	Ns_ = obj.Ns_;
        Ns_BL = obj.Ns_BL;
        NsBL = obj.NsBL;
        Ns_NT = obj.Ns_NT;
        NsNT = obj.NsNT;
    
    	dx = obj.dx;
    	du = obj.du;
    	dxu = obj.dxu;
    	dxux = obj.dxux;
    	dpc = obj.dpc;
    
    	dimRed = obj.dimRed;
    	dimRedBL = obj.dimRed;
        dimRedShoot = obj.dimRedShoot;
        dimRedPath = obj.dimRedPath;

        cumCG = obj.cumCG;
        cumRestart = obj.cumRestart;

    	z = obj.z;
    	zNex = obj.zNex;
    
    	gAL = obj.gAL;
    
    	hALshoot = obj.hALshoot;
    	hALpath = obj.hALpath;
        hALshootRed = obj.hALshootRed;
        hALpathRed = obj.hALpathRed;
    
    	zmin = obj.zmin;
    	zmax = obj.zmax;

        freeVars = obj.freeVars;
        indFreeVars = obj.indFreeVars;
        varStatus = obj.varStatus;
        ivarFreeState = obj.ivarFreeState;
        ivarFreeInput = obj.ivarFreeInput;
        ivarFreeSlack = obj.ivarFreeSlack;
        ivarFreeGroupShoot = obj.ivarFreeGroupShoot;
        nfreeState = obj.nfreeState;
        nfreeInput = obj.nfreeInput;
        nfreeSlack = obj.nfreeSlack;
        nfreeGroupShoot = obj.nfreeGroupShoot;

        maxCumCG = obj.maxCumCG;
    	maxIt = obj.maxIt;
    	pr_reg_coef = obj.pr_reg_coef;
   		activTol = obj.activTol;
        bwd = obj.bwd;

    	/* */
    	if (zRed != NULL)
    		delete [] zRed;
    	if (obj.zRed != NULL){
    		zRed_ = new double[dZ];
    		std::copy(obj.zRed,obj.zRed+obj.dZ,zRed_);
    		zRed = zRed_;
    	}
    	else {
    		zRed = NULL;
    	}

    	/* */
    	if (rRed != NULL)
    		delete [] rRed;
    	if (obj.rRed != NULL){
    		rRed_ = new double[dZ];
    		std::copy(obj.rRed,obj.rRed+obj.dZ,rRed_);
    		rRed = rRed_;
    	}
    	else {
    		rRed = NULL;
    	}

    	/* */
    	if (yRed != NULL)
    		delete [] yRed;
    	if (obj.yRed != NULL){
    		yRed_ = new double[dZ]; 
    		std::copy(obj.yRed,obj.yRed+obj.dZ,yRed_);
    		yRed = yRed_;
    	}
    	else {
    		yRed = NULL;
    	}

    	/* */
    	if (pRed != NULL)
    		delete [] pRed;
    	if (obj.pRed != NULL){
    		pRed_ = new double[dZ];
    		std::copy(obj.pRed,obj.pRed+obj.dZ,pRed_);
    		pRed = pRed_;
    	}
    	else {
    		pRed = NULL;
    	}

    	/* */
    	if (xRed != NULL)
    		delete [] xRed;
    	if (obj.xRed != NULL){
    		xRed_ = new double[dZ]; 
    		std::copy(obj.xRed,obj.xRed+obj.dZ,xRed_);
    		xRed = xRed_;
   	 	} 	
   	 	else {
   	 		xRed = NULL;
   	 	}	

    	/* */
    	if (zRedmin != NULL)
    		delete [] zRedmin;
    	if (obj.zRedmin != NULL){
    		zRedmin_ = new double[dZ];
    		std::copy(obj.zRedmin,obj.zRedmin+obj.dZ,zRedmin_);
    		zRedmin = zRedmin_;
    	}
    	else {
    		zRedmin = NULL;
    	}
    
    	/* */
    	if (zRedmax != NULL)
    		delete [] zRedmax;
    	if (obj.zRedmax != NULL){
    		zRedmax_ = new double[dZ];
    		std::copy(obj.zRedmax,obj.zRedmax+obj.dZ,zRedmax_);	
    		zRedmax = zRedmax_;
    	}
    	else {
    		zRedmax = NULL;
    	}
	}
	return *this;
}

/**
* Initialization method
*/
void PreconRefine::init(){

    assert(dZ>=1);

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
    
    /* zRedmin */
    if (zRedmin != NULL)
        delete [] zRedmin;
    zRedmin = new double[dZ];

    /* zRedmax */
    if (zRedmax != NULL)
        delete [] zRedmax;
    zRedmax = new double[dZ];

    /* Init. preconditioner */
    precon.setDimFull(dZ);
    precon.setReducedHessians(hALshootRed,hALpathRed);
    precon.setDimFull(dZ);
    precon.setNs_(Ns_);
    precon.setShootDim(dZshoot,dxux,dhS,dHshoot,dGshoot);
    precon.setSlackDim(dZpath,dpc,dhP,dHpath,dGpath);
    precon.setDiagRegCoef(pr_reg_coef);
    precon.setFreeVars(freeVars);
    precon.setVarStatus(varStatus);
    precon.setIndexFreeVars(ivarFreeState,ivarFreeInput,ivarFreeSlack,ivarFreeGroupShoot);
    precon.setNumberFreeVars(nfreeState,nfreeInput,nfreeSlack,nfreeGroupShoot);
    precon.init();
}

/**
* Assemble reduced hessian Z'BZ in a separated blockwise manner
*/
void PreconRefine::assembleReducedHessian(){

    int i,j,n,nfree,irw;
    int *nfgS,*ivfgS;
    double *hS,*hSr;
#if DPC
    int *nfP,*ivfP;
    double *hP,*hPr;

    hP = hALpath;
    hPr = hALpathRed;
    ivfP = ivarFreeSlack;
    nfP = nfreeSlack;
#endif

    hS = hALshoot;
    hSr = hALshootRed;
    ivfgS = ivarFreeGroupShoot;
    nfgS = nfreeGroupShoot;
    for (n=0;n<Ns_;++n){    
        /* Shooting group */
        nfree = *(nfgS+n);
        for (i=0;i<nfree;++i){
            irw = *(ivfgS+i)*dxux;
            for (j=0;j<nfree;++j)
                *hSr++ = *(hS+irw+*(ivfgS+j));
        }
        hS += dhS;
        ivfgS += dxux;
#if DPC
        /* Slack group */
        nfree = *(nfP+n);
        for (i=0;i<nfree;++i){
            irw = *(ivfP+i)*dpc;
            for (j=0;j<nfree;++j)
                *hPr++ = *(hP+irw+*(ivfP+j));
        }   
        hP += dhP;
        ivfP += nfree;
#endif
    }
}

/**
* Extract active-set at candidate point zNex and compute residual -Zm'*(g+B(z)*(zNex-z))
*/
void PreconRefine::extractActiveSet(){

    int k,ngcur,ngpre,nx,idx,nu,nfr,iu;
    int *iifv;
    double azn,al,au; 
    double *zznS,*zzlS,*zzuS;
    double *rr;
#if DPC
    int ns;
    int *ivfP;
    double *zznP,*zzlP,*zzuP
#endif

    /* Set counters of # free variables to zero */
    memset(nfreeState,0,NsNT);
    memset(nfreeInput,0,Ns_NT);
    memset(nfreeGroupShoot,0,Ns_NT);
#if DPC
    memset(nfreeSlack,0,Ns_NT);
#endif

    rr = rRed;
    /** Shooting loop */
    dimRedShoot = 0;
    iifv = indFreeVars;
    zznS = zNex;
    zzlS = zmin;
    zzuS = zmax;
    for (k=0;k<dZshoot;++k){
        azn = *zznS++;
        al = *zzlS++;  
        au = *zzuS++;
        /*** Lower bound active */
        if (ABS(azn-al)<activTol){
            *(freeVars+k) = 0;
            *(varStatus+k) = -1;            
            continue;
        }
        /*** Upper bound active */
        if (ABS(azn-au)<activTol){
            *(freeVars+k) = 0;
            *(varStatus+k) = 1;
            continue;
        }
        /* No active bound, free variable */
        *rr++ = -*(gAL+k)-*(prdHss+k);
        /* Update variable status */
        *(freeVars+k) = 1;
        *(varStatus+k) = 0;
        *iifv++ = k;
        ++dimRedShoot;
        /* 
        * Map free variable to state, input & group 
        */
        /* State number */
        nx = k/dxu; 
        /* Input number */
        nu = nx;
        /* Current group number */
        ngcur = nx;
        /* Previous group number */
        ngpre = nx-1;
        /* Compute index in current state or input */
        idx = k%dxu;
        switch (idx/dx){
            case 0:
                /* Index corresponds to state variable */
                nfreeState[nx]++;
//                *ivfX++ = idx;
                if (ngpre>=0){
                    nfr = nfreeGroupShoot[ngpre];
                    *(ivarFreeGroupShoot+ngpre*dxux+nfr) = dxu+idx;
                    nfreeGroupShoot[ngpre]++;
                }
                if (ngcur<=Ns_-1){
                    nfr = nfreeGroupShoot[ngcur];
                    *(ivarFreeGroupShoot+ngcur*dxux+nfr) = idx;
                    nfreeGroupShoot[ngcur]++;
                }
                break;
            case 1:
                /* Index corresponds to input variable */
                iu = idx%dx;
                nfreeInput[nu]++;
//                *ivfU++ = iu;
                nfr = nfreeGroupShoot[ngcur];
                *(ivarFreeGroupShoot+ngcur*dxux+nfr) = dx+iu;
                nfreeGroupShoot[ngcur]++;
                break;
            default:
//                std::cerr<<"ERROR<ActivityDetector>: problem in extracting active-set from shooting block, idx/dx= "<<idx/dx<<", idx= "<<idx<<", dx= "<<dx<<std::endl;
                return;
        }
    }
    /** 
    * Slacks loop
    */
    dimRedPath = 0;
#if DPC
    ivfP = ivarFreeSlack;
    zznP = zNex+dZshoot;
    zzlP = zzlS+dZshoot;
    zzuP = zzuS+dZshoot;
    for (k=dZshoot;k<dZ;++k){
        azn = *zznP++;
        al = *zzlP++;
        au = *zzuP++;
        /*** Lower bound active */
        if (ABS(azn-al)<activTol){
            *(varStatus+k) = -1;
            *(freeVars+k) = 0;
            continue;
        }
        /*** Upper bound active */
        if (ABS(azn-au)<activTol){
            *(varStatus+k) = 1;
            *(freeVars+k) = 0;
            continue;
        }
        /* No active bound, free variable */
        *rr++ = -*(gAL+k)-*(prdHss+k);
        /* Update variable status */
        *(varStatus+k) = 0;
        *(freeVars+k) = 1;
        *iifv++ = k;
        ++dimRedPath;
        /* Map variable to slack group */
        /* Slack number */
        ns = (k-dZshoot)/dpc;
        /* Compute index in current slack */
        nfreeSlack[ns]++;
        *ivfP++ = (k-dZshoot)%dpc;
    }
#endif
    /* Gather total number of free variables */
    dimRed = dimRedShoot+dimRedPath;
    dimRedBL = dimRed*SZDBL;
}

/**
* Find candidate by optimizing over null-space of active constraints 
* with restart when problem bound is hit during search.
*/
void PreconRefine::findCandidateRestart(const double tol2,const double trustRadius){

	int it,n,i,j,k,idf;
    int iu,nx,nu,ns;
    int nfr,nfrg,nfrx,nfru;
    int ngcur,ngpre,idx;
    int nwlb,nwub;
    int *vidf,*iifv,*ivfX,*ivfU,*ivfP;
    int *nfX,*nfU,*nfgS;
    double aux,auxr,auxy;
    double curv,px,prod;
    double rh1,rh2,n2res,beta;
    double a1,a2,ax,al,au;
    double zzr;
    double *zrm,*zrM;
    double *zmm,*zMM;
    double *zzS,*znS;
    double *zr,*rr,*yr,*pr,*xr;
    double *xrS,*xxrS;
    double *prS,*pprS;
    double *hs,*hsr;
    double *phsS,*pphsS;
#if DPC
    int nfrp;
    int *nfP;
    double *prP,*pprP;
    double *hpr,*xrP;
    double *phsP;
    double *hp;
    double *zzP,*znP;
#endif

    /* Extract active-set */
    extractActiveSet();
    maxIt = dimRed;
    /* Assemble reduced hessian */
    assembleReducedHessian();

    /* Build preconditioner and apply it to reduced residual by solving M.yRed = rRed */
    precon.setDimRed(dimRed,dimRedShoot,dimRedPath);
#if JACO_PRE 
    /* Jacobi preconditioner */
    precon.buildJacobi();
    precon.applyJacobi(yRed,rRed);
#endif
#if BAND_PRE
    /* Banded preconditioner */
    precon.buildBand();
    precon.applyBand(yRed,rRed);
#endif

    /* Init. optimizer and reduced bounds
    *  Compute squared 2-norm of model residual n2res and preconditioned residual rh2 
    */
    rh1 = 1.;
    rh2 = 0.;
    n2res = 0.;
    rr = rRed; // Residual 
    yr = yRed; // Preconditioned residual
    zr = zRed;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
    zrm = zRedmin;
    zrM = zRedmax;
    vidf = indFreeVars;
    while (vidf != indFreeVars+dimRed){
        idf = *vidf++;
        *zrm++ = *(zmin+idf);
        *zrM++ = *(zmax+idf);
        *zr++ = *(zNex+idf);
        auxr = *rr++;
        rh2 += *yr++*auxr;
        n2res += auxr*auxr;
    }
    /* Init. CG direction */
    memset(pRed,0,dimRedBL);
	/**
	** Safeguarded TPCG main loop
	**/
	it = 0;
#if COUNT_CG
    cumCG = 0;
    cumRestart = 0;
#endif
	while (1){
		/* Check tolerance and # iterations, stop if satisfied */
        if ((n2res<=tol2)||(it>=maxIt)){
//        if ((n2res<=tol2)||(it>=100)){
            zr = zRed;
            vidf = indFreeVars;
            while (vidf != indFreeVars+dimRed){
                idf = *vidf++;
                *(zNex+idf) = *zr++;
            }
#if PRINT_CG            
            std::cout<<"SUBSPACE-REFINER: tolerance or max. num. iters hit."<<std::endl;
            std::cout<<"n2res="<<n2res<<", tol2="<<tol2<<std::endl;
            std::cout<<"it="<<it<<", maxIt="<<maxIt<<std::endl;
#endif
            break;
        }
        
        /* Step-size computation */
        beta = rh2/rh1;
        
        /* Conjugate direction update pr <- yr + beta*pr */
        pr = pRed;
        yr = yRed;
        while (pr != pRed+dimRed){
            aux = *pr;
            *pr++ = beta*aux+*yr++;
        }
        
        /**
        * Product with reduced AL hessian xRed <- Hred*pRed
        */
        memset(xRed,0,dimRedBL);

        hsr = hALshootRed;
        nfX = nfreeState;
        nfU = nfreeInput;
        nfgS = nfreeGroupShoot;
        prS = pRed;
        xrS = xRed;
#if DPC 
        hpr = hALpathRed;       
        nfP = nfreeSlack;
        prP = pRed+dimRedShoot;
        xrP = xRed+dimRedShoot;
#endif
        for (n=0;n<Ns_;++n){
            /* Reduced shooting block */
            nfrg = *nfgS++;
            nfrx = *nfX++;
            nfru = *nfU++;
            xxrS = xrS;
            for (i=0;i<nfrg;++i){
                aux = 0.;
                pprS = prS;
                for (j=0;j<nfrg;++j)
                    aux += *hsr++**pprS++;
                *xxrS++ += aux; 
            }
            xrS += nfrx+nfru;
            prS += nfrx+nfru;
#if DPC
            /* Reduced slack block */
            nfrp = *nfP++;
            for (i=0;i<nfrp;++i){
                aux = 0.;
                pprP = prP;
                for (j=0;j<nfrp;++j)
                    aux += *hpr++**pprP++;
                *xrP++ = aux;
            }   
            prP += nfrp;
#endif
        }

        /**  1) Compute curvature
             2) Find largest step-size such that (trust-region or NLP) bound is hit */
        a1 = POSINF;
        curv = 0.;
        zrm = zRedmin;
        zrM = zRedmax;
        pr = pRed;
        zr = zRed;
        xr = xRed;
        while (xr != xRed+dimRed){
            px = *pr;
            curv += *xr*px;
            zzr = *zr;
            
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
            vidf = indFreeVars;
            while (vidf != indFreeVars+dimRed){
                idf = *vidf++;
                *(zNex+idf) = *zr+*pr*a1;
                ++zr;
                ++pr;
            }
#if PRINT_CG
            std::cout<<"SUBSPACE-REFINER: nonpositive curvature after "<<it<<" iters, curv="<<curv<<std::endl;
            std::cout<<"n2res= "<<n2res<<", tol2= "<<tol2<<std::endl;
#endif           
            break;
        }
    
        a2 = rh2/curv;

        /** If CG step-size too large, constraint hit:
            - case 1: trust-region bound => stop
            - case 2: NLP bounds => restart on updated subspace **/
        if (a2>a1){
            /** Update zNex <- zRed+a1*pRed and extract newly active NLP bounds 
                (don't consider trust-region bounds)*/ 
            nwlb = 0;
            nwub = 0;
            zr = zRed;
            pr = pRed;
            vidf = indFreeVars;
            while (vidf != indFreeVars+dimRed){
                idf = *vidf++;
                ax = *pr++*a1+*zr++;
                al = *(zmin+idf);
                au = *(zmax+idf);
                /* NLP lower bound becomes active */
                if (ABS(ax-al)<activTol){
                    ++nwlb;
                    *(freeVars+idf) = 0;
                    *(varStatus+idf) = -1;
                    *(zNex+idf) = al;
                    continue;
                }
                /* NLP upper bound becomes active */
                if (ABS(ax-au)<activTol){
                    ++nwub; 
                    *(freeVars+idf) = 0;
                    *(varStatus+idf) = 1;
                    *(zNex+idf) = au;
                    continue;
                }
                /* Variable remains free */
                *(zNex+idf) = ax;
            }

            /* If no NLP bound hit, abort (we've hit trust-region bound) */
            if ((nwlb==0)&&(nwub==0)){
#if PRINT_CG
                std::cout<<"SUBSPACE-REFINER: TR constraint hit after "<<it<<" iters."<<std::endl;
                std::cout<<"n2res= "<<n2res<<", tol2= "<<tol2<<std::endl;
#endif
                break;
            }

            /** 
            * If NLP bound has become active, prepare for restart on new subspace 
            * with updated freeVars, varStatus and zNex 
            */
            /* Recompute prdHss = H*(zNex-z) */
            memset(prdHss,0,dZBL); 

            phsS = prdHss;
            hs = hALshoot;
            zzS = z;
            znS = zNex;
#if DPC
            phsP = prdHss+dZshoot;
            hp = hALpath;
            zzP = z+dZshoot;
            znP = zNex+dZshoot;
#endif
            for (n=0;n<Ns_;++n){
                /* Full shooting block */
                pphsS = phsS;
                for (i=0;i<dxux;++i){
                    prod = 0.;
                    for (j=0;j<dxux;++j)
                        prod += *hs++*(*(znS+j)-*(zzS+j));
                    *pphsS++ += prod;
                }   
                phsS += dxu;
                zzS += dxu;
                znS += dxu;
#if DPC
                /* Full slack block */
                for (i=0;i<dpc;++i){
                    prod = 0.;
                    for (j=0;j<dpc;++j)
                        prod += *hp++*(*(znP+j)-*(zzP+j));
                    *phsP++ = prod;
                }
                zzP += dpc;
                znP += dpc;
#endif
            }

            /** Update rRed, indices of free variables in state, input and group */
            /** Init. counters of free state, input and group variables to 0 */ 
            memset(nfreeState,0,NsNT);
            memset(nfreeInput,0,Ns_NT);
            memset(nfreeGroupShoot,0,Ns_NT);
#if DPC
            memset(nfreeSlack,0,Ns_NT);
#endif

            /** Shooting loop */
            dimRedShoot = 0;
            iifv = indFreeVars;
//            ivfX = ivarFreeState;
//            ivfU = ivarFreeInput;
            rr = rRed;  
            for (k=0;k<dZshoot;++k){
                /* Map free variable to state, input & group */
                if (*(varStatus+k)==0){
                    *rr++ = -*(gAL+k)-*(prdHss+k);
                    *iifv++ = k;
                    ++dimRedShoot;

                    /* State number */
                    nx = k/dxu; 
                    /* Input number */
                    nu = nx;
                    /* Current group number */
                    ngcur = nx;
                    /* Previous group number */
                    ngpre = nx-1;
                    /* Compute index in current state or input */
                    idx = k%dxu;
                    switch (idx/dx){
                        case 0:
                            /* Index corresponds to state variable */
                            nfreeState[nx]++;
//                            *ivfX++ = idx;
                            if (ngpre>=0){
                                nfr = nfreeGroupShoot[ngpre];
                                *(ivarFreeGroupShoot+ngpre*dxux+nfr) = dxu+idx;
                                nfreeGroupShoot[ngpre]++;
                            }
                            if (ngcur<=Ns_-1){
                                nfr = nfreeGroupShoot[ngcur];
                                *(ivarFreeGroupShoot+ngcur*dxux+nfr) = idx;
                                nfreeGroupShoot[ngcur]++;
                            }
                            break;
                        case 1:
                            /* Index corresponds to input variable */
                            iu = idx%dx;
                            nfreeInput[nu]++;
//                            *ivfU++ = iu;
                            nfr = nfreeGroupShoot[ngcur];
                            *(ivarFreeGroupShoot+ngcur*dxux+nfr) = dx+iu;
                            nfreeGroupShoot[ngcur]++;
                            break;
                        default:
                        //    std::cerr<<"ERROR<ActivityDetector>: problem in extracting active-set from shooting block, aux/dx= "<<aux/dx<<std::endl;
                            return;
                    }

                }
            }
            /** Slacks loop */
            dimRedPath = 0;
#if DPC
            ivfP = ivarFreeSlack;
            for (k=dZshoot;k<dZ;++k){
                /* Map free variable to slack group */
                if (*(varStatus+k)==0){
                    *rr++ = -*(gAL+k)-*(prdHss+k);
                    *iifv++ = k;
                    ++dimRedPath;
                    /* Slack number */
                    ns = (k-dZshoot)/dpc;
                    /* Compute index in current slack */
                    nfreeSlack[ns]++;
                    *ivfP++ = (k-dZshoot)%dpc;
                }
            }
#endif
            dimRed = dimRedShoot+dimRedPath;
            dimRedBL = dimRed*SZDBL;

            /* Build preconditioner from updated reduced hessian */
            precon.setDimRed(dimRed,dimRedShoot,dimRedPath);
            assembleReducedHessian();
#if JACO_PRE
            precon.buildJacobi();
            precon.applyJacobi(yRed,rRed);
#endif
#if BAND_PRE
            precon.buildBand();
            precon.applyBand(yRed,rRed);
#endif
            /* Reinit. residuals and bounds */
            rh1 = 1.;
            rh2 = 0.;
            n2res = 0.;
            rr = rRed; 
            yr = yRed;
            zr = zRed;
            zrm = zRedmin;
            zrM = zRedmax;
            vidf = indFreeVars;
            /* This loop could be done inside the applyJacobi() method !!! */
            while (vidf != indFreeVars+dimRed){
                idf = *vidf++;
                *zrm++ = *(zmin+idf);
                *zrM++ = *(zmax+idf);
                *zr++ = *(zNex+idf);
                auxr = *rr++;
                rh2 += *yr++*auxr;
                n2res += auxr*auxr;
            }
            /* Init. CG direction */
            memset(pRed,0,dimRedBL);
            /* Restart */
#if COUNT_CG
            ++cumRestart;
#endif
#if PRINT_CG
            std::cout<<"Restart after "<<it<<" TPCG iters."<<std::endl;
            std::cout<<"!!!! PRECON-REFINE, RESTART !!!!"<<std::endl;
#endif
            it = 0;
            continue;
        }
        /** 
        * PCG update 
        */
        zr = zRed;
        pr = pRed;
        xr = xRed;
        rr = rRed;
        while (zr != zRed+dimRed){
            *zr += *pr*a2;
            *rr -= *xr*a2;
            ++zr;
            ++xr;
            ++pr;
            ++rr;
        }
        /* Apply preconditioner, ie solve M.yRed = rRed */
#if JACO_PRE
        precon.applyJacobi(yRed,rRed);
#endif
#if BAND_PRE        
        precon.applyBand(yRed,rRed);
#endif
        /* Update residuals */
        rh1 = rh2;

        rr = rRed;
        yr = yRed;
        rh2 = 0.;
        n2res = 0.;
        /* This loop could be done inside the applyJacobi() or applyBand() methods !!! */
        while (rr != rRed+dimRed){
            auxr = *rr++;
            auxy = *yr++;
            rh2 += auxy*auxr;
            n2res += auxr*auxr;
        }
        /* Go to next PCG iteration */
		++it;
#if COUNT_CG
        ++cumCG;
#endif
	}
}

/**
* Real-time version of subspace refinement (abort if cumulative PCG larger than max)
* Find candidate by optimizing over null-space of active constraints 
* with restart when problem bound is hit during search.
*/
void PreconRefine::findCandidateRestartRT(const double tol2,const double trustRadius,const int cumCGprim){

    int it,n,i,j,k,idf;
    int iu,nx,nu,ns;
    int nfr,nfrg,nfrx,nfru;
    int ngcur,ngpre,idx;
    int nwlb,nwub;
    int *vidf,*iifv,*ivfX,*ivfU,*ivfP;
    int *nfX,*nfU,*nfgS;
    double aux,auxr,auxy;
    double curv,px,prod;
    double rh1,rh2,n2res,beta;
    double a1,a2,ax,al,au;
    double zzr;
    double *zrm,*zrM;
    double *zmm,*zMM;
    double *zzS,*znS;
    double *zr,*rr,*yr,*pr,*xr;
    double *xrS,*xxrS;
    double *prS,*pprS;
    double *hs,*hsr;
    double *phsS,*pphsS;
#if DPC
    int nfrp;
    int *nfP;
    double *prP,*pprP;
    double *hpr,*xrP;
    double *phsP;
    double *hp;
    double *zzP,*znP;
#endif

    /* Extract active-set */
    extractActiveSet();
    maxIt = dimRed;
    /* Assemble reduced hessian */
    assembleReducedHessian();

    /* Build preconditioner and apply it to reduced residual by solving M.yRed = rRed */
    precon.setDimRed(dimRed,dimRedShoot,dimRedPath);
#if JACO_PRE 
    /* Jacobi preconditioner */
    precon.buildJacobi();
    precon.applyJacobi(yRed,rRed);
#endif
#if BAND_PRE
    /* Banded preconditioner */
    precon.buildBand();
    precon.applyBand(yRed,rRed);
#endif

    /* Init. optimizer and reduced bounds
    *  Compute squared 2-norm of model residual n2res and preconditioned residual rh2 
    */
    rh1 = 1.;
    rh2 = 0.;
    n2res = 0.;
    rr = rRed; // Residual 
    yr = yRed; // Preconditioned residual
    zr = zRed;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
    zrm = zRedmin;
    zrM = zRedmax;
    vidf = indFreeVars;
    while (vidf != indFreeVars+dimRed){
        idf = *vidf++;
        *zrm++ = *(zmin+idf);
        *zrM++ = *(zmax+idf);
        *zr++ = *(zNex+idf);
        auxr = *rr++;
        rh2 += *yr++*auxr;
        n2res += auxr*auxr;
    }
    /* Init. CG direction */
    memset(pRed,0,dimRedBL);
    /**
    ** Safeguarded TPCG main loop
    **/
    it = 0;
#if COUNT_CG
    cumCG = 0;
    cumRestart = 0;
#endif
    while (1){
        /* Check tolerance and # iterations, stop if satisfied */
        /* For RT: abort if total cumulative PCG count above max */
        if ((n2res<=tol2)||(it>=maxIt)||((cumCGprim+cumCG)>maxCumCG)){
//        if ((n2res<=tol2)||(it>=100)){
            zr = zRed;
            vidf = indFreeVars;
            while (vidf != indFreeVars+dimRed){
                idf = *vidf++;
                *(zNex+idf) = *zr++;
            }
#if PRINT_CG            
            std::cout<<"SUBSPACE-REFINER: tolerance or max. num. iters hit."<<std::endl;
            std::cout<<"n2res="<<n2res<<", tol2="<<tol2<<std::endl;
            std::cout<<"it="<<it<<", maxIt="<<maxIt<<std::endl;
#endif
            break;
        }
        
        /* Step-size computation */
        beta = rh2/rh1;
        
        /* Conjugate direction update pr <- yr + beta*pr */
        pr = pRed;
        yr = yRed;
        while (pr != pRed+dimRed){
            aux = *pr;
            *pr++ = beta*aux+*yr++;
        }
        
        /**
        * Product with reduced AL hessian xRed <- Hred*pRed
        */
        memset(xRed,0,dimRedBL);

        hsr = hALshootRed;
        nfX = nfreeState;
        nfU = nfreeInput;
        nfgS = nfreeGroupShoot;
        prS = pRed;
        xrS = xRed;
#if DPC 
        hpr = hALpathRed;       
        nfP = nfreeSlack;
        prP = pRed+dimRedShoot;
        xrP = xRed+dimRedShoot;
#endif
        for (n=0;n<Ns_;++n){
            /* Reduced shooting block */
            nfrg = *nfgS++;
            nfrx = *nfX++;
            nfru = *nfU++;
            xxrS = xrS;
            for (i=0;i<nfrg;++i){
                aux = 0.;
                pprS = prS;
                for (j=0;j<nfrg;++j)
                    aux += *hsr++**pprS++;
                *xxrS++ += aux; 
            }
            xrS += nfrx+nfru;
            prS += nfrx+nfru;
#if DPC
            /* Reduced slack block */
            nfrp = *nfP++;
            for (i=0;i<nfrp;++i){
                aux = 0.;
                pprP = prP;
                for (j=0;j<nfrp;++j)
                    aux += *hpr++**pprP++;
                *xrP++ = aux;
            }   
            prP += nfrp;
#endif
        }

        /**  1) Compute curvature
             2) Find largest step-size such that (trust-region or NLP) bound is hit */
        a1 = POSINF;
        curv = 0.;
        zrm = zRedmin;
        zrM = zRedmax;
        pr = pRed;
        zr = zRed;
        xr = xRed;
        while (xr != xRed+dimRed){
            px = *pr;
            curv += *xr*px;
            zzr = *zr;
            
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
            vidf = indFreeVars;
            while (vidf != indFreeVars+dimRed){
                idf = *vidf++;
                *(zNex+idf) = *zr+*pr*a1;
                ++zr;
                ++pr;
            }
#if PRINT_CG
            std::cout<<"SUBSPACE-REFINER: nonpositive curvature after "<<it<<" iters, curv="<<curv<<std::endl;
            std::cout<<"n2res= "<<n2res<<", tol2= "<<tol2<<std::endl;
#endif           
            break;
        }
    
        a2 = rh2/curv;

        /** If CG step-size too large, constraint hit:
            - case 1: trust-region bound => stop
            - case 2: NLP bounds => restart on updated subspace **/
        if (a2>a1){
            /** Update zNex <- zRed+a1*pRed and extract newly active NLP bounds 
                (don't consider trust-region bounds)*/ 
            nwlb = 0;
            nwub = 0;
            zr = zRed;
            pr = pRed;
            vidf = indFreeVars;
            while (vidf != indFreeVars+dimRed){
                idf = *vidf++;
                ax = *pr++*a1+*zr++;
                al = *(zmin+idf);
                au = *(zmax+idf);
                /* NLP lower bound becomes active */
                if (ABS(ax-al)<activTol){
                    ++nwlb;
                    *(freeVars+idf) = 0;
                    *(varStatus+idf) = -1;
                    *(zNex+idf) = al;
                    continue;
                }
                /* NLP upper bound becomes active */
                if (ABS(ax-au)<activTol){
                    ++nwub; 
                    *(freeVars+idf) = 0;
                    *(varStatus+idf) = 1;
                    *(zNex+idf) = au;
                    continue;
                }
                /* Variable remains free */
                *(zNex+idf) = ax;
            }

            /* If no NLP bound hit, abort (we've hit trust-region bound) */
            if ((nwlb==0)&&(nwub==0)){
#if PRINT_CG
                std::cout<<"SUBSPACE-REFINER: TR constraint hit after "<<it<<" iters."<<std::endl;
                std::cout<<"n2res= "<<n2res<<", tol2= "<<tol2<<std::endl;
#endif
                break;
            }

            /** 
            * If NLP bound has become active, prepare for restart on new subspace 
            * with updated freeVars, varStatus and zNex 
            */
            /* Recompute prdHss = H*(zNex-z) */
            memset(prdHss,0,dZBL); 

            phsS = prdHss;
            hs = hALshoot;
            zzS = z;
            znS = zNex;
#if DPC
            phsP = prdHss+dZshoot;
            hp = hALpath;
            zzP = z+dZshoot;
            znP = zNex+dZshoot;
#endif
            for (n=0;n<Ns_;++n){
                /* Full shooting block */
                pphsS = phsS;
                for (i=0;i<dxux;++i){
                    prod = 0.;
                    for (j=0;j<dxux;++j)
                        prod += *hs++*(*(znS+j)-*(zzS+j));
                    *pphsS++ += prod;
                }   
                phsS += dxu;
                zzS += dxu;
                znS += dxu;
#if DPC
                /* Full slack block */
                for (i=0;i<dpc;++i){
                    prod = 0.;
                    for (j=0;j<dpc;++j)
                        prod += *hp++*(*(znP+j)-*(zzP+j));
                    *phsP++ = prod;
                }
                zzP += dpc;
                znP += dpc;
#endif
            }

            /** Update rRed, indices of free variables in state, input and group */
            /** Init. counters of free state, input and group variables to 0 */ 
            memset(nfreeState,0,NsNT);
            memset(nfreeInput,0,Ns_NT);
            memset(nfreeGroupShoot,0,Ns_NT);
#if DPC
            memset(nfreeSlack,0,Ns_NT);
#endif

            /** Shooting loop */
            dimRedShoot = 0;
            iifv = indFreeVars;
//            ivfX = ivarFreeState;
//            ivfU = ivarFreeInput;
            rr = rRed;  
            for (k=0;k<dZshoot;++k){
                /* Map free variable to state, input & group */
                if (*(varStatus+k)==0){
                    *rr++ = -*(gAL+k)-*(prdHss+k);
                    *iifv++ = k;
                    ++dimRedShoot;

                    /* State number */
                    nx = k/dxu; 
                    /* Input number */
                    nu = nx;
                    /* Current group number */
                    ngcur = nx;
                    /* Previous group number */
                    ngpre = nx-1;
                    /* Compute index in current state or input */
                    idx = k%dxu;
                    switch (idx/dx){
                        case 0:
                            /* Index corresponds to state variable */
                            nfreeState[nx]++;
//                            *ivfX++ = idx;
                            if (ngpre>=0){
                                nfr = nfreeGroupShoot[ngpre];
                                *(ivarFreeGroupShoot+ngpre*dxux+nfr) = dxu+idx;
                                nfreeGroupShoot[ngpre]++;
                            }
                            if (ngcur<=Ns_-1){
                                nfr = nfreeGroupShoot[ngcur];
                                *(ivarFreeGroupShoot+ngcur*dxux+nfr) = idx;
                                nfreeGroupShoot[ngcur]++;
                            }
                            break;
                        case 1:
                            /* Index corresponds to input variable */
                            iu = idx%dx;
                            nfreeInput[nu]++;
//                            *ivfU++ = iu;
                            nfr = nfreeGroupShoot[ngcur];
                            *(ivarFreeGroupShoot+ngcur*dxux+nfr) = dx+iu;
                            nfreeGroupShoot[ngcur]++;
                            break;
                        default:
                        //    std::cerr<<"ERROR<ActivityDetector>: problem in extracting active-set from shooting block, aux/dx= "<<aux/dx<<std::endl;
                            return;
                    }

                }
            }
            /** Slacks loop */
            dimRedPath = 0;
#if DPC
            ivfP = ivarFreeSlack;
            for (k=dZshoot;k<dZ;++k){
                /* Map free variable to slack group */
                if (*(varStatus+k)==0){
                    *rr++ = -*(gAL+k)-*(prdHss+k);
                    *iifv++ = k;
                    ++dimRedPath;
                    /* Slack number */
                    ns = (k-dZshoot)/dpc;
                    /* Compute index in current slack */
                    nfreeSlack[ns]++;
                    *ivfP++ = (k-dZshoot)%dpc;
                }
            }
#endif
            dimRed = dimRedShoot+dimRedPath;
            dimRedBL = dimRed*SZDBL;

            /* Build preconditioner from updated reduced hessian */
            precon.setDimRed(dimRed,dimRedShoot,dimRedPath);
            assembleReducedHessian();
#if JACO_PRE
            precon.buildJacobi();
            precon.applyJacobi(yRed,rRed);
#endif
#if BAND_PRE
            precon.buildBand();
            precon.applyBand(yRed,rRed);
#endif
            /* Reinit. residuals and bounds */
            rh1 = 1.;
            rh2 = 0.;
            n2res = 0.;
            rr = rRed; 
            yr = yRed;
            zr = zRed;
            zrm = zRedmin;
            zrM = zRedmax;
            vidf = indFreeVars;
            /* This loop could be done inside the applyJacobi() method !!! */
            while (vidf != indFreeVars+dimRed){
                idf = *vidf++;
                *zrm++ = *(zmin+idf);
                *zrM++ = *(zmax+idf);
                *zr++ = *(zNex+idf);
                auxr = *rr++;
                rh2 += *yr++*auxr;
                n2res += auxr*auxr;
            }
            /* Init. CG direction */
            memset(pRed,0,dimRedBL);
            /* Restart */
#if COUNT_CG
            ++cumRestart;
#endif
#if PRINT_CG
            std::cout<<"Restart after "<<it<<" TPCG iters."<<std::endl;
            std::cout<<"!!!! PRECON-REFINE, RESTART !!!!"<<std::endl;
#endif
            it = 0;
            continue;
        }
        /** 
        * PCG update 
        */
        zr = zRed;
        pr = pRed;
        xr = xRed;
        rr = rRed;
        while (zr != zRed+dimRed){
            *zr += *pr*a2;
            *rr -= *xr*a2;
            ++zr;
            ++xr;
            ++pr;
            ++rr;
        }
        /* Apply preconditioner, ie solve M.yRed = rRed */
#if JACO_PRE
        precon.applyJacobi(yRed,rRed);
#endif
#if BAND_PRE        
        precon.applyBand(yRed,rRed);
#endif
        /* Update residuals */
        rh1 = rh2;

        rr = rRed;
        yr = yRed;
        rh2 = 0.;
        n2res = 0.;
        /* This loop could be done inside the applyJacobi() or applyBand() methods !!! */
        while (rr != rRed+dimRed){
            auxr = *rr++;
            auxy = *yr++;
            rh2 += auxy*auxr;
            n2res += auxr*auxr;
        }
        /* Go to next PCG iteration */
        ++it;
#if COUNT_CG
        ++cumCG;
#endif
    }
}


/**
** For debugging
**/ 
/* Displays shooting part of AL hessian */
void PreconRefine::displayHalShoot(){

    int n,i;
    double *hal;

    std::cout<<"PRECOND-REFINE<hALshoot>: "<<std::endl;    
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

/* Displays residuals */
void PreconRefine::displayResiduals(){

    double *rr,*yy;

    rr = rRed;
    yy = yRed;
/*    std::cout<<"Reduced residuals: "<<std::endl;
    while (rr != rRed+dimRed)
        std::cout<<"rRed= "<<*rr++<<", yRed= "<<*yy++<<std::endl; */

    std::cout<<"rRed= "<<std::endl;
    while (rr != rRed+dimRed)
        std::cout<<*rr++<<std::endl;
    std::cout<<"yRed= "<<std::endl;
    while (yy != yRed+dimRed)
        std::cout<<*yy++<<std::endl;
}

/* Displays conjugate direction */
void PreconRefine::displayConjugateDirection(){

    double *pp;

    pp = pRed;
    std::cout<<"Conjugate direction: "<<std::endl;
    while (pp != pRed+dimRed)
        std::cout<<"pRed= "<<*pp++<<std::endl;
}

/* Displays reduced optimizer */
void PreconRefine::displayReducedOptimizer(){

    double *zz;

    zz = zRed;
    std::cout<<"Reduced optimizer: "<<std::endl;
    while (zz != zRed+dimRed)
        std::cout<<"zRed= "<<*zz++<<std::endl;
}

/* Displays reduced product */
void PreconRefine::displayReducedProduct(){

    double *xx;

    xx = xRed;  
    std::cout<<"Reduced product: "<<std::endl;
    while (xx != xRed+dimRed)
        std::cout<<"xRed= "<<*xx++<<std::endl;
}
