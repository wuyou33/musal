/*********
* Implementation of activity detection mechanism for bound-constrained quadratic model
**************/

#include <iostream>
#include <string.h>

#include "ActivityDetector.h"
#include "userData.h"

/** Default constructor */
ActivityDetector::ActivityDetector(){

    zL = NULL;
    zU = NULL;
    
    zC = NULL;
    z = NULL;
    zdbg = NULL;

    grad = NULL;

    nfree = -1;
    nfreeS = -1;
    nfreeP = -1;
    err = 0;

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

    rRed = NULL;
    prdHss = NULL;
    prdHssDbg = NULL;
}

/** Copy constructor */
ActivityDetector::ActivityDetector(const ActivityDetector& obj) : qmod(obj.qmod){

    dx = obj.dx;
    dx_ = obj.dx_;
    dpc = obj.dpc;
    du = obj.du;
    dxu = obj.dxu;
    dxu_ = obj.dxu_;
    dxux = obj.dxux;
    Ns_ = obj.Ns_;
    dZ = obj.dZ;
    dZshoot = obj.dZshoot;
    dGshoot = obj.dGshoot;
    
    zL = obj.zL;
    zU = obj.zU;
    
    zC = obj.zC;
    z = obj.z;

    grad = obj.grad;
    
    err = obj.err;
    nfree = obj.nfree;
    nfreeS = obj.nfreeS;
    nfreeP = obj.nfreeP;

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

    rRed = obj.rRed;
    prdHss = obj.prdHss;

    betaBck = obj.betaBck;
    cstBck = obj.cstBck;
    maxItBck = obj.maxItBck;

    betaInter = obj.betaInter;
    betaExtra = obj.betaExtra;
    maxInterIter = obj.maxInterIter;
    maxExtraIter = obj.maxExtraIter;

    activTol = obj.activTol;
}

/** Assignment operator */
ActivityDetector& ActivityDetector::operator= (const ActivityDetector& obj){

    if (this != &obj){

        qmod = obj.qmod;

        dx = obj.dx;
        dx_ = obj.dx_;
        dpc = obj.dpc;
        du = obj.du;
        dxu = obj.dxu;
        dxu_ = obj.dxu_;
        dxux = obj.dxux;
        Ns_ = obj.Ns_;
        dZ = obj.dZ;
        dZshoot = obj.dZshoot;
        dGshoot = obj.dGshoot;
    
        zL = obj.zL;
        zU = obj.zU;
    
        zC = obj.zC;
        z = obj.z;

        grad = obj.grad;
    
        err = obj.err;
        nfree = obj.nfree;
        nfreeS = obj.nfreeS;
        nfreeP = obj.nfreeP;

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

        rRed = obj.rRed;
        prdHss = obj.prdHss;

        betaBck = obj.betaBck;
        cstBck = obj.cstBck;
        maxItBck = obj.maxItBck;

        betaInter = obj.betaInter;
        betaExtra = obj.betaExtra;
        maxInterIter = obj.maxInterIter;
        maxExtraIter = obj.maxExtraIter;

        activTol = obj.activTol;
    }

    return *this;
}

/** Destructor */
ActivityDetector::~ActivityDetector(){

    if (zdbg != NULL)
        delete [] zdbg;

    if (prdHssDbg != NULL)
        delete [] prdHssDbg;
}

/** Allocate extra memory for debugging */
void ActivityDetector::allocateDebug(){

    zdbg = new double[dZ];
    prdHssDbg = new double[dZ];
}

/** Extract active-set at Cauchy point and compute reduced residual & gradient  */
void ActivityDetector::extractActiveSet(double* n2res,double* n2gra){

    int k,ngcur,ngpre,nx,idx,nu,nfr,iu;
    int *iifv;
    double azc,al,au; 
    double nnr,nng;
    double auxr,auxg;
    double *zzcS;
    double *zzlS,*zzuS;
    double *rr;
#if DPC
    int ns;
    int *ivfP;
    double *zzcP,*zzlP,*zzuP
#endif

    /* Set counters of # free variables to zero */
    memset(nfreeState,0,NsNT);
    memset(nfreeInput,0,Ns_NT);
    memset(nfreeGroupShoot,0,Ns_NT);
#if DPC
    memset(nfreeSlack,0,Ns_NT);
#endif

    rr = rRed;
    nnr = 0.;
    nng = 0.;
    /** 
    * Shooting loop 
    */
    nfreeS = 0;
    iifv = indFreeVars;
//    ivfX = ivarFreeState;
//    ivfU = ivarFreeInput;  
    zzcS = zC;
    zzlS = zL;
    zzuS = zU;
    for (k=0;k<dZshoot;++k){
        azc = *zzcS++;
        al = *zzlS++;  
        au = *zzuS++;
        /*** Lower bound active */
        if (ABS(azc-al)<activTol){
            *(freeVars+k) = 0;
            *(varStatus+k) = -1;            
            continue;
        }
        /*** Upper bound active */
        if (ABS(azc-au)<activTol){
            *(freeVars+k) = 0;
            *(varStatus+k) = 1;
            continue;
        }
        /*** No active bound, free variable */
        auxg = *(grad+k);
        /* For squared 2-norm of reduced gradient */
        nng += auxg*auxg;
        auxr = auxg+*(prdHss+k);
        /* For squared 2-norm of reduced residual */
        nnr += auxr*auxr;
        /* Update reduced residual */
        *rr++ = -auxr;
        /* Update variable status */
        *(freeVars+k) = 1;
        *(varStatus+k) = 0;
        *iifv++ = k;
        ++nfreeS;
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
    nfreeP = 0;
#if DPC
    ivfP = ivarFreeSlack;
    zzcP = zC+dZshoot;
    zzlP = zL+dZshoot;
    zzuP = zU+dZshoot;
    for (k=dZshoot;k<dZ;++k){
        azc = *zzcP++;
        al = *zzlP++;
        au = *zzuP++;
        /*** Lower bound active */
        if (ABS(azc-al)<activTol){
            *(varStatus+k) = -1;
            *(freeVars+k) = 0;
            continue;
        }
        /*** Upper bound active */
        if (ABS(azc-au)<activTol){
            *(varStatus+k) = 1;
            *(freeVars+k) = 0;
            continue;
        }
        /*** No bound active, free variable */
        auxg = *(grad+k);
        /* For squared 2-norm of reduced gradient */
        nng += auxg*auxg;
        auxr = auxg+*(prdHss+k);
        /* For squared 2-norm of reduced residual */
        nnr += auxr*auxr;
        /* Update reduced residual */
        *rr++ = -auxr;
        /* Update variable status */
        *(varStatus+k) = 0;
        *(freeVars+k) = 1;
        *iifv++ = k;
        ++nfreeP;
        /* Map variable to slack group */
        /* Slack number */
        ns = (k-dZshoot)/dpc;
        /* Compute index in current slack */
        nfreeSlack[ns]++;
        *ivfP++ = (k-dZshoot)%dpc;
    }
#endif
    /* Gather total number of free variables */
    nfree = nfreeS+nfreeP;
    /* Gather squared 2-norm of reduced residual */
    *n2res = nnr;
    /* Gather squared 2-norm of reduced gradient */
    *n2gra = nng;
}

/**
 * 1. Compute Cauchy point via Armijo-backtracking projected search
 * 2. Extract active set and map it to partially separable structure
 * - On output, n2res contains squared 2-norm of residual -Zc'*(gAL+B(z)*(zc-z))
 * - On output, n2g contains squared 2-norm of reduced gradient Zc'*gAL
 * - On output, norIstp contains ||zC-z||_\infty
 * - trustRadius is current trust region radius in NonlinOCPsolver 
 */
void ActivityDetector::findProjectedSearch(double* alpha,double* vmod,double* norIstp,const double trustRadius){
     
    int i,k;
    double rhs,lhs,gg,zz;
    double pax,ax,au,al; 
    double nIs,alphas;
    double *phs;

    /**
    ** 1. Backtracking loop
    **/
    i = 0;
    alphas = *alpha;
    while (1){
        /* - zC <- P(z-ss*g)
           - compute rhs of sufficient decrease condition */
        rhs = 0.;
        nIs = -1.;
        phs = prdHss;
        for (k=0;k<dZ;++k){
            zz = *(z+k);
            gg = *(grad+k);
            ax = zz-alphas*gg;
            al = *(zL+k);
            au = *(zU+k);
            pax = PROJ(ax,al,au);
            rhs += gg*(pax-zz);
            nIs = MAX(nIs,ABS(pax-zz));
            *(zC+k) = pax;
            /* Initialize product B(z)*(zC-z) to 0. */
            *phs++ = 0.;
        }
        /* Evaluate 2nd (quadratic) term of model objective */
        lhs = rhs+qmod.evalQuadTermCauchy();
        rhs *= cstBck;
        /* Evaluate lhs of backtracking condition, model value at Cauchy point */
        //    lhs = qmod.evalModelCauchy();
        
        /* Evaluate rhs of backtracking condition */
    /*    rhs = 0.;
        for (k=0;k<dZ;++k)
            rhs += *(grad+k)*(*(zC+k)-*(z+k));
        rhs *= cstBck; */

        /* Stop if sufficient decrease and containment in trust region */
        if ((nIs<=trustRadius)&&(lhs<rhs)){
//            std::cout<<"BREAK<ActivityDetector>: lhs="<<lhs<<", rhs="<<rhs<<std::endl;
//            std::cout<<"BREAK<ActivityDetector>: nIs="<<nIs<<std::endl;
//            displayInfNormGrad();
//            displayInfNormCauchy();
            break; //search = (lhs>=rhs);
        }

        /* If max # backtracking iters reached, output error and stop */
        if (i>=maxItBck){
#if DETEC_CAUCH_ERR
            std::cerr<<"ERROR<ActivityDetector>: Cauchy point not found after max. # backtracking iterations."<<std::endl;
            std::cout<<"alphas="<<alphas<<std::endl;
            std::cout<<"ERROR<ActivityDetector>: lhs="<<lhs<<", rhs="<<rhs<<std::endl;
            std::cout<<"ERROR<ActivityDetector>: trustRadius="<<trustRadius<<std::endl;
            std::cout<<"ERROR<ActivityDetector>: nIs="<<nIs<<std::endl;
            err = 1;
            displayInfNormGrad();
            displayInfNormCauchy();
#endif
            break;
        }        
        ++i;
        /* Shrink step-size */
        alphas *= betaBck;
    }
//    std::cout<<"Stop search after "<<i<<"iter with lhs="<<lhs<<", rhs= "<<rhs<<std::endl;
    *vmod = lhs;
    *alpha = alphas;
    *norIstp = nIs;
}

/**
* 1. Compute Cauchy point via interpolate & extrapolate procedure
* 2. Extract active set at Cauchy point and map to multiple-shooting structure
* - On output, n2res contains squared 2-norm of residual -Zc'*(gAL+B(z)*(zc-z))
* - On output, n2g contains squared 2-norm of reduced gradient Zc'*gAL
* - On output, norIstp contains ||zC-z||_\infty
* - trustRadius is current trust region radius in NonlinOCPsolver
*/
void ActivityDetector::findInterpExtrap(double* alpha,double* vmod,double* norIstp,const double trustRadius){

    int k;
    int interp,iex,iin,search;
    double alphas,alphass,alphamax,gg;
    double nIs,nnIs,gts,stp;
    double az,azc,azl,azu; 
    double lhs,llhs;
    double *gal;
    double *zz,*zzc;
    double *zl,*zu;

    /** 
    ** - Compute maximum break-point smax on line P(z-s*g)
    ** - Compute projected step P(z-s*g)-z and decide whether to interpolate or extrapolate 
    **/
    alphas = *alpha;
    alphamax = NEGINF;
    gal = grad;
    zz = z;
    zzc = zC;
    zl = zL;
    zu = zU;
    nIs = -1.;
    gts = 0.;
    while (zz != z+dZ){
        gg = *gal++;
        az = *zz++;
        azl = *zl++;
        azu = *zu++;
        azc = PROJ(az-alphas*gg,azl,azu);
        *zzc++ = azc;
        stp = azc-az;
        nIs = MAX(nIs,ABS(stp));
        gts += gg*stp;
        if (gg<0){
            alphamax = MAX(alphamax,(azl-az)/gg);
            continue;
        }
        if (gg>0){
            alphamax = MAX(alphamax,(azu-az)/gg);
            continue;
        }
    }

    /* Decide to interpolate if outside trust region or sufficient decrease condition not satisfied */
    if (nIs>trustRadius){
        interp = 1; 
//        std::cout<<"ActivityDetector: NOT IN TRUST REGION."<<std::endl;
//        std::cout<<"ActivityDetector, nIs="<<nIs<<", trustRadius="<<trustRadius<<std::endl;
    }
    else {
        lhs = qmod.evalModelCauchy();
        interp = (lhs>=cstBck*gts);
    }

    /**
    ** Interpolation or extrapolation procedure
    **/
    err = 0;
    switch (interp){
        case 0:
            /* Extrapolation: increase step-size as long as inside trust region and sufficient decrease satisfied */
//            std::cout<<"EXTRAPOLATE"<<std::endl;
            iex = 0;
            search = 1;
            alphass = alphas;
            nnIs = nIs;
            llhs = lhs;
            while (search && (alphas<=alphamax)){
                /* Increase step-size */
                alphas *= betaExtra;
                /* Compute candidate for Cauchy point P(z-s*g) */
                zz = z;
                zzc = zC;
                zl = zL;
                zu = zU;
                gal = grad;
                nIs = -1.;
                gts = 0.;
                while (zz != z+dZ){
                    gg = *gal++;
                    az = *zz++;
                    azc = PROJ(az-alphas*gg,*zl,*zu);
                    *zzc++ = azc;
                    stp = azc-az;
                    nIs = MAX(nIs,ABS(stp));
                    gts += gg*stp;
                    ++zl;
                    ++zu;
                }
                /* Check containment in trust region and sufficient decrease */
                if (nIs<=trustRadius){
                    lhs = qmod.evalModelCauchy();
                    search = (lhs<cstBck*gts);
                    if (search){
                        alphass = alphas;
                        nnIs = nIs;
                        llhs = lhs;
                    }
                //    sss = search*ss+(1-search)*sss;
                //    nnIs = search*nIs+(1-search)*nnIs;
                }
                else {
                    search = 0;
                }
                ++iex;
            }
 //           std::cout<<"iex="<<iex<<std::endl;
            /* Recover last successful iterate */
            alphas = alphass;
            nIs = nnIs;
            lhs = llhs;
            /* Compute Cauchy point P(z-s*g) at last successful iterate */
            zz = z;
            zzc = zC;
            zl = zL;
            zu = zU;
            gal = grad;
            while (zz != z+dZ){
                gg = *gal++;
                az = *zz++;
                *zzc++ = PROJ(az-alphas*gg,*zl,*zu);
                ++zl;
                ++zu;
            }
            break;
        case 1:
            /* Interpolation: decrease step-size until inside trust region and sufficient decrease satisfied */
//            std::cout<<"INTERPOLATE"<<std::endl;
//            search = 1;
            iin = 0;
            while (1){
                if (iin>=maxInterIter){
                    err = 1;
                    return;
                }
                /* Compute candidate for Cauchy point P(z-s*g) */
                zz = z;
                zzc = zC;
                zl = zL;
                zu = zU;
                gal = grad;
                nIs = -1.;
                gts = 0.;
                while (zz != z+dZ){
                    gg = *gal++;
                    az = *zz++;
                    azc = PROJ(az-alphas*gg,*zl,*zu);
                    *zzc++ = azc;
                    stp = azc-az;
                    nIs = MAX(nIs,ABS(stp));
                    gts += gg*stp;
                    zl++;
                    zu++;
                }
                lhs = qmod.evalModelCauchy();
                /* Check containment in trust region and sufficient decrease */
                if ((nIs<=trustRadius)&&(lhs<cstBck*gts))                    
                    break;
                //    search = (lhs>=cstBck*gts);
                /* Shrink step-size */
                alphas *= betaInter;
                ++iin;
            } 
            break;
        default:
            std::cout<<"ACTIVITY-DETECTOR ERROR: cannot decide whether to interpolate or extrapolate."<<std::endl;
            break;
    }
    /* Save Cauchy step-size */
    *alpha = alphas;
    /* Save infinity norm of Cauchy step */
    *norIstp = nIs;
    /* Save model value at Cauchy point */
    *vmod = lhs;
}

/**
** For debugging
**/
/* Check model value along projected arc */
void ActivityDetector::checkDescent(){

    int ii;
    double ss;
    double *zz,*zzd;
    double *zl,*zu;
    double *gg;

//    std::cout<<"<<<ActivityDetector>>>, check descent:"<<std::endl;
    ss = 1.;
    ii =0;
    while (ii<maxItBck){
        zz = z;
        zzd = zdbg;
        gg = grad;
        zl = zL;
        zu = zU;
        while (zz != z+dZ){
//           std::cout<<"gg="<<*gg<<std::endl;
            *zzd = PROJ(*zz-ss**gg,*zl,*zu);
//            std::cout<<"*zzd="<<*zzd<<std::endl;
//            *zzd = *zz-ss**gg;
            ++zzd;
            ++zz;
            ++gg;
            ++zl;
            ++zu;
        }
        std::cout<<"<<<ActivityDetector>>> model for step = "<<ss<<", "<<qmod.evalModel(zdbg,prdHssDbg)<<std::endl;
//        std::cout<<"<<<ActivityDetector>>>, model at base point: "<<qmod.evalModel(z)<<std::endl;
        ss *= betaBck;
        ++ii;
    }
}

void ActivityDetector::displayCauchy(){

    double *zc;

    zc = zC;
    std::cout<<"Cauchy point:"<<std::endl;
    while (zc != zC+dZ)
        std::cout<<*zc++<<std::endl;
}

void ActivityDetector::displayInfNormCauchy(){

    double nIzc;
    double *zzc,*zz,*zl,*zu;

    nIzc = -1.;
    zzc = zC;
    zz = z;
    zl = zL;
    zu = zU;
    while (zzc != zC+dZ){
//        std::cout<<"*zl="<<*zl++<<", *zz="<<*zz++<<", *zu="<<*zu++<<std::endl;
        if ((*zz>*zu)||(*zz<*zl)){
            std::cout<<"MIIIIIAOOOOOOOOO"<<std::endl;
            std::cout<<"*zl="<<*zl<<", *zz="<<*zz<<", *zu="<<*zu<<std::endl;
        }
        nIzc = MAX(nIzc,ABS(*zzc));
        ++zzc;
        ++zz;
        ++zl;
        ++zu;
    }
    std::cout<<"<<<ActivityDetector>>> Infinity-norm of Cauchy point:"<<nIzc<<std::endl;
}

void ActivityDetector::displayInfNormGrad(){

    double nIg;
    double *gg;

    nIg = -1.;
    gg = grad;
    while (gg != grad+dZ){
        nIg = MAX(nIg,ABS(*gg));
        ++gg;
    }
    std::cout<<"<<<ActivityDetector>>> Infinity-norm of model gradient:"<<nIg<<std::endl;
}

void ActivityDetector::displayGradient(){

    double *gg;

    gg = grad;
    std::cout<<"<<ActivityDetector>> Model gradient="<<std::endl;
    while (gg != grad+dZ)
        std::cout<<*gg++<<std::endl;
    qmod.displayGradient();
}

void ActivityDetector::displayProductHessian(){

    double *pp;

    pp = prdHss;
    std::cout<<"Hessian product:"<<std::endl;
    while (pp != prdHss+dZ)
        std::cout<<*pp++<<std::endl;
}

void ActivityDetector::displayFreeStates(){

    int i,n;

    std::cout<<"Number of free vars per state:"<<std::endl;
    for (n=0;n<Ns_+1;++n)
        std::cout<<nfreeState[n]<<std::endl;

    std::cout<<"Indices of free states:"<<std::endl;
    for (i=0;i<dx*(Ns_+1);++i)
        std::cout<<ivarFreeState[i]<<std::endl;
}

void ActivityDetector::displayFreeInputs(){

    int i,n;

    std::cout<<"Number of free vars per input:"<<std::endl;
    for (n=0;n<Ns_;++n)
        std::cout<<nfreeInput[n]<<std::endl;

    std::cout<<"Indices of free inputs:"<<std::endl;
    for (i=0;i<du*Ns_;++i)
        std::cout<<ivarFreeInput[i]<<std::endl;
}

void ActivityDetector::displayFreeGroups(){

    int i,n;

    std::cout<<"Number of free vars per group:"<<std::endl;
    for (n=0;n<Ns_;++n)
        std::cout<<nfreeGroupShoot[n]<<std::endl;

    std::cout<<"Indices in free groups"<<std::endl;
    for (i=0;i<dxux*Ns_;++i)
        std::cout<<ivarFreeGroupShoot[i]<<std::endl;
}

void ActivityDetector::displayFreeIndices(){

    int i;
    std::cout<<"Indices of free variables: "<<std::endl;
    for (i=0;i<nfree;++i)
        std::cout<<indFreeVars[i]<<std::endl;
}