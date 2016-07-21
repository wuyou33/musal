//
//  QuadraticModel.cpp
//  MusAL
//
//  Created by Jean on 4/14/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//
#include <string.h>

#include "QuadraticModel.h"
#include "userData.h"

/**
 * Default constructor
 */
QuadraticModel::QuadraticModel(){
    
    z = NULL;
    zC = NULL;
    zS = NULL;
    
    grad = NULL;
    
    hessShoot = NULL;
    hessPath = NULL;

    prdHss = NULL;
}

/**
 * Copy constructor
*/
QuadraticModel::QuadraticModel(const QuadraticModel& obj){
    
    dZ = obj.dZ;
    dZBL = obj.dZBL;
    dZshoot = obj.dZshoot;
    
    z = obj.z;
    zC = obj.zC;
    zS = obj.zS;
    
    grad = obj.grad;
    
    hessShoot = obj.hessShoot;
    hessPath = obj.hessPath;

    prdHss = obj.prdHss;
    
    dxdu = obj.dxdu;
    dxdudx = obj.dxdudx;
    dpc = obj.dpc;
    ddShoot = obj.ddShoot;
    ddPath = obj.ddPath;
    dHshoot = obj.dHshoot;
    
    Ns_ = obj.Ns_;
}

/**
 * Assignment operator
 */
QuadraticModel& QuadraticModel::operator=(const QuadraticModel& obj){
    
    if (this != &obj){
    
        dZ = obj.dZ;
        dZBL = obj.dZBL;
        dZshoot = obj.dZshoot;
        
        z = obj.z;
        zC = obj.zC;
        zS = obj.zS;
        
        grad = obj.grad;
        
        hessShoot = obj.hessShoot;
        hessPath = obj.hessPath;

        prdHss = obj.prdHss;
        
        dxdu = obj.dxdu;
        dxdudx = obj.dxdudx;
        dpc = obj.dpc;
        ddShoot = obj.ddShoot;
        ddPath = obj.ddPath;
        dHshoot = obj.dHshoot;
        
        Ns_ = obj.Ns_;
    }
    
    return *this;
}

/**
 * Destructor
 */
QuadraticModel::~QuadraticModel(){
    
}

/**
* Evaluate model at y, vmod <- gAL'*(y-z) + 0.5*(y-z)'*B(z)*(y-z)
* and update hessProd <- B(z)*(y-z)
*/
double QuadraticModel::evalModel(double* y,double* hessProd){

    int n,i,j; 
    double vmodl,vmodq,vmod,dyz;
    double prod,aux;
    double *gn,*yy,*zz;
    double *phsS;
    double *hsn;
#if DPC
    double *yyp,*zzp;
    double *phsP;
    double *hpn;
#endif

    /* Evaluate gAL'*(y-z) and initialize hessProd to 0 */
    vmodl = 0.;
    gn = grad;
    yy = y;
    zz = z;
    phsS = hessProd;
    while (gn != grad+dZ){
        vmodl += *gn++*(*yy++-*zz++);
        *phsS++ = 0.;
    }

    /* Evaluate 0.5*(y-z)'*B(z)*(y-z) as well as product B(z)*(y-z) */
    hsn = hessShoot;
    phsS = hessProd;
    yy = y;
    zz = z;
    vmodq = 0.;
#if DPC
    phsP = phsS+dZshoot;
    hpn = hessPath;
    yyp = y+dZshoot;
    zzp = z+dZshoot;
#endif
    for (n=0;n<Ns_;++n){
        /* Shooting block of model hessian */
        for (i=0;i<dxdudx;++i){
            dyz = *(yy+i)-*(zz+i);
            /* Hessian inner product */
            prod = 0.;
            for (j=0;j<dxdudx;++j){
                aux = *hsn*(*(yy+j)-*(zz+j));
                prod += aux;
                vmodq += aux*dyz;
                ++hsn;
            }
            *(phsS+i) += prod;
        }
        yy += dxdu;
        zz += dxdu;
        phsS += dxdu;
#if DPC
        /* Path constraint block of model hessian */
        for (i=0;i<dpc;++i){
            dyz = *(yyp+i)-*(zzp+i);
            /* Hessian cross product */
            prod = 0.;
            for (j=0;j<dpc;++j){
                aux = *hpn*(*(yyp+j)-*(zzp+j));
                prod += aux;
                vmodq += aux*dyz;
                ++hpn;
            }
            *phsP++ = prod;
        }
        yyp += dpc;
        zzp += dpc;
#endif
    }
    /* Add-up linear & quadratic part and return */
//    std::cout<<"vmodl="<<vmodl<<std::endl;
//    std::cout<<"vmodq="<<vmodq<<std::endl;
    vmod = vmodl+0.5*vmodq;

    return vmod;
}

/**
* Evaluate model at candidate point vmod <- gal'*(zS-z) + 0.5*(zS-z)'*hal*(zS-z)
*/
double QuadraticModel::evalModelCandidate(){
 
    int k,n,i,j;
    double vmod,vmodl,vmodq,dzz;
    double *zs,*zz;
    double *gn,*hsn;
    double *phsS,prod,aux;
#if DPC
    double *zsp,*zzp;
    double *hpn;
    double *phsP;
#endif

    /* Evaluate gAL'*(zS-z) and initialize prdHss to 0 */
    vmodl = 0.;
    gn = grad;
    zs = zS;
    zz = z;
    phsS = prdHss;
    while (gn != grad+dZ){
//    for (k=0;k<dZ;++k){
        vmodl += *gn++*(*zs++-*zz++);
        *phsS++ = 0.;
//        ++gn;
//        ++zs;
//        ++zz;
    }
    
    /* Evaluate 0.5*(zS-z)'*hAL*(zS-z) as well as product hAL*(zS-z) */
    hsn = hessShoot;
    phsS = prdHss;
    zs = zS;
    zz = z;
    vmodq = 0.;
#if DPC
    hpn = hessPath;
    zsp = zS+dZshoot;
    zzp = z+dZshoot;
    phsP = prdHss+dZshoot;
#endif
    for (n=0;n<Ns_;++n){
        /* Shooting block of approx AL hessian */
        for (i=0;i<dxdudx;++i){
            dzz = *(zs+i)-*(zz+i);
            /* Hessian cross product */
            prod = 0.;
            for (j=0;j<dxdudx;++j){
                aux = *hsn*(*(zs+j)-*(zz+j));
                prod += aux;
                vmodq += aux*dzz;
                ++hsn;
            }
            *(phsS+i) += prod;
        }
        zs += dxdu;
        zz += dxdu;
        phsS += dxdu;

#if DPC
        /* Path constraint block of approx AL hessian */
        for (i=0;i<dpc;++i){
            dzz = *(zsp+i)-*(zzp+i);
            /* Hessian cross product */
            prod = 0.;
            for (j=0;j<dpc;++j){
                aux = *hpn*(*(zsp+j)-*(zzp+j));
                prod += aux;
                vmodq += aux*dzz;
                ++hpn;
            }
            *phsP++ = prod;
        }
        zsp += dpc;
        zzp += dpc;
#endif
    }
    vmod = vmodl+0.5*vmodq;
    return vmod;
}

/**
 * Evaluate model at Cauchy point vmod <- gal'*(zC-z) + 0.5*(zC-z)'*hal*(zC-z)
 */
double QuadraticModel::evalModelCauchy(){

    int n,i,j;
    double vmod,vmodl,vmodq,dzz;
    double *zc,*zz;
    double *gn,*hsn;
    double *phsS,prod,aux;
#if DPC
    double *zcp,*zzp;
    double *hpn;
    double *phsP;
#endif
    /* Evaluate gAL'*(zC-z) and initialize prdHss to 0 */
    gn = grad;
    zc = zC;
    zz = z;
    phsS = prdHss;
    vmodl = 0.;
    while (gn != grad+dZ){
//    for (k=0;k<dZ;++k){
        vmodl += *gn++*(*zc++-*zz++);
        *phsS++ = 0.;
//        ++gn;
//        ++zs;
//        ++zz;
    }
    /* Evaluate 0.5*(zC-z)'*hAL*(zC-z) as well as product hAL*(zC-z) */
    hsn = hessShoot;
    phsS = prdHss;
    zc = zC;
    zz = z;
    vmodq = 0.;
#if DPC
    hpn = hessPath;
    zcp = zC+dZshoot;
    zzp = z+dZshoot;
    phsP = prdHss+dZshoot;
#endif
    for (n=0;n<Ns_;++n){
        /* Shooting block of approx AL hessian */
        for (i=0;i<dxdudx;++i){
            dzz = *(zc+i)-*(zz+i);
            /* Hessian cross product */
            prod = 0.;
            for (j=0;j<dxdudx;++j){
                aux = *hsn*(*(zc+j)-*(zz+j));
                prod += aux;
                vmodq += aux*dzz;
                ++hsn;
            }
            *(phsS+i) += prod;
        }
        zc += dxdu;
        zz += dxdu;
        phsS += dxdu;
#if DPC
        /* Path constraint block of approx AL hessian */
        for (i=0;i<dpc;++i){
            dzz = *(zcp+i)-*(zzp+i);
            /* Hessian cross product */
            prod = 0.;
            for (j=0;j<dpc;++j){
                aux = *hpn*(*(zcp+j)-*(zzp+j));
                prod += aux;
                vmodq += aux*dzz;
                ++hpn;
            }
            *phsP++ = prod;
        }
        zcp += dpc;
        zzp += dpc;
#endif
    }
    vmod = vmodl+0.5*vmodq;
    return vmod;
}

/* Evaluate quadratic term  0.5*(zC-z)'*B(z)*(zC-z) and product B(z)*(zC-z)
CAUTION: prdHss must be initialized to 0. */
double QuadraticModel::evalQuadTermCauchy(){

    int i,j,n;
    double vmodq,dzz,prod,aux;
    double *zc,*zz;
    double *gn,*hsn,*phsS;
#if DPC
    double *zcp,*zzp;
    double *hpn;
    double *phsP;
#endif

    hsn = hessShoot;
    phsS = prdHss;
    zz = z;
    zc = zC;
    vmodq = 0.;
#if DPC
    hpn = hessPath;
    zcp = zc+dZshoot;
    zzp = zz+dZshoot;
    phsP = prdHss+dZshoot;
#endif
    for (n=0;n<Ns_;++n){
        /* Shooting block */
        for (i=0;i<dxdudx;++i){
            dzz = *(zc+i)-*(zz+i);
            prod = 0.;
            for (j=0;j<dxdudx;++j){
                aux = *hsn++*(*(zc+j)-*(zz+j));
                prod += aux;
                vmodq += aux*dzz;
            }
            *(phsS+i) += prod;
        }
        zz += dxdu;
        zc += dxdu;
        phsS += dxdu;
#if DPC
        /* Path slacks block */
        for (i=0;i<dpc;++i){
            dzz = *(zcp+i)-*(zzp+i);
            prod = 0.;
            for (j=0;j<dpc;++j){
                aux = *hpn++*(*(zcp+j)-*(zzp+j));
                prod += aux;
                vmodq += aux*dzz;
            }
            *phsP++ += prod;
        }
        zzp += dpc;
        zcp += dpc;
#endif
    }
    return 0.5*vmodq;
}

/* Evaluate quadratic term 0.5*(y-z)'*B(z)*(y-z) and product B(z)*(y-z)
CAUTION: prdHss must be initialized to 0. */
double QuadraticModel::evalQuadTerm(double* y){

    int i,j,n;
    double vmodq,dyz,prod,aux;
    double *yy,*zz;
    double *gn,*hsn,*phsS;
#if DPC
    double *yyp,*zzp;
    double *hpn;
    double *phsP;
#endif

    hsn = hessShoot;
    phsS = prdHss;
    zz = z;
    yy = y;
    vmodq = 0.;
#if DPC
    hpn = hessPath;
    zzp = zz+dZshoot;
    yyp = yy+dZshoot;
    phsP = phsS+dZshoot;
#endif
    for (n=0;n<Ns_;++n){
        /* Shooting block */
        for (i=0;i<dxdudx;++i){        
            dyz = *(yy+i)-*(zz+i);
            prod = 0.;
            for (j=0;j<dxdudx;++j){
                aux = *hsn++*(*(yy+j)-*(zz+j));
                prod += aux;
                vmodq += aux*dyz;
            }
            *(phsS+i) += prod;
        }
        yy += dxdu;
        zz += dxdu;
        phsS += dxdu;
#if DPC
        /* Path slacks block */
        for (i=0;i<dpc;++i){
            dyz = *(yyp+i)-*(zzp+i);
            prod = 0.;
            for (j=0;j<dpc;++j){
                aux = *hpn++*(*(yyp+j)-*(zzp+j));
                prod += aux;
                vmodq += aux*dyz;
            }
            *(phsP+i) += prod;
        }   
        yyp += dpc;
        zzp += dpc;
#endif
    }
    return 0.5*vmodq;
}

/* (v-w)'*H*(v-w) */
double QuadraticModel::multHshootDD(double* v,double* w){
    
    int n,irw,i,j;
    double res,vsi,wsi,dvw;
    double *hsn,*vs,*ws;
#if DPC
    double vpi,wpi;
    double *hpn,*vp,*wp;
#endif    
    hsn = hessShoot;    
    vs = v;
    ws = w;
#if DPC
    hpn = hessPath;
    vp = v+dZshoot;
    wp = w+dZshoot;
#endif
    res = 0.;
    for (n=0;n<Ns_;++n){  
        /* Shooting block */
        irw = 0;
        for (i=0;i<dxdudx;++i){
            vsi = *(vs+i);
            wsi = *(ws+i);
            dvw = vsi-wsi;
            for (j=0;j<dxdudx;++j){
                res += *(hsn+irw+j)*dvw*(*(vs+j)-*(ws+j));
            }
            irw += dxdudx;
        }
        vs += dxdu;
        ws += dxdu;
        hsn += ddShoot;
#if DPC
        /* Path-constraint block */
        irw = 0;
        for (i=0;i<dpc;++i){
            vpi = *(vp+i);
            wpi = *(wp+i);
            dvw = vpi-wpi;
            for (j=0;j<dpc;++j){
                res += *(hpn+irw+j)*dvw*(*(vp+j)-*(wp+j));
            }
            irw += dpc;
        }
        vp += dpc;
        wp += dpc;
        hpn += ddPath;
#endif 
    }
    return res;
}

/**
** For debugging
**/
/* Displays model gradient */
void QuadraticModel::displayGradient(){

    double *gg;

    gg = grad;
    while (gg != grad+dZ)
        std::cout<<"<<QuadraticModel>>, *gg="<<*gg++<<std::endl;
}

/* Displays hessian of quadratic model */
void QuadraticModel::displayHessian(){

    int n,i,j;
    double *hshoot;
#if DPC
    double *hslack;
#endif

    hshoot = hessShoot;
#if DPC
    hslack = hessPath;
#endif
    for (n=0;n<Ns_;++n){
        
        std::cout<<"Hessian shooting block "<<n<<": "<<std::endl;
        for (i=0;i<dxdudx;++i){
            for (j=0;j<dxdudx;++j){
                std::cout<<*hshoot<<" ";
                ++hshoot;
            }
            std::cout<<""<<std::endl;
        }
#if DPC
        std::cout<<"Hessian slack block"<<std::endl;
        for (i=0;i<dpc;++i){
            for (j=0;j<dpc;++j){
                std::cout<<*hslack<<" ";
                ++hslack;
            }
            std::cout<<""<<std::endl;
        }
#endif
    }
}
