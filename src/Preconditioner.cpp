#include <iostream>
#include <string.h>
#include <math.h>

#include "userData.h"
#include "macros.h"
#include "Preconditioner.h"

/**
** Default constructor
**/
Preconditioner::Preconditioner(){

	pMbShoot = NULL;
	pMbPath = NULL;
	precD = NULL;
	
	hALshootRed = NULL;
	hALpathRed = NULL;

	freeVars = NULL;
    varStatus = NULL;

    ivarFreeState = NULL;
    ivarFreeInput = NULL;
    ivarFreeSlack = NULL;
    ivarFreeGroupShoot = NULL;

    nfreeState = NULL;
    nfreeInput = NULL;
    nfreeSlack = NULL;
    nfreeGroupShoot = NULL;
}

/** 
** Destructor
**/
Preconditioner::~Preconditioner(){

	if (pMbShoot != NULL)
		delete [] pMbShoot;
#if DPC
	if (pMbPath != NULL)
		delete [] pMbPath;
#endif
	if (precD != NULL)
		delete [] precD;
}

/**
** Copy constructor
**/
Preconditioner::Preconditioner(const Preconditioner& obj){

	jacRegCoef = obj.jacRegCoef;
	bw = obj.bw;
	bbw = obj.bbw;
	pivTol = obj.pivTol;

	Ns_ = obj.Ns_;
	dxux = obj.dxux;
	dpc = obj.dpc;

	dZ = obj.dZ;
	dZshoot = obj.dZshoot;
	dZpath = obj.dZpath;
	dGshoot = obj.dGshoot;
	dGpath = obj.dGpath;
	dHshoot = obj.dHshoot;
	dHshootBL = obj.dHshootBL;	
	dHpath = obj.dHpath;
	dHpathBL = obj.dHpathBL;

	dRed = obj.dRed;
	dRedBL = obj.dRedBL;
	dRedShoot = obj.dRedShoot;
	dRedPath = obj.dRedPath;

	freeVars = obj.freeVars;
	varStatus = obj.varStatus;

	ivarFreeState = obj.ivarFreeState;
    ivarFreeInput = obj.ivarFreeInput;
    ivarFreeSlack = obj.ivarFreeSlack;
    ivarFreeGroupShoot = obj.ivarFreeGroupShoot;

    nfreeState = obj.nfreeState;
    nfreeInput = obj.nfreeInput;
    nfreeSlack = obj.nfreeSlack;
    nfreeGroupShoot = obj.nfreeGroupShoot;

    pMbShoot = new double[dHshoot];
    std::copy(obj.pMbShoot,obj.pMbShoot+obj.dHshoot,pMbShoot);
#if DPC
    pMbPath = new double[dHpath];
    std::copy(obj.pMbPath,obj.pMbPath+obj.dHpath,pMbPath);
#endif
    precD = new double[dZ];
    std::copy(obj.precD,obj.precD+obj.dZ,precD);
}

/**
** Assignment operator
**/
Preconditioner& Preconditioner::operator= (const Preconditioner& obj){

	double *precD_,*pMbShoot_,*pMbPath_;

	if (this != &obj){
		
		Ns_ = obj.Ns_;
		dxux = obj.dxux;
		dpc = obj.dpc;

		dZ = obj.dZ;
		dZshoot = obj.dZshoot;
		dZpath = obj.dZpath;
		dGshoot = obj.dGshoot;
		dHshoot = obj.dHshoot;
		dHshootBL = obj.dHshootBL;
		dGpath = obj.dGpath;
		dHpath = obj.dHpath;
		dHpathBL = obj.dHpathBL;

		dRed = obj.dRed;
		dRedBL = obj.dRedBL;
		dRedShoot = obj.dRedShoot;
		dRedPath = obj.dRedPath;

		jacRegCoef = obj.jacRegCoef;
		bw = obj.bw;
		bbw = obj.bbw;
		pivTol = obj.pivTol;

		freeVars = obj.freeVars;
    	varStatus = obj.varStatus;

    	ivarFreeState = obj.ivarFreeState;
    	ivarFreeInput = obj.ivarFreeInput;
    	ivarFreeSlack = obj.ivarFreeSlack;
    	ivarFreeGroupShoot = obj.ivarFreeGroupShoot;

    	nfreeState = obj.nfreeState;
    	nfreeInput = obj.nfreeInput;
    	nfreeSlack = obj.nfreeSlack;
    	nfreeGroupShoot = obj.nfreeGroupShoot;

		/* pMbShoot */
		if (pMbShoot != NULL)
			delete [] pMbShoot;
		pMbShoot_ = new double[dHshoot];
		std::copy(obj.pMbShoot,obj.pMbShoot+obj.dHshoot,pMbShoot_);
		pMbShoot = pMbShoot_;
#if DPC
		/* pMbPath */
		if (pMbPath != NULL)
			delete [] pMbPath;
		pMbPath_ = new double[dHpath];
		std::copy(obj.pMbPath,obj.pMbPath+obj.dHpath,pMbPath_);
		pMbPath = pMbPath_;
#endif
		/* precD */
		if (precD != NULL)
			delete [] precD;
		precD_ = new double[dZ];
		std::copy(obj.precD,obj.precD+obj.dZ,precD_);
		precD = precD_;
	}
	return *this;
}

/**
* Initialize data
*/
void Preconditioner::init(){

	/* pMbShoot */
	if (pMbShoot != NULL)
		delete [] pMbShoot;
	pMbShoot = new double[dHshoot];
	memset(pMbShoot,0,dHshootBL);
#if DPC
	/* pMbPath */
	if (pMbPath != NULL)
		delete [] pMbPath;
	pMbPath = new double[dHpath];
	memset(pMbPath,0,dHpathBL);
#endif
	/* precD */
	if (precD != NULL)
		delete [] precD;
	precD = new double[dZ];
	memset(precD,0,dZ*SZDBL);
}	

/**
** Builds positive definite Jacobi (diagonal) preconditioner from shifted diagonal of reduced hessian approx
**/
void Preconditioner::buildJacobi(){

	int i,n,nfr;
	int *nfgS,*nfP;
	double px;
	double *prDs,*prDp,*pprD;
	double *hsr,*hpr;

	/* Set first dRed preconditioner elements to 0 */
	memset(precD,0,dRedBL);

	/* Extract diagonal of reduced hessian approx */
	prDs = precD;
	prDp = precD+dRedShoot;
	hsr = hALshootRed;
	hpr = hALpathRed;
	nfgS = nfreeGroupShoot;
	nfP = nfreeSlack;
	for (n=0;n<Ns_;++n){
		/* Reduced shooting block n */
		nfr = *nfgS++;
		pprD = prDs;
		for (i=0;i<nfr;++i){
			*pprD++ += *(hsr+i);
			hsr += nfr;
		}
		prDs += *(nfreeState+n)+*(nfreeInput+n);
#if DPC
		/* Reduced slacks block n */
		nfr = *nfP++;
		for (i=0;i<nfr;++i){
			*prDp++ = *(hpr+i);
			hpr += nfr;
		}
#endif
	}
	/* Create positive-definite preconditioner from reduced hessian diagonal */
	prDs = precD;
	while (prDs != precD+dRed){
		px = *prDs;
		if (px<=POSZER)
			*prDs = jacRegCoef;
		prDs++;
//		*prDs++ = 1.;
//		*prDs++ = MAX(px,jacRegCoef);
	}
}

/**
** Applies Jacobi (diagonal) preconditioner 
**/
void Preconditioner::applyJacobi(double* vout,double* vin){
	
	double *prD;
	double *vi,*vo;

	vi = vin;
	vo = vout;
	prD = precD;
	while (prD != precD+dRed)
		*vo++ = *vi++/(*prD++);
}

/**
* Assembles and factorizes band preconditioner (band format) 
**/
void Preconditioner::buildBand(){

	int b,n,i,j,k,pdef,m,ipjm;
	int nfr;
	double *hsr,*hhsr,*hpr,*hhpr;
	double *pMbS,*pMbP,*ppMbS,*ppMbP;
	double *diag,*offdiag,*off,*offl,*offr;
	double gam,aux,tau1,tau2,dg,offsum,offdg;

//	std::cout<<"bbw="<<bbw<<std::endl;

	/*
	* Extract diagonal bands of reduced AL hessian 
	*/
	memset(pMbShoot,0,bbw*dRedShootBL);
	pMbS = pMbShoot;
#if DPC
	memset(pMbPath,0,bbw*dRedPathBL);
	pMbP = pMbPath;
#endif
	for (b=0;b<bbw;++b){
//		std::cout<<"Build band "<<b<<std::endl;
		hsr = hALshootRed+b;
		hpr = hALpathRed+b;
		for (n=0;n<Ns_;++n){
			/* Reduced shooting block n */
			nfr = *(nfreeGroupShoot+n);
			ppMbS = pMbS;
			hhsr = hsr;
			for (i=0;i<nfr-b;++i){
				*ppMbS++ += *(hhsr+i);
				hhsr += nfr;
			}
			hsr += nfr*nfr;
			pMbS += *(nfreeState+n)+*(nfreeInput+n);
#if DPC
			/* Reduced slacks block n */
			nfr = *(nfreeSlack+n);
			ppMbP = pMbP;
			hhpr = hpr;
			for (i=0;i<nfr-b;++i){
				*ppMbP++ = *(hhpr+i);
				hhpr += nfr;
			}
			hpr += nfr*nfr;
			pMbP += nfr;
#endif
		}
		pMbS += *(nfreeState+Ns_);
//		pMbS += *(nfreeState+Ns_)-b;
//		pMbP += -b;
	}

/*	std::cout<<"Reduced hessian: "<<std::endl;
	displayReducedHessian();
	std::cout<<"Before LDL' decomposition of band preconditioner: "<<std::endl;
	displayBand(); */

	/*
	* Compute band LDL' decomposition of shooting preconditioner (based on Lancelot's band module) 
	*/
	/* Check positive-definiteness */
	pdef = 1;
	gam = -1.;
	pMbS = pMbShoot;
	while (pMbS != pMbShoot+dRedShoot){
		aux = *pMbS++;
		if (aux<=POSZER)
			pdef = 0;
		gam = MAX(gam,ABS(aux));
	}
	/* Factorization loop */
	tau1 = gam*pow(pivTol,0.333); // Pivot tolerance
	tau2 = tau1;
	diag = pMbShoot;
	offdiag = diag+dRedShoot;
	/* Iterate over columns */
	for (i=0;i<dRedShoot;++i){ 
		dg = *(diag+i);
		m = MIN(bw,dRedShoot-i-1);
		if (pdef){
			off = offdiag;
			for (j=0;j<m;++j){
				aux = *(off+i);
				if ((*(diag+i+j+1)-aux*aux/dg)<=tau1){
					pdef = 0;
					break;
				}
				off += dRedShoot;
			}
		}
		/* If not positive definite, perturb diagonal to make it positive */
		if (pdef==0){
			offsum = 0.;
			off = offdiag;
			for (j=0;j<m;++j){
				offsum += ABS(*(off+i));
				off += dRedShoot;
			}
			offsum = MAX(POSZER,-dg+MAX(offsum,tau2));
			*(diag+i) = dg+offsum;
		}	
		/* Continue factorization */
		dg = *(diag+i);
		off = offdiag;
		for (j=0;j<m;++j){
			offdg = *(off+i);			
		//	ipjm = j;
			offl = off;
			offr = offdiag;
			for (k=0;k<j-1;++k){	
//				*(offdiag+(ipjm-k)*dRedShoot+i+k) -= offdg**(offdiag+k*dRedShoot+i);
				*(offl+i+k) -= offdg**(offr+i);
				offl -= dRedShoot;
				offr += dRedShoot;
			}
			offdg /= dg;
			*(diag+i+j+1) -= offdg**(off+i);
			*(off+i) = offdg;
			off += dRedShoot; 
		}
	}

#if DPC
	/* 
	* Compute band LDL' decomposition of slack preconditioner (based on Lancelot's band module)
	*/
	/* Check positive-definiteness */
	pdef = 1;
	gam = -1.;
	pMbP = pMbPath;
	while (pMbP != pMbPath+dRedPath){
		aux = *pMbP++;
		if (aux<=POSZER)
			pdef = 0;
		gam = MAX(gam,ABS(aux));
	}
	/* Factorization loop */
	tau1 = gam*pow(pivTol,0.333); // Pivot tolerance
	tau2 = tau1;
	diag = pMbPath;
	offdiag = diag+dRedPath;
	/* Iterate over columns */
	for (i=0;i<dRedPath;++i){ 
		dg = *(diag+i);
		m = MIN(bw,dRedPath-i-1);
		if (pdef){
			for (j=0;j<m;++j){
				aux = *(offdiag+j*dRedPath+i);
				if ((*(diag+i+j+1)-aux*aux/dg)<=tau1){
					pdef = 0;
					break;
				}
			}
		}
		/* If not positive definite, perturb diagonal to make it positive */
		if (pdef==0){
			offsum = 0.;
			for (j=0;j<m;++j)
				offsum += ABS(*(offdiag+j*dRedPath+i));
			offsum = MAX(POSZER,-dg+MAX(offsum,tau2));
			*(diag+i) = dg+offsum;
		}	
		/* Continue factorization */
		dg = *(diag+i);
		for (j=0;j<m;++j){
			offdg = *(offdiag+j*dRedPath+i);			
			ipjm = j;
			for (k=0;k<j-1;++k)	
				*(offdiag+(ipjm-k)*dRedPath+i+k) -= offdg**(offdiag+k*dRedPath+i);
			offdg /= dg;
			*(diag+i+j+1) -= offdg**(offdiag+j*dRedPath+i);
			*(offdiag+j*dRedPath+i) = offdg; 
		}
	}
#endif

//	std::cout<<"After LDL' decomposition of band preconditioner: "<<std::endl;
//	displayBand();
}

/**
* Applies band preconditioner (forward-backward solve using band factor), based on Lancelot's band module
**/
void Preconditioner::applyBand(double* vout,double* vin){

	int i,j,n,m;
	double aux;
	double *diag,*offdiag;
	double *vv,*off;

	memcpy(vout,vin,dRedBL); // Check if PCG could be modified to avoid this !!!!

	/**
	* Apply shooting band preconditioner
	*/
	vv = vout;
	diag = pMbShoot;
	offdiag = diag+dRedShoot;
	/* Forward and diagonal solve */
	for (i=0;i<dRedShoot;++i){
		m = MIN(bw,dRedShoot-i-1);
		aux = *(vv+i);
		off = offdiag;
		for (j=0;j<m;++j){
			*(vv+i+j+1) -= *(off+i)*aux;
			off += dRedShoot;
	//		off += dRedShoot-j-1; 
		}
		*(vv+i) /= *(diag+i);
	}
	/* Backward solve */
	for (i=dRedShoot-1;i>=0;--i){
		m = MIN(bw,dRedShoot-i-1);
		aux = *(vv+i);
		off = offdiag;
		for (j=0;j<m;++j){
			aux -= *(off+i)**(vv+i+j+1);
			off += dRedShoot;
	//		off += dRedShoot-j-1;
		}
		*(vv+i) = aux;
	}

#if DPC
	/**
	* Apply slack band preconditioner
	*/
	vv = vout+dRedShoot;
	diag = pMbPath;
	offdiag = diag+dRedPath;
	/* Forward and diagonal solve */
	for (i=0;i<dRedPath;++i){
		m = MIN(bw,dRedPath-i-1);
		aux = *(vv+i);
	//	off = offdiag;
		for (j=0;j<m;++j){
			*(vv+i+j+1) -= *(offdiag+j*dRedPath+i)*aux;
	//		off += dRedPath-j-1;
		}
		*(vv+i) /= *(diag+i);
	}
	/* Backward solve */
	for (i=dRedPath-1;i>=0;--i){
		m = MIN(bw,dRedPath-i-1);
		aux = *(vv+i);
	//	off = offdiag;
		for (j=0;j<m;++j){
			aux -= *(offdiag+j*dRedPath+i)**(vv+i+j+1);
	//		off += dRedPath-j-1;
		}
		*(vv+i) = aux;
	}
#endif
}

/* Displays band preconditioner */
void Preconditioner::displayBand(){

	int b,i;
	double *pMbS,*pMbP;

	pMbS = pMbShoot;
	std::cout<<"Shooting band preconditioner: "<<std::endl;
	for (b=0;b<bbw;++b){
		std::cout<<"Band "<<b<<std::endl;
		for (i=0;i<dRedShoot;++i)
			std::cout<<*pMbS++<<" ";
		std::cout<<""<<std::endl;
	}

	pMbP = pMbPath;
	std::cout<<"Slacks band preconditioner: "<<std::endl;
	for (b=0;b<bbw;++b){
		std::cout<<"Band "<<b<<std::endl;
		for (i=0;i<dRedPath;++i)
			std::cout<<*pMbP++<<" ";
		std::cout<<""<<std::endl;
	}
}

/* Displays Jacobi preconditioner */
void Preconditioner::displayJacobi(){
	double *prD;

	prD = precD;
	std::cout<<"Jacobi preconditioner:"<<std::endl;
	while (prD != precD+dRed)
		std::cout<<*prD++<<std::endl;
}

/* Displays reduced hessian */
void Preconditioner::displayReducedHessian(){

	int n,i,j,nfr;
	int *nfrS,*nfrP;
	double *hsr,*hpr;

	hsr = hALshootRed;
	hpr = hALpathRed;
	nfrS = nfreeGroupShoot;
	nfrP = nfreeSlack;
	for (n=0;n<Ns_;++n){
		nfr = *nfrS++;
		std::cout<<"Shooting block "<<n<<" , "<<nfr<<" free vars: "<<std::endl;
		for (i=0;i<nfr;++i){
			for (j=0;j<nfr;++j)
				std::cout<<*hsr++<<" ";
			std::cout<<""<<std::endl;
		}
		std::cout<<"Slack block "<<n<<": "<<std::endl;
		nfr = *nfrP++;
		for (i=0;i<nfr;++i){
			for (j=0;j<nfr;++j)
				std::cout<<*hpr++<<" ";
			std::cout<<""<<std::endl;
		}
	}
}