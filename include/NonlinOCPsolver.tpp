#include <iostream>
#include <iomanip>
#include <string.h>
#include <algorithm>
#include <time.h>

#include "NonlinOCPsolver.h"
#include "macros.h"
#include "userData.h"

/**
** Default constructor
**/
template <class dynT,class costT,class pconT,class mayerT> NonlinOCPsolver<dynT,costT,pconT,mayerT>::NonlinOCPsolver(){
	
	z = NULL;
    zS = NULL;
    negstp = NULL;

    econ = NULL;
	mu = NULL;
	uopt = NULL;
#ifdef BIOREC
    uxAux = NULL;
    xAux = NULL;
    pcAux = NULL;
    pcAuxx = NULL;
#endif

    zmin = NULL;
    zmax = NULL;
    
    gAL = NULL;
    gALshoot = NULL;
    gALshootS = NULL;
    gALpath = NULL;
    gALpathS = NULL;
	dffGshoot = NULL;
	dffGpath = NULL;

    hALshoot = NULL;
    hALshootRed = NULL;
    hALpath = NULL;
    hALpathRed = NULL;
    prdHss = NULL;
    
    indFreeVars = NULL;
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

    secVecShoot = NULL;
    secVecPath = NULL;

    errPrim = 0;
}

/**
** Destructor
**/
template <class dynT,class costT,class pconT,class mayerT> NonlinOCPsolver<dynT,costT,pconT,mayerT>::~NonlinOCPsolver(){

//    std::cout<<"~NonlinOCPsolver()"<<std::endl;

	if (z != NULL)
		delete [] z;

    if (zS != NULL)
        delete [] zS;
    
    if (econ != NULL)
        delete [] econ;
     	
	if (uopt != NULL)
		delete [] uopt;

#ifdef BIOREC
    if (uxAux != NULL)
        delete [] uxAux;

    if (xAux != NULL)
        delete [] xAux;

    if (pcAux != NULL)
        delete [] pcAux;
    
    if (pcAuxx != NULL)
        delete [] pcAuxx;
#endif

    if (zmin != NULL)
        delete [] zmin;
    
    if (zmax != NULL)
        delete [] zmax;
     
    if (gAL != NULL)
        delete [] gAL;
    
    if (gALshoot != NULL)    
        delete [] gALshoot;

    if (gALshootS != NULL)
        delete [] gALshootS;

    if (gALpathS != NULL)
        delete [] gALpathS;

	if (dffGshoot != NULL)
		delete [] dffGshoot;

	if (dffGpath != NULL)
		delete [] dffGpath;

    if (hALshoot != NULL)
        delete [] hALshoot;

    if (hALshootRed != NULL)
        delete [] hALshootRed;
    
    if (hALpath != NULL)
		delete [] hALpath;

    if (hALpathRed != NULL)
        delete [] hALpathRed;

    if (prdHss != NULL)
        delete [] prdHss;

    if (freeVars != NULL)
        delete [] freeVars;

    if (indFreeVars != NULL)
        delete [] indFreeVars;

    if (varStatus != NULL)
        delete [] varStatus;

    if (ivarFreeState != NULL)
        delete [] ivarFreeState;

    if (ivarFreeInput != NULL)
        delete [] ivarFreeInput;

    if (ivarFreeSlack != NULL)
        delete [] ivarFreeSlack;
   
    if (nfreeState != NULL)
        delete [] nfreeState;
    
    if (nfreeInput != NULL)
        delete [] nfreeInput;
    
    if (nfreeSlack != NULL)
        delete [] nfreeSlack;
    
    if (nfreeGroupShoot != NULL)
        delete [] nfreeGroupShoot;

    if (ivarFreeGroupShoot != NULL)
        delete [] ivarFreeGroupShoot;

    if (negstp != NULL)
        delete [] negstp;
     
    if (mu != NULL)
        delete [] mu;
     
    if (secVecShoot != NULL)
        delete [] secVecShoot;
     
    if (secVecPath != NULL)
        delete [] secVecPath;

//    std::cout<<"end of solver dest"<<std::endl;
}	  

/**
** Copy constructor
**/ 
template <class dynT,class costT,class pconT,class mayerT> 
NonlinOCPsolver<dynT,costT,pconT,mayerT>::NonlinOCPsolver(const NonlinOCPsolver<dynT,costT,pconT,mayerT>& obj) : shooter(obj.shooter),
                                                                                                                 activDetector(obj.activDetector),
                                                                                                                 precSubRefiner(obj.precSubRefiner),
                                                                                                                 precProjRefiner(obj.precProjRefiner),
                                                                                                                 qmod(obj.qmod){

    sini = obj.sini;

    pivTol = obj.pivTol;

    eta0 = obj.eta0;
    eta1 = obj.eta1;
    eta2 = obj.eta2;

    sig1 = obj.sig1;
    sig2 = obj.sig2;
    sig3 = obj.sig3;

    maxPit = obj.maxPit;
    maxDit = obj.maxDit;
    maxCumCG = obj.maxCumCG;

    mulPen = obj.mulPen;
    halfRho = obj.halfRho;

    Ns_ = obj.Ns_;

    dZ = obj.dZ;
    dZBL = obj.dZBL;
    dZshoot = obj.dZshoot;
    dZshootBL = obj.dZshootBL;
    dZpath = obj.dZpath;
    dZpathBL = obj.dZpathBL;

    dLshoot = obj.dLshoot;
    dLpath = obj.dLpath;
    dL = obj.dL;
    dLBL = obj.dLBL;

    dGshoot = obj.dGshoot;
    dGshootBL = obj.dGshootBL;
    dGshootNT = obj.dGshootNT;
	dGpath = obj.dGpath;
	dGpathBL = obj.dGpathBL;
    dGpathNT = obj.dGpathNT;

	dhS = obj.dhS;
	dhP = obj.dhP;
    dHshoot = obj.dHshoot;
    dHpath = obj.dHpath;
    
    du = obj.du;
    dubl = obj.dubl;
	dx = obj.dx;
    dxbl = obj.dxbl;
	dxu = obj.dxu;
	dxux = obj.dxux;
	dpc = obj.dpc;
    
    kkt = obj.kkt;
    kkt2 = obj.kkt2;
    kktTolAbs = obj.kktTolAbs;
    nec = obj.nec;
    nec2 = obj.nec2;
    necTolAbs = obj.necTolAbs;
    cgTol2 = obj.cgTol2;
    trRad = obj.trRad;
    norIstp = obj.norIstp;
    nor2stp = obj.nor2stp;

    cumCG = obj.cumCG;
    cumRestart = obj.cumRestart;
    errPrim = obj.errPrim;

    z = new double[dZ];
    std::copy(obj.z,obj.z+obj.dZ,z);

    zS = new double[dZ];
    std::copy(obj.zS,obj.zS+obj.dZ,zS);

	negstp = new double[dZ];
	std::copy(obj.negstp,obj.negstp+obj.dZ,dZ);
    
    econ = new double[dL];
    std::copy(obj.econ,obj.econ+obj.dL,econ);

    mu = new double[dL];
    std::copy(obj.mu,obj.mu+obj.dL,mu);
    
    uopt = new double[du];
    std::copy(obj.uopt,obj.uopt+obj.du,uopt);

#ifdef BIOREC
    uxAux = new double[dxu];
    std::copy(obj.uxAux,obj.uxAux+obj.dxu,uxAux);

    xAux = new double[dx];
    std::copy(obj.xAux,obj.xAux+obj.dx,xAux);

    pcAux = new double[dpc];
    std::copy(obj.pcAux,obj.pcAux+obj.dpc,pcAux);

    pcAuxx = new double[dpc];
    std::copy(obj.pcAuxx,obj.pcAuxx+obj.dpc,pcAuxx);
#endif

    zmin = new double[dZ];
    std::copy(obj.zmin,obj.zmin+obj.dZ,zmin);
    
    zmax = new double[dZ];
    std::copy(obj.zmax,obj.zmax+obj.dZ,zmax);
    
    gAL = new double[dZ];
    std::copy(obj.gAL,obj.gAL+obj.dZ,gAL);

    gALshoot = new double[dGshoot];
    std::copy(obj.gALshoot,obj.gALshoot+obj.dGshoot,gALshoot);

    gALshootS = new double[dGshoot];
    std::copy(obj.gALshootS,obj.gALshootS+obj.dGshoot,gALshootS);

    gALpath = obj.gALpath;

    gALpathS = new double[dGpath];
    std::copy(obj.gALpathS,obj.gALpathS+obj.dGpath,gALpathS);

	dffGshoot = new double[dGshoot];
	std::copy(obj.dffGshoot,obj.dffGshoot+dGshoot,dffGshoot);
	
	dffGpath = new double[dGpath];
	std::copy(obj.dffGpath,obj.dffGpath+dGpath,dffGpath);

    hALshoot = new double[dHshoot];
    std::copy(obj.hALshoot,obj.hALshoot+obj.dHshoot,hALshoot);

    hALshootRed = new double[dHshoot];
    std::copy(obj.hALshootRed,obj.hALshootRed+obj.dHshoot,hALshootRed);

    hALpath = new double[dHpath];
    std::copy(obj.hALpath,obj.hALpath+obj.dHpath,hALpath);

    hALpathRed = new double[dHpath];
    std::copy(obj.hALpathRed,obj.hALpathRed+obj.dHpath,hALpathRed);

    prdHss = new double[dZ];
    std::copy(obj.prdHss,obj.prdHss+obj.dZ,prdHss);
    
    freeVars = new int[dZ];
    std::copy(obj.freeVars,obj.freeVars+obj.dZ,freeVars);

    indFreeVars = new int[dZ];
    std::copy(obj.indFreeVars,obj.indFreeVars+obj.dZ,indFreeVars);

    varStatus = new int[dZ];
    std::copy(obj.varStatus,obj.varStatus+obj.dZ,varStatus);

    ivarFreeState = new int[(Ns_+1)*dx];
    std::copy(obj.ivarFreeState,obj.ivarFreeState+(obj.Ns_+1)*obj.dx,ivarFreeState);

    ivarFreeInput = new int[Ns_*du];
    std::copy(obj.ivarFreeInput,obj.ivarFreeInput+obj.Ns_*obj.du,ivarFreeInput);

    ivarFreeSlack = new int[Ns_*dpc];
    std::copy(obj.ivarFreeSlack,obj.ivarFreeSlack+obj.Ns_*obj.dpc,ivarFreeSlack);

    ivarFreeGroupShoot = new int[Ns_*dxux];
    std::copy(obj.ivarFreeGroupShoot,obj.ivarFreeGroupShoot+obj.Ns_*obj.dxux,ivarFreeGroupShoot);

    nfreeState = new int[Ns_+1];
    std::copy(obj.nfreeState,obj.nfreeState+obj.Ns_+1,nfreeState);

    nfreeInput = new int[Ns_];
    std::copy(obj.nfreeInput,obj.nfreeInput+obj.Ns_,nfreeInput);

    nfreeSlack = new int[Ns_];
    std::copy(obj.nfreeSlack,obj.nfreeSlack+obj.Ns_,nfreeSlack);

    nfreeGroupShoot = new int[Ns_];
    std::copy(obj.nfreeGroupShoot,obj.nfreeGroupShoot+obj.Ns_,nfreeGroupShoot);

    secVecShoot = new double[dxux];
    std::copy(obj.secVecShoot,obj.secVecShoot+obj.dxux,secVecShoot);

    secVecPath = new double[dpc];
    std::copy(obj.secVecPath,obj.secVecPath+obj.dpc,secVecPath);

}

/**
* Assignment operator
*/
template <class dynT,class costT,class pconT,class mayerT> 
NonlinOCPsolver<dynT,costT,pconT,mayerT>& NonlinOCPsolver<dynT,costT,pconT,mayerT>::operator=(const NonlinOCPsolver<dynT,costT,pconT,mayerT>& obj){

    double *z_,*negstp_,*zS_,*mu_,*econ_,*uopt_;
#ifdef BIOREC
    double *uxAux_,*pcAux_;
    double *xAux_,*pcAuxx_;
#endif
    double *zmin_,*zmax_;
    double *gAL_,*gALshoot_,*gALshootS_,*dffGshoot_;
    double *gALpathS_,*dffGpath_;
    double *hALshoot_,*hALpath_,*hALshootRed_,*hALpathRed_;
    double *prdHss_;
    int *indFreeVars_,*freeVars_,*varStatus_;
    int *ivarFreeState_,*ivarFreeInput_,*ivarFreeSlack_,*ivarFreeGroupShoot_;
    int *nfreeState_,*nfreeInput_,*nfreeSlack_,*nfreeGroupShoot_;
    double *secVecShoot_,*secVecPath_;

    if (this != &obj){
        shooter = obj.shooter;
        activDetector = obj.activDetector;
        precSubRefiner = obj.precSubRefiner;
        precProjRefiner = obj.precProjRefiner;
        qmod = obj.qmod;

        sini = obj.sini;

        pivTol = obj.pivTol;

        eta0 = obj.eta0;
        eta1 = obj.eta1;
        eta2 = obj.eta2;
        
        sig1 = obj.sig1;
        sig2 = obj.sig2;
        sig3 = obj.sig3;
        
        maxPit = obj.maxPit;
        maxDit = obj.maxDit;
        maxCumCG = obj.maxCumCG;

        mulPen = obj.mulPen;
        halfRho = obj.halfRho;

        trRad = obj.trRad;
        
        kkt = obj.kkt;
        kkt2 = obj.kkt2;
        kktTolAbs = obj.kktTolAbs;
        nec = obj.nec;
        nec2 = obj.nec2;
        necTolAbs = obj.necTolAbs;
        cgTol2 = obj.cgTol2;
        norIstp = obj.norIstp;
        nor2stp = obj.nor2stp;

        cumCG = obj.cumCG;
        cumRestart = obj.cumRestart;
        errPrim = obj.errPrim;

        Ns_ = obj.Ns_;
        Ns = obj.Ns;

        dZshoot = obj.dZshoot;
        dZshootBL = obj.dZshootBL;
        dZpath = obj.dZpath;
        dZpathBL = obj.dZpathBL;
        dZ = obj.dZ;
        dZBL = obj.dZBL;
        
        dLshoot = obj.dLshoot;
        dLpath = obj.dLpath;
        dL = obj.dL;
        dLBL = obj.dLBL;
        
        du = obj.du;
        dubl = obj.dubl;
        dx = obj.dx;
        dxbl = obj.dxbl;
        dxu = obj.dxu;
        dxux = obj.dxux;
        dpc = obj.dpc;

        dGshoot = obj.dGshoot;
        dGshootBL = obj.dGshootBL;
        dGshootNT = obj.dGshootNT;

        dGpath = obj.dGpath;
        dGpathBL = obj.dGpathBL;
        dGpathNT = obj.dGpathNT;

        dhS = obj.dhS;
        dhP = obj.dhP;
        dHshoot = obj.dHshoot;
        dHpath = obj.dHpath;

        /* z */
        if (z != NULL)
            delete [] z;
        if (obj.z != NULL){
            z_ = new double[dZ];
            std::copy(obj.z,obj.z+obj.dZ,z_);
            z = z_;
        }
        else {
            z = NULL;
        }

        /* negstp */
        if (negstp != NULL)
            delete [] negstp;
        if (obj.negstp != NULL){
            negstp_ = new double[dZ];
            std::copy(obj.negstp,obj.negstp+obj.dZ,negstp_);
            negstp = negstp_;
        }
        else {
            negstp = NULL;
        }

        /* zS */
        if (zS != NULL)
            delete [] zS;
        if (obj.zS != NULL){
            zS_ = new double[dZ];
            std::copy(obj.zS,obj.zS+obj.dZ,zS_);
            zS = zS_;
        }
        else {
            zS = NULL;
        }

        /* mu */
        if (mu != NULL)
            delete [] mu;
        if (obj.mu != NULL){
            mu_ = new double[dL];
            std::copy(obj.mu,obj.mu+obj.dL,mu_);
            mu = mu_;
        }
        else {
            mu = NULL;
        }

        /* econ */
        if (econ != NULL)
            delete [] econ;
        if (obj.econ != NULL){
            econ_ = new double[dL];
            std::copy(obj.econ,obj.econ+obj.dL,econ_);
            econ = econ_;
        }
        else {
            econ = NULL;
        }

        /* uopt */
        if (uopt != NULL)
            delete [] uopt;
        if (obj.uopt != NULL){
            uopt_ = new double[du];
            std::copy(obj.uopt,obj.uopt+obj.du,uopt_);
            uopt = uopt_;
        }
        else {
            uopt = NULL;
        }

#ifdef BIOREC
        /* uxAux */
        if (uxAux != NULL)
            delete [] uxAux;
        if (obj.uxAux != NULL){
            uxAux_ = new double[dxu];
            std::copy(obj.uxAux,obj.uxAux+obj.dxu,uxAux_);
            uxAux = uxAux_;
        }   
        else {
            uxAux = NULL;
        }

        /* xAux */
        if (xAux != NULL)
            delete [] xAux;
        if (obj.xAux != NULL){
            xAux_ = new double[dx];
            std::copy(obj.xAux,obj.xAux+obj.dx,xAux_);
            xAux = xAux_;
        }
        else {
            xAux = NULL;
        }

        /* pcAux */
        if (pcAux != NULL)
            delete [] pcAux;
        if (obj.pcAux != NULL){
            pcAux_ = new double[dpc];
            std::copy(obj.pcAux,obj.pcAux+obj.dpc,pcAux_);
            pcAux = pcAux_;
        }
        else {
            pcAux = NULL;
        }

        /* pcAuxx */
        if (pcAuxx != NULL)
            delete [] pcAuxx;
        if (obj.pcAuxx != NULL){
            pcAuxx_ = new double[dpc];
            std::copy(obj.pcAuxx,obj.pcAuxx+obj.dpc,pcAuxx_);
            pcAuxx = pcAuxx_;
        }
        else {
            pcAuxx = NULL;
        }
#endif

        /* zmin */
        if (zmin != NULL)
            delete [] zmin;
        if (obj.zmin != NULL){    
            zmin_ = new double[dZ];
            std::copy(obj.zmin,obj.zmin+obj.dZ,zmin_);
            zmin = zmin_;
        }
        else {
            zmin = NULL;
        }
        
        /* zmax */
        if (zmax != NULL)
            delete [] zmax;
        if (obj.zmax != NULL){
            zmax_ = new double[dZ];
            std::copy(obj.zmax,obj.zmax+obj.dZ,zmax_);
            zmax = zmax_;
        }
        else {
            zmax = NULL;
        }
        
        /* gAL */
        if (gAL != NULL)
            delete [] gAL;

        if (obj.gAL != NULL){
            gAL_ = new double[dZ];
            std::copy(obj.gAL,obj.gAL+obj.dZ,gAL_);
            gAL = gAL_;
        }
        else {
            gAL = NULL;
        }

        /* gALshoot */
        if (gALshoot != NULL)
            delete [] gALshoot;
        if (obj.gALshoot != NULL){
            gALshoot_ = new double[dGshoot];
            std::copy(obj.gALshoot,obj.gALshoot+obj.dGshoot,gALshoot_);
            gALshoot = gALshoot_;
        }
        else {
            gALshoot = NULL;
        }

        /* gALshootS */
        if (gALshootS != NULL)
            delete [] gALshootS;
        if (obj.gALshootS != NULL){
            gALshootS_ = new double[dGshoot];
            std::copy(obj.gALshootS,obj.gALshootS+obj.dGshoot,gALshootS_);
            gALshootS = gALshootS_;
        }
        else {
            gALshootS = NULL;
        }

        /* dffGshoot */
        if (dffGshoot != NULL)
            delete [] dffGshoot;
        if (obj.dffGshoot != NULL){    
            dffGshoot_ = new double[dGshoot];
            std::copy(obj.dffGshoot,obj.dffGshoot+obj.dGshoot,dffGshoot_);
            dffGshoot = dffGshoot_;
        } 
        else {
            dffGshoot = NULL;
        }

        /* gALpath */
        gALpath = obj.gALpath;

        /* gALpathS */
        if (gALpathS != NULL)
            delete [] gALpathS;
        if (obj.gALpathS != NULL){
            gALpathS_ = new double[dGpath];
            std::copy(obj.gALpathS,obj.gALpathS+obj.dGpath,gALpathS_);
            gALpathS = gALpathS_;
        }
        else {
            gALpathS = NULL;
        }

        /* dffGpath */
        if (dffGpath != NULL)
            delete [] dffGpath;
        if (obj.dffGpath != NULL){ 
            dffGpath_ = new double[dGpath];
            std::copy(obj.dffGpath,obj.dffGpath+obj.dGpath,dffGpath_);
            dffGpath = dffGpath_;
        }
        else {
            dffGpath = NULL;
        }

        /* hALshoot */
        if (hALshoot != NULL)
            delete [] hALshoot;
        if (obj.hALshoot != NULL){
            hALshoot_ = new double[dHshoot];
            std::copy(obj.hALshoot,obj.hALshoot+obj.dHshoot,hALshoot_);
            hALshoot = hALshoot_;
        }
        else {
            hALshoot = NULL;
        }

        /* hALshootRed */
        if (hALshootRed != NULL)
            delete [] hALshootRed;
        if (obj.hALshootRed != NULL){
            hALshootRed_ = new double[dHshoot];
            std::copy(obj.hALshootRed,obj.hALshootRed+obj.dHshoot,hALshootRed_);
            hALshootRed = hALshootRed_;
        }
        else {
            hALshootRed = NULL;
        }

        /* hALpath */
        if (hALpath != NULL)
            delete [] hALpath;

        if (obj.hALpath != NULL){
            hALpath_ = new double[dHpath];
            std::copy(obj.hALpath,obj.hALpath+obj.dHpath,hALpath_);
            hALpath = hALpath_;
        }
        else {
            hALpath = NULL;
        }

        /* hALpathRed */
        if (hALpathRed != NULL)
            delete [] hALpathRed;
        if (obj.hALpathRed != NULL){
            hALpathRed_ = new double[dHpath];
            std::copy(obj.hALpathRed,obj.hALpathRed+obj.dHpath,hALpathRed_);
            hALpathRed = hALpathRed_;
        }
        else {
            hALpathRed = NULL;
        }  

        /* prdHss */
        if (prdHss != NULL)
            delete [] prdHss;
        if (obj.prdHss != NULL){
            prdHss_ = new double[dZ];
            std::copy(obj.prdHss,obj.prdHss+obj.dZ,prdHss_);
            prdHss = prdHss_;
        }
        else {
            prdHss = NULL;
        }

        /* secVecShoot */
        if (secVecShoot != NULL)
            delete [] secVecShoot;
        if (obj.secVecShoot != NULL){
            secVecShoot_ = new double[dxux];
            std::copy(obj.secVecShoot,obj.secVecShoot+obj.dxux,secVecShoot_);
            secVecShoot = secVecShoot_;
        }
        else {
            secVecShoot = NULL;
        }

        /* secVecPath */
        if (secVecPath != NULL)
            delete [] secVecPath;
        if (obj.secVecPath != NULL){
            secVecPath_ = new double[dpc];
            std::copy(obj.secVecPath,obj.secVecPath+obj.dpc,secVecPath_);
            secVecPath = secVecPath_;
        }
        else {  
            secVecPath = NULL;
        }

        /* freeVars */
        if (freeVars != NULL)
            delete [] freeVars;
        if (obj.freeVars != NULL){
            freeVars_ = new int[dZ];
            std::copy(obj.freeVars,obj.freeVars+obj.dZ,freeVars_);
            freeVars = freeVars_;
        }
        else {
            freeVars = NULL;
        }

        /* indFreeVars */
        if (indFreeVars != NULL)
            delete [] indFreeVars;
        if (obj.indFreeVars != NULL){
            indFreeVars_ = new int[dZ];
            std::copy(obj.indFreeVars,obj.indFreeVars+obj.dZ,indFreeVars_);
            indFreeVars = indFreeVars_;
        }
        else {
            indFreeVars = NULL;
        }

        /* varStatus */
        if (varStatus != NULL)
            delete [] varStatus;
        if (obj.varStatus != NULL){
            varStatus_ = new int[dZ];
            std::copy(obj.varStatus,obj.varStatus+obj.dZ,varStatus_);
            varStatus = varStatus_;
        }
        else {
            varStatus = NULL;
        }

        /* ivarFreeState */
        if (ivarFreeState != NULL)
            delete [] ivarFreeState;
        if (obj.ivarFreeState != NULL){
            ivarFreeState_ = new int[Ns*dx];
            std::copy(obj.ivarFreeState,obj.ivarFreeState+obj.Ns*obj.dx,ivarFreeState_);   
            ivarFreeState = ivarFreeState_;
        }
        else {
            ivarFreeState = NULL;
        }

        /* ivarFreeInput */
        if (ivarFreeInput != NULL)
            delete [] ivarFreeInput;
        if (obj.ivarFreeInput != NULL){
            ivarFreeInput_ = new int[Ns_*du];
            std::copy(obj.ivarFreeInput,obj.ivarFreeInput+obj.Ns_*obj.du,ivarFreeInput_);   
            ivarFreeInput = ivarFreeInput_;
        }
        else {
            ivarFreeInput = NULL;
        }

        /* ivarFreeSlack */
        if (ivarFreeSlack != NULL)
            delete [] ivarFreeSlack;
        if (obj.ivarFreeSlack != NULL){
            ivarFreeSlack_ = new int[Ns_*dpc];
            std::copy(obj.ivarFreeSlack,obj.ivarFreeSlack+obj.Ns_*obj.dpc,ivarFreeSlack_);   
            ivarFreeSlack = ivarFreeSlack_;
        }
        else {
            ivarFreeSlack = NULL;
        }

        /* ivarFreeGroupShoot */
        if (ivarFreeGroupShoot != NULL)
            delete [] ivarFreeGroupShoot;
        if (obj.ivarFreeGroupShoot != NULL){
            ivarFreeGroupShoot_ = new int[Ns_*dxux];
            std::copy(obj.ivarFreeGroupShoot,obj.ivarFreeGroupShoot+obj.Ns_*obj.dxux,ivarFreeGroupShoot_);   
            ivarFreeGroupShoot = ivarFreeGroupShoot_;
        }
        else {
            ivarFreeGroupShoot = NULL;
        }

        /* nfreeState */
        if (nfreeState != NULL)
            delete [] nfreeState;
        if (obj.nfreeState != NULL){
            nfreeState_ = new int[Ns];   
            std::copy(obj.nfreeState,obj.nfreeState+obj.Ns,nfreeState_);
            nfreeState = nfreeState_;
        }
        else {
            nfreeState = NULL;
        }

        /* nfreeInput */
        if (nfreeInput != NULL)
            delete [] nfreeInput;
        if (obj.nfreeInput != NULL){
            nfreeInput_ = new int[Ns_];   
            std::copy(obj.nfreeInput,obj.nfreeInput+obj.Ns_,nfreeInput_);
            nfreeInput = nfreeInput_;
        }
        else {
            nfreeInput = NULL;
        }

        /* nfreeSlack */
        if (nfreeSlack != NULL)
            delete [] nfreeSlack;
        if (obj.nfreeSlack != NULL){
            nfreeSlack_ = new int[Ns_];   
            std::copy(obj.nfreeSlack,obj.nfreeSlack+obj.Ns_,nfreeSlack_);
            nfreeSlack = nfreeSlack_;
        }
        else {
            nfreeSlack = NULL;
        }

        /* nfreeGroupShoot */
        if (nfreeGroupShoot != NULL)
            delete [] nfreeGroupShoot;
        if (obj.nfreeGroupShoot != NULL){
            nfreeGroupShoot_ = new int[Ns_];   
            std::copy(obj.nfreeGroupShoot,obj.nfreeGroupShoot+obj.Ns_,nfreeGroupShoot_);
            nfreeGroupShoot = nfreeGroupShoot_;
        }
        else {
            nfreeGroupShoot = NULL;
        }
    }

    return *this;
}

/**
* Initialize members
*/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::init(){

	int i,n,ibgx,indx,ibgu,indu;
	double *xmin,*xmax,*umin,*umax;
    double *xmaxT;
    double *zl,*zu,*zzM;
    double *hp,*hs,*hhs,*hhp;

    /* Initialize shooter */
    shooter.init();
    shooter.displayShootingParameters();
#if SCAL_VAR
    /* Compute state & input scaling transformations */
    shooter.computeStateScaling();
    shooter.computeInputScaling();
#endif

    /* Get stage dimensions */
    du = shooter.getDu();
    dubl = du*SZDBL;
    dx = shooter.getDx();
    dxbl = dx*SZDBL;
    dpc = shooter.getDpc();

    /* Shooting group dimensions */
    dxu = dx+du;
    dxux = dxu+dx;
    
    /* Initialize dimensions */
    Ns_ = shooter.getNs_();
    Ns = shooter.getNs();
	dZshoot = shooter.getDprimShoot();
    dZshootBL = dZshoot*SZDBL;
	dZpath = shooter.getDprimPath();
    dZpathBL = dZpath*SZDBL;
    dZ = shooter.getDprim();
    dZBL = dZ*SZDBL;
    dZNT = dZ*SZINT;
    dLshoot = shooter.getDduaShoot();
    dL = shooter.getDdua();
    dLBL = dL*SZDBL;
    dGshoot = Ns_*dxux;
    dGshootBL = dGshoot*SZDBL;
    dGshootNT = dGshoot*SZINT;
    dGpath = Ns_*dpc;
    dGpathBL = dGpath*SZDBL;
    dGpathNT = dGpath*SZINT;
    dhS = dxux*dxux;
    dhP = dpc*dpc;
    dHshoot = Ns_*dhS;
    dHpath = Ns_*dhP;

    /* Allocate primal optimizer */
    if (z != NULL)
        delete [] z;
    z = new double[dZ];
    memset(z,0,dZBL);

	/* Allocate candidate for next iterate */
    if (zS != NULL)
        delete [] zS;
    zS = new double[dZ];
    memset(zS,0,dZBL);
    
    /* Allocate difference between iterates */
    if (negstp != NULL)
        delete [] negstp;
    negstp = new double[dZ];
    memset(negstp,0,dZBL);

	/* Allocate dual optimizer */
    if (mu != NULL)
        delete [] mu;
    mu = new double[dL];
    
    /* Allocate equality constraints */
    if (econ != NULL)
        delete [] econ;
    econ = new double[dL];
    
    /* Allocate optimal input */
    if (uopt != NULL)
        delete [] uopt;
    uopt = new double[du];

#ifdef BIOREC
    /* Allocate auxiliary vector of dimension dxu (primal) */
    if (uxAux != NULL)
        delete [] uxAux;
    uxAux = new double[dxu];

    /* Allocate auxiliary vector of dimension dx (dual) */
    if (xAux != NULL)
        delete [] xAux;
    xAux = new double[dx];

    /* Allocate auxiliary vector of dimension dpc (primal) */
    if (pcAux != NULL)
        delete [] pcAux;
    pcAux = new double[dpc];
    
    /* Allocate auxiliary vector of dimension dpc (dual) */
    if (pcAuxx != NULL)
        delete [] pcAuxx;
    pcAuxx = new double[dpc];
#endif

    /* Allocate primal bounds (fixed and varying) */
    if (zmin != NULL)
        delete [] zmin;
    zmin = new double[dZ];

    if (zmax != NULL)
        delete [] zmax;
    zmax = new double[dZ];
    
    /* Initialize bounds on slacks for path-constraints */    
    memset(zmin+dZshoot,0,dZpath*SZDBL);
    zzM = zmax+dZshoot;
    while (zzM != zmax+dZ)
        *zzM++ = POSINF;

    /* Initialize bounds zmin & zmax */
	umin = shooter.getUmin();	
	umax = shooter.getUmax();
	xmin = shooter.getXmin();
	xmax = shooter.getXmax();
    xmaxT = shooter.getXmaxT();
	zl = zmin;
    zu = zmax;    
    for (n=0;n<Ns_;++n){
    	std::copy(xmin,xmin+dx,zl);
    	std::copy(xmax,xmax+dx,zu);
    	zl += dx;
    	zu += dx;
    	std::copy(umin,umin+du,zl);
    	std::copy(umax,umax+du,zu);
    	zl += du;
    	zu += du;
    }
    std::copy(xmin,xmin+dx,zl);
//    std::copy(xmax,xmax+dx,zu);
	std::copy(xmaxT,xmaxT+dx,zu);    

    /* Initialize AL gradients */
    if (gAL != NULL)
        delete [] gAL;
    gAL = new double[dZ];

    gALpath = gAL+dZshoot;

    if (gALpathS != NULL)
        delete [] gALpathS;
    gALpathS = new double[dGpath];
    
    if (gALshoot != NULL)
        delete [] gALshoot;
    gALshoot = new double[dGshoot];

    if (gALshootS != NULL)
        delete [] gALshootS;
    gALshootS = new double[dGshoot];
    
    /* Initialize differences of AL gradients */
    if (dffGshoot != NULL)
        delete [] dffGshoot;
    dffGshoot = new double[dGshoot];

    if (dffGpath != NULL)
        delete [] dffGpath;
    dffGpath = new double[dGpath];
    
    /* Initialize AL hessian estimates (diagonal scaling matrices for every stage) */
    if (hALshoot != NULL)
        delete [] hALshoot;
    hALshoot = new double[dHshoot];
    memset(hALshoot,0,dHshoot*SZDBL);
    hs = hALshoot;
    for (n=0;n<Ns_;++n){
        hhs = hs;
        for (i=0;i<dxux;++i){
            *hhs = sr_shoot_ini;
            hhs += dxux+1;
        }
        hs += dhS;
    }

    if (hALpath != NULL)
        delete [] hALpath;
    hALpath = new double[dHpath];
    memset(hALpath,0,dHpath*SZDBL);
    hp = hALpath;
    for (n=0;n<Ns_;++n){
        hhp = hp;
        for (i=0;i<dpc;++i){
            *hhp = sr_path_ini;
            hhp += dpc+1;
        }
        hp += dhP;
    }

    /* Initialize reduced hessian approximation */
    if (hALshootRed != NULL)
        delete [] hALshootRed;
    hALshootRed = new double[dHshoot];

    if (hALpathRed != NULL)
        delete [] hALpathRed;
    hALpathRed = new double[dHpath];

    /* Initialize prdHss */
    if (prdHss != NULL)
        delete [] prdHss;
    prdHss = new double[dZ];
    memset(prdHss,0,dZBL);

    /* Allocate secant vectors */
    if (secVecShoot != NULL)
        delete [] secVecShoot;
    secVecShoot = new double[dxux];

    if (secVecPath != NULL)
        delete [] secVecPath;
    secVecPath = new double[dpc];

    /* Initialize vector of free/active variables */
    if (freeVars != NULL)
        delete [] freeVars;
    freeVars = new int[dZ];
    memset(freeVars,0,dZNT);

    /* Initialize vector of indices of free/active variables */
    if (indFreeVars != NULL)
        delete [] indFreeVars;
    indFreeVars = new int[dZ];
    for (i=0;i<dZ;++i)
        indFreeVars[i] = -1;

    /* Initialize status vector (-1: var fixed to lower bound, 0: var free, 1: var fixed to upper bound)*/
    if (varStatus != NULL)
        delete [] varStatus;
    varStatus = new int[dZ];
    memset(varStatus,0,dZNT);

    /* Initialize arrays of free variables per state, input, slack and shooting group */
    if (ivarFreeState != NULL)
        delete [] ivarFreeState;
    ivarFreeState = new int[Ns*dx];
    for (i=0;i<dx*Ns;++i)
        ivarFreeState[i] = -1;

    if (ivarFreeInput != NULL)
        delete [] ivarFreeInput;
    ivarFreeInput = new int[Ns_*du]; 
    for (i=0;i<du*Ns_;++i)
        ivarFreeInput[i] = -1;

    if (ivarFreeSlack != NULL)
        delete [] ivarFreeSlack;
    ivarFreeSlack = new int[Ns_*dpc];
    for (i=0;i<dpc*Ns_;++i)
        ivarFreeSlack[i] = -1;

    if (ivarFreeGroupShoot != NULL)
        delete [] ivarFreeGroupShoot;
    ivarFreeGroupShoot = new int[Ns_*dxux];
    for (i=0;i<dxux*Ns_;++i)
        ivarFreeGroupShoot[i] = -1;

    if (nfreeState != NULL)
        delete [] nfreeState;
    nfreeState = new int[Ns];
    for (i=0;i<Ns;++i)
        nfreeState[i] = -1;

    if (nfreeInput != NULL)
        delete [] nfreeInput;
    nfreeInput = new int[Ns_];
    for (i=0;i<Ns_;++i)
        nfreeInput[i] = -1;

    if (nfreeSlack != NULL)
        delete [] nfreeSlack;
    nfreeSlack = new int[Ns_];
    for (i=0;i<Ns_;++i)
        nfreeSlack[i] = -1;

    if (nfreeGroupShoot != NULL)
        delete [] nfreeGroupShoot;
    nfreeGroupShoot = new int[Ns_];
    for (i=0;i<Ns_;++i)
        nfreeGroupShoot[i] = -1;

    /**
    ** Initialize local quadratic model qmod 
    **/
    qmod.setDimensions(dZ,dZshoot,dxu,dxux,dpc);
    qmod.setNumberShoots(Ns_);
	qmod.setIterate(z);
    qmod.setCauchy(zS);
	qmod.setCandidate(zS);
	qmod.setGradient(gAL);	
	qmod.setHessianShoot(hALshoot);
	qmod.setHessianPath(hALpath);
    qmod.setProductHessian(prdHss);

    /**
    ** Initialize pointers for shooter, activDetector, precSubRefiner and precProjRefiner 
    **/
    /** Shooter initialization (members must have been allocated in Simulator) */
    shooter.setPrimal(z);
    shooter.setCandidate(zS);
    shooter.setDual(mu);
    shooter.setEcon(econ);
    shooter.setGradient(gAL);
    shooter.setGradientShooting(gALshoot);
    shooter.setGradientShootingCandidate(gALshootS);
    shooter.setGradientPath(gALpath);
    shooter.setGradientPathCandidate(gALpathS);
    shooter.setDiffGradientPath(dffGpath);
#if REFIN_PCG_RES
    /** Initialization of preconditoned subspace refiner */
    precSubRefiner.setQuadraticModel(qmod);
    precSubRefiner.setDimension(dZ);
    precSubRefiner.setDimensionShoot(dx,du,dZshoot,dHshoot,dhS,dGshoot);
    precSubRefiner.setDimensionPath(dpc,dZpath,dHpath,dhP,dGpath);
    precSubRefiner.setGradient(gAL);
    precSubRefiner.setProductHessian(prdHss);
    precSubRefiner.setHessians(hALshoot,hALpath);
    precSubRefiner.setReducedHessians(hALshootRed,hALpathRed);
    precSubRefiner.setIterate(z);
    precSubRefiner.setCandidate(zS);
    precSubRefiner.setNLPbounds(zmin,zmax);
    precSubRefiner.setNumIntervals(shooter.getNs_());
    precSubRefiner.setRegCoefPrec(pr_reg_coef);
    precSubRefiner.setFreeVars(freeVars);
    precSubRefiner.setIndFreeVars(indFreeVars);
    precSubRefiner.setVarStatus(varStatus);
    precSubRefiner.setIndexFreeVars(ivarFreeState,ivarFreeInput,ivarFreeSlack,ivarFreeGroupShoot);
    precSubRefiner.setNumberFreeVars(nfreeState,nfreeInput,nfreeSlack,nfreeGroupShoot);
    precSubRefiner.init();
#endif
#if REFIN_MOR_TOR
    /** Initialize preconditioned projected subspace refiner */
    precProjRefiner.setQuadraticModel(qmod);
    precProjRefiner.setDimension(dZ);
    precProjRefiner.setDimensionShoot(dx,du,dZshoot,dHshoot,dhS,dGshoot);
    precProjRefiner.setDimensionPath(dpc,dZpath,dHpath,dhP,dGpath);
    precProjRefiner.setGradient(gAL);
    precProjRefiner.setProductHessian(prdHss);
    precProjRefiner.setHessians(hALshoot,hALpath);
    precProjRefiner.setReducedHessians(hALshootRed,hALpathRed);
    precProjRefiner.setIterate(z);
    precProjRefiner.setMinorIterate(zS);
    precProjRefiner.setNLPbounds(zmin,zmax);
    precProjRefiner.setNumIntervals(shooter.getNs_());
    precProjRefiner.setRegCoefPrec(pr_reg_coef);
    precProjRefiner.setFreeVars(freeVars);
    precProjRefiner.setIndFreeVars(indFreeVars);
    precProjRefiner.setVarStatus(varStatus);
    precProjRefiner.setIndexFreeVars(ivarFreeState,ivarFreeInput,ivarFreeSlack,ivarFreeGroupShoot);
    precProjRefiner.setNumberFreeVars(nfreeState,nfreeInput,nfreeSlack,nfreeGroupShoot);
    precProjRefiner.init();
#endif
    /** Activity detector initialization */
    activDetector.setQuadraticModel(qmod);
    activDetector.setNs_(Ns_);
    activDetector.setDimensions(dZ,dZshoot,dGshoot);
    activDetector.setStateDimension(dx);
    activDetector.setInputDimension(du);
    activDetector.setGroupDimensions(dxu,dxux,dpc);
    activDetector.setIterate(z);
    activDetector.setCauchy(zS);
    activDetector.setGradient(gAL);
    activDetector.setProductHessian(prdHss);
    activDetector.setNLPbounds(zmin,zmax);
#if REFIN_PCG_RES
    activDetector.setReducedResidual(precSubRefiner.getReducedResidual());
#endif
#if REFIN_MOR_TOR
    activDetector.setReducedResidual(precProjRefiner.getReducedResidual());
#endif
    activDetector.setFreeVars(freeVars);
    activDetector.setIndFreeVars(indFreeVars);
    activDetector.setVarStatus(varStatus);
    activDetector.setIndexFreeVars(ivarFreeState,ivarFreeInput,ivarFreeSlack,ivarFreeGroupShoot);
    activDetector.setNumberFreeVars(nfreeState,nfreeInput,nfreeSlack,nfreeGroupShoot);
    activDetector.allocateDebug();
}

/**
** Set bounds on terminal state equal to xest
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::setPeriodicTerminal(){

/*    int i,izTerm;
    double *xest,*zm,*zM;

    izTerm = shooter.getIndexLastNode();
    xest = shooter.getXest();

    zm = zmin+izTerm;
    zM = zmax+izTerm;
    zm[0] = xest[0]-1.;
    zM[0] = xest[0]+1.;
    zm[1] = xest[1]-1.;
    zM[1] = xest[1]+1.;
    zm[2] = xest[2]-1.;
    zM[2] = xest[2]+1.; */
}

/**
** Cold-start primal-dual optimizer
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::coldStartPrimalDual(){
    
//    std::cout<<"00000000 >> COLD-START"<<std::endl;

    /* Cold start primal iterate z */
/*#if UNSCAL_VAR
    memset(z,0,dZBL);
#endif
#if SCAL_VAR
    for (int i=0;i<dZ;++i) 
        z[i] = 0.5;
#endif */
    double *zz;
    double *zzm,*zzM;
/*#if DPER
    int i,n,izTerm;
    double *sini,*sfin;
    double *xest;
#endif */
    zz = z;
    zzm = zmin;
    zzM = zmax;
    while (zz != z+dZ){
        if ((*zzm<=NEGINF)||(*zzM>=POSINF)){
            *zz++ = 0.;
        }
        else {
            *zz++ = 0.5*(*zzm+*zzM);
        }
        zzm++;
        zzM++;
    }
//#if DPER
    /* Init. first & final shooting nodes enforcing periodicity */
/*    izTerm = shooter.getIndexLastNode();
    xest = shooter.getXest();
    sini = z;
    sfin = sini+izTerm;
    for (i=0;i<dx;++i){
        sini[i] = 0.; //xest[i];
        sfin[i] = 0.; //xest[i];
    } */
    /* Init. rest of shooting nodes and inputs */
//    zz = z+dxu;
//    zz = z;
//    while (zz != z+izTerm){
        /* Node */
/*        *zz++ = 5.8;
        *zz++ = 19.;
        *zz++ = 20.5;
        *zz++ = 0.;
        *zz++ = 0.; */
        /* Input */
/*        *zz++ = 33.;
    } 
#endif
*/

    /* Cold-start primal candidate zS */
    memset(zS,0,dZBL);

    /* Cold-start Lagrange multiplier */
    memset(mu,0,dLBL);
}

/**
** Bioreactor: initialise on economically optimal trajectory 
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::initBioReactor(){

    /* Primal */
    z[0] = 5.479699008;
    z[1] = 17.18171851;
    z[2] = 19.28396712;
    z[3] = -7.881400456e-11;
    z[4] = -5.689121156e-10;
    z[5] = 29.84158637;
    z[6] = 5.607988223;
    z[7] = 16.54676905;
    z[8] = 19.61674354;
    z[9] = 1.492079319;
z[10] = 0.2771502613;
z[11] = 29.98179851;
z[12] = 5.741442557;
z[13] = 16.03629316;
z[14] = 19.9987356;
z[15] = 2.991169244;
z[16] = 0.5608833314;
z[17] = 29.8288396;
z[18] = 5.874122422;
z[19] = 15.53497323;
z[20] = 20.40602318;
z[21] = 4.482611224;
z[22] = 0.8512813084;
z[23] = 29.65393399;
z[24] = 6.003741315;
z[25] = 15.03970785;
z[26] = 20.82458135;
z[27] = 5.965307924;
z[28] = 1.14824439;
z[29] = 29.52052222;
z[30] = 6.128340379;
z[31] = 14.56842171;
z[32] = 21.24301239;
z[33] = 7.441334034;
z[34] = 1.451571468;
z[35] = 30.67406054;
z[36] = 6.23530053;
z[37] = 14.53827928;
z[38] = 21.62619742;
z[39] = 8.975037061;
z[40] = 1.760778467;
z[41] = 33.31942926;
    z[42] = 6.294542966;
    z[43] = 15.35545389;
    z[44] = 21.89384621;
    z[45] = 10.64100852;
    z[46] = 2.074293107;
    z[47] = 33.44307774;
    z[48] = 6.299080147;
    z[49] = 16.05494844;
    z[50] = 22.01118844;
    z[51] = 12.31316241;
    z[52] = 2.389323263;
    z[53] = 35.5344222;
    z[54] = 6.250633173;
    z[55] = 17.30425979;
    z[56] = 21.97045169;
    z[57] = 14.08988352;
    z[58] = 2.703309665;
    z[59] = 36.08450532;
    z[60] = 6.153251348;
    z[61] = 18.50113316;
    z[62] = 21.77155043;
    z[63] = 15.89410879;
    z[64] = 3.013573457;
    z[65] = 35.3457809;
    z[66] = 6.032244747;
    z[67] = 19.24549784;
    z[68] = 21.47071073;
    z[69] = 17.66139783;
    z[70] = 3.318249937;
    z[71] = 36.88782857;
    z[72] = 5.896441424;
    z[73] = 20.35946786;
    z[74] = 21.09617258;
    z[75] = 19.50578926;
    z[76] = 3.616546561;
    z[77] = 37.47287845;
    z[78] = 5.745394538;
    z[79] = 21.4542836;
    z[80] = 20.6527613;
    z[81] = 21.37943318;
    z[82] = 3.907641287;
    z[83] = 37.72386005;
    z[84] = 5.587383205;
    z[85] = 22.42556323;
    z[86] = 20.16494597;
    z[87] = 23.26562619;
    z[88] = 4.190972408;
    z[89] = 36.63629576;
    z[90] = 5.437319247;
    z[91] = 22.87396783;
    z[92] = 19.67433909;
    z[93] = 25.09744097;
    z[94] = 4.466519209;
    z[95] = 34.12845093;
    z[96] = 5.318597675;
    z[97] = 22.46369858;
    z[98] = 19.24651628;
    z[99] = 26.80386352;
    z[100] = 4.73523604;
    z[101] = 29.79704633;
    z[102] = 5.263588645;
    z[103] = 20.79843604;
    z[104] = 18.97371612;
    z[105] = 28.29371584;
    z[106] = 4.999450913;
    z[107] = 28.7135592;
    z[108] = 5.285036883;
    z[109] = 19.15948931;
    z[110] = 18.90949496;
    z[111] = 29.7293938;
    z[112] = 5.262867472;
    z[113] = 28.92430806;
    z[114] = 5.365878358;
    z[115] = 17.91505089;
    z[116] = 19.03175202;
    z[117] = 31.1756092;
    z[118] = 5.528940529;
    z[119] = 29.84529621;
    z[120] = 5.479699007;
    z[121] = 17.18171851;
    z[122] = 19.28396712;
    z[123] = 32.66787401;
    z[124] = 5.8;

    /* Optimal input */
    uopt[0] = 29.84158637;

    /* Dual */
mu[0] = 0.;
mu[1] = 0.;
mu[2] = 0.;
mu[3] = 5.447808693e-09;
mu[4] = -0.529999452;
mu[5] = -0.04583317938;
mu[6] = -1.292885798e-08;
mu[7] = 0.02083328138;
mu[8] = 1.150962659e-08;
mu[9] = -0.5299994597;
mu[10] = -0.04583319329;
mu[11] = -1.548210449e-08;
mu[12] = 0.02083332379;
mu[13] = 1.55190083e-08;
mu[14] = -0.5299994669;
mu[15] = -0.04583322981;
mu[16] = -1.651256909e-08;
mu[17] = 0.02083333059;
mu[18] = 1.6809798e-08;
mu[19] = -0.52999948;
mu[20] = -0.04583328376;
mu[21] = -1.326503352e-08;
mu[22] = 0.02083332338;
mu[23] = 1.713971187e-08;
mu[24] = -0.5299994826;
mu[25] = -0.04583332951;
mu[26] = -2.496207685e-08;
mu[27] = 0.02083332272;
mu[28] = 4.759437289e-09;
mu[29] = -0.5299994809;
mu[30] = -0.04583335007;
mu[31] = -2.594635617e-09;
mu[32] = 0.0208333147;
mu[33] = 2.493116824e-10;
mu[34] = -0.5299994958;
mu[35] = -0.04583332465;
mu[36] = 1.305373587e-08;
mu[37] = 0.02083331643;
mu[38] = -3.845013197e-09;
mu[39] = -0.5299995022;
mu[40] = -0.0458332867;
mu[41] = 2.19456453e-08;
mu[42] = 0.02083331639;
mu[43] = -7.426859128e-09;
mu[44] = -0.5299995052;
mu[45] = -0.04583324964;
mu[46] = 1.616911049e-08;
mu[47] = 0.02083332196;
mu[48] = -9.102230081e-09;
mu[49] = -0.5299995081;
mu[50] = -0.04583323797;
mu[51] = -7.112532785e-10;
mu[52] = 0.02083333308;
mu[53] = -7.854605855e-09;
mu[54] = -0.5299995098;
mu[55] = -0.04583325426;
mu[56] = -1.426077034e-08;
mu[57] = 0.02083334327;
mu[58] = -5.031530748e-09;
mu[59] = -0.5299995077;
mu[60] = -0.04583331054;
mu[61] = -2.882831751e-08;
mu[62] = 0.02083335557;
mu[63] = -5.393197e-09;
mu[64] = -0.5299994958;
mu[65] = -0.04583338589;
mu[66] = -2.001776522e-08;
mu[67] = 0.020833361;
mu[68] = -5.922373703e-09;
mu[69] = -0.5299994858;
mu[70] = -0.04583344297;
mu[71] = 5.170619488e-09;
mu[72] = 0.020833372;
mu[73] = -5.363887112e-09;
mu[74] = -0.5299994826;
mu[75] = -0.04583349746;
mu[76] = 2.329922921e-08;
mu[77] = 0.02083340049;
mu[78] = -5.1034732e-09;
mu[79] = -0.5299994866;
mu[80] = -0.04583361726;
mu[81] = 1.75873538e-08;
mu[82] = 0.02083343692;
mu[83] = -3.439915019e-09;
mu[84] = -0.5299994928;
mu[85] = -0.04583379395;
mu[86] = 2.101430141e-09;
mu[87] = 0.02083345961;
mu[88] = 1.304201191e-09;
mu[89] = -0.5299994993;
mu[90] = -0.04583398266;
mu[91] = -2.129070253e-08;
mu[92] = 0.02083346923;
mu[93] = 5.187494878e-09;
mu[94] = -0.5299994906;
mu[95] = -0.04583416503;
mu[96] = -4.721538716e-08;
mu[97] = 0.02083349318;
mu[98] = -5.224798372e-09;
mu[99] = -0.5299994685;
mu[100] = -0.04583447427;
mu[101] = -2.386038034e-08;
mu[102] = 0.02083352355;
mu[103] = 3.491607004e-09;
mu[104] = -0.5299995238;
mu[105] = -0.04583418677;
mu[106] = 1.173557251e-07;
mu[107] = 0.02083374279; 
}

/**
** Shift primal & dual variables in an NMPC-like fashion
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::shiftPrimalDual(){
    

//#ifdef BIOREC
    /* Copy first input q_0 and state s_1  */
//    std::copy(z+dx,z+dxu,uxAux);
//#if DPC
    /* Copy slack part of primal sequence */
//#endif
    /* Copy first dual */
//    std::copy(mu,mu+dx,xAux);
//#if DPC
    /* Copy first path dual*/
//#endif
//#endif


    /* Shift shooting-part of primal sequence */
    std::copy(z+dxu,z+dxu*(Ns_-1)+dx,z);
#if DPC
    /* Shift slacks-part of primal sequence */  
    std::copy(z+dZshoot+dpc,z+dZ,z+dZshoot);
#endif 
    /* Shift shooting-part of dual sequence */
    std::copy(mu+dx,mu+dx+dLshoot,mu);
#if DPC
    /* Shift slacks-part of dual sequence */
    std::copy(mu+dx+dLshoot+dpc,mu+dL,mu+dx+dLshoot);
#endif


//#ifdef BIOREC
    /* Put first at last position */
//    std::copy(uxAux,uxAux+dxu,z+Ns_*dxu-du);
//#if DPC
    /* Put first slack at last position */
//#endif
    /* Put first shooting dual at last position */
//    std::copy(xAux,xAux+dx,mu+dLshoot);
//#if DPC
    /* Put first path dual at last position */
//#endif
//#endif
}

/**
** Set initial hessian approximations
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::initHessianApprox(){

    int i,n;
    double *hs,*hhs,*hp,*hhp;

    /* Set AL hessian of shooting part to diagonal scaling matrix */
    memset(hALshoot,0,dHshoot*SZDBL);
    hs = hALshoot;
    for (n=0;n<Ns_;++n){
        hhs = hs;
        for (i=0;i<dxux;++i){
            *hhs = sr_shoot_ini;
//            *hhs = 2.*halfRhoIni;
            hhs += dxux+1;
        }
        hs += dhS;
    }
#if DPC
    /* Set AL hessian of slacks part to diagonal scaling matrix */
    memset(hALpath,0,dHpath*SZDBL);
    hp = hALpath;
    for (n=0;n<Ns_;++n){
        hhp = hp;
        for (i=0;i<dpc;++i){
            *hhp = sr_path_ini;
//            *hhp = halfRhoIni;
            hhp += dpc+1;
        }
        hp += dhP;
    }
#endif    
}

/**
** Evaluate primal KKT satisfaction ||P(x-g(x))-x||
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::evalKKTprim(){

	int k;
    double dd,az,al,au,ag;

    kkt2 = 0.;
    for (k=0;k<dZ;++k){
        az = *(z+k);
        ag = *(gAL+k);
        al = *(zmin+k);
        au = *(zmax+k);
        dd = PROJ(az-ag,al,au)-az;
        kkt2 += dd*dd;
    /*    dd = az-ag;
        if (dd>au){
            dd = au-az;
            kkt2 += dd*dd;
            continue;
        }
        if (dd<al){
            dd = al-az;
            kkt2 += dd*dd;
            continue;
        } 
        dd -= az; */
    }
    kkt = sqrt(kkt2);
}

/** 
** Evaluate 2-norm of equality constraints
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::evalEconNorm(){

    double aux,*e,*end;

    e = econ;
    nec2 = 0.;
    while (e != econ+dL){
        aux = *e++; 
        nec2 += aux*aux;
    }
    nec = sqrt(nec2);
}    

/**
** Compute AL gradient at current primal-dual iterate and initial penalty parameter
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::evalALgradient(){

    double val;
    
#ifndef PEROCP
    val = shooter.shootStatesFromIterate(halfRhoIni);
    shooter.shootAdjointsFromIterate();
#else
    std::cout<<"ERROR: no gradient evaluation for periodic OCP."<<std::endl;
#endif
}

/**
** Assemble reduced hessian approximation in a blockwise manner
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::assembleReducedHessian(){

    int i,j,n,nfree,irw;
    int *nfgS,*nfP,*ivfgS,*ivfP;
    double *hs,*hp,*hSr,*hPr;

    nfgS = nfreeGroupShoot;
    ivfgS = ivarFreeGroupShoot;
    hs = hALshoot;
    hSr = hALshootRed;
#if DPC
    nfP = nfreeSlack;
    ivfP = ivarFreeSlack;
    hp = hALpath;
    hPr = hALpathRed;
#endif
    for (n=0;n<Ns_;++n){    
        /* Shooting group */
        nfree = *(nfgS+n);
        for (i=0;i<nfree;++i){
            irw = *(ivfgS+i)*dxux;
            for (j=0;j<nfree;++j)
                *hSr++ = *(hs+irw+*(ivfgS+j));
        }
        hs += dhS;
        ivfgS += dxux;
#if DPC
        /* Slack group */
        nfree = *(nfP+n);
        for (i=0;i<nfree;++i){
            irw = *(ivfP+i)*dpc;
            for (j=0;j<nfree;++j)
                *hPr++ = *(hp+irw+*(ivfP+j));
        }   
        hp += dhP;
        ivfP += nfree;
#endif
    }
}

/**
** Full primal loop: trust-region SR1 on augmented Lagrangian.
** Stopped if KKT tolerance or max. # iterations hit. 
** Input: - kkt tolerance from Lancelot's outer loop
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::primalFull(const double kktLan){

	int it,i,j,n,testIterate,testDecrease;
    double vmod,vcan,val,actRed,preRed;
    double aux,ssz,alpha;
    double auxz,auxg;
    double srtio,srtio2,ndz,ndz2,ndg,ndg2;
    double n2r,n2g,gts;
    double *zz;
	double *zm,*zM;
	double *pz,*dz,*pzs;
	double *dzs,*dzp,*ddzs,*ddzp,*dddzs,*dddzp;
	double *dgp,*pgp,*pgpS,*ddgp,*dddgp;
	double *dgs,*pgs,*pgsS,*ddgs,*dddgs;
    double *g;
	double *phs,*pphs,*php,*pphp;

    /**
    ** Initialization phase  
    **/
//    std::cout<<"Start shooting from iterate..."<<std::endl;
    /* Evaluate augmented Lagrangian (objective+constraints) val */
//    val = shooter.shootStatesFromIterate(halfRho);
    /* Evaluate full gradient of augmented Lagrangian gAL & gALshoot */
//    shooter.shootAdjointsFromIterate();
#ifdef PEROCP
//    std::cout<<"Periodic shooting from iterate"<<std::endl;
    val = shooter.shootFromIteratePeriodic(halfRho);
#else
    val = shooter.shootFromIterate(halfRho);
#endif
//    std::cout<<"...done."<<std::endl;

	/* Get KKT satisfaction */
    evalKKTprim();

    /* Initialize squared tolerance for refinement phase kkt2*min(0.01,...)
    	- kkt for local superlinear convergence 
    	- kkt2 for quadratic local convergence */
#if QUAD_CV
	cgTol2 = kkt2*MIN(0.01,kkt2);
#endif
#if SUPL_CV
    cgTol2 = kkt2*MIN(0.01,kkt);
#endif

	/* Initialize trust-region radius: coef*kkt */
	trRad = 0.1*kkt;

    /**
    ** Trust-region iterations 
    **/
	it = 0;
    numPrimIt = 0;
#if PRINT_TR
    std::cout<<"val initial= "<<std::setprecision(10)<<val<<std::endl;
    std::cout<<"kkt initial= "<<std::setprecision(10)<<kkt<<std::endl;
#endif
#if COUNT_CG
    cumCG = 0;
    cumRestart = 0;
#endif 
	/* Set initial step-size */
    ssz = sini;
    while (1){
#if PRINT_TR
        std::cout<<"T-R ITERATION "<<it+1<<std::endl;
        std::cout<<"T-R radius= "<<trRad<<std::endl;
        std::cout<<std::setprecision(5)<<"val= "<<val<<std::endl;
        std::cout<<std::setprecision(5)<<"kkt= "<<kkt<<", kktLan= "<<kktLan<<std::endl;
//        displayMembers();
//        std::cout<<"Running step-size= "<<ssz<<std::endl;
#endif
        /* Check KKT and # iterations */
        if ((kkt<=kktLan)||(it>=maxPit)){
#if PRINT_TR
            std::cout<<"Break after "<<it<<" iterations with kkt: "<<kkt<<std::endl;
            std::cout<<"val final: "<<std::setprecision(10)<<val<<std::endl;
            std::cout<<"Max. primal iters: "<<maxPit<<std::endl;
#endif
        	break;
        }
        /* Check size of trust region radius */
        if (trRad<=trRadTol){
#if PRINT_TR
            std::cout<<"Trust-region radius too small, abort."<<std::endl;
            std::cout<<"trRad="<<trRad<<"<= trRadTol="<<trRadTol<<std::endl;
#endif 
            break;
        }   
		/**
        ** Generate candidate for next iterate inside trust-region 
        **/
        /** 1. Cauchy phase */
#if PRINT_TR
        std::cout<<"Start Cauchy phase..."<<std::endl;
#endif
#if CCH_PRJ_SRCH
        /* Use backtracking projected search */
//        std::cout<<"!!!!! Initial step-size= !!!!"<<ssz<<std::endl;
        activDetector.findProjectedSearch(&ssz,&vmod,&norIstp,trRad);
        if (activDetector.getError()){
//            activDetector.checkDescent();
            errPrim = 1;
            std::cout<<"ERROR IN CAUCHY, interruption at TR iter "<<it<<" with kkt="<<kkt<<std::endl;
            return;
        }
//        ssz = sini;
#endif
#if CCH_INTER_EXTRA
        /* Use interpolation-extrapolation strategy (see TRON) */
        activDetector.findInterpExtrap(&ssz,&vmod,&norIstp,trRad);
        if (activDetector.getError()){
            errPrim = 1;
            std::cout<<"ERROR IN CAUCHY, interruption at TR iter "<<it<<" with kkt="<<kkt<<std::endl;
            return;
        }
//        ssz = sini;
#endif
        /** 2. Refinement phase */
#if PRINT_TR            
        std::cout<<"Start refinement phase..."<<std::endl;
#endif
#if REFIN_PCG_RES
        /* PCG on active face at Cauchy point with restart if problem bound hit during search */
        precSubRefiner.findCandidateRestart(cgTol2,trRad);
#if COUNT_CG
//        std::cout<<"Add "<<precSubRefiner.getCumCG()<<" PCG iters."<<std::endl;
        cumCG += precSubRefiner.getCumCG();
        cumRestart += precSubRefiner.getCumRestart();
#endif
#endif
#if REFIN_MOR_TOR
        /* Projected searches along directions generated by PCG */
        precProjRefiner.findCandidate(vmod,cgTol2,trRad);
#if COUNT_CG
//        std::cout<<"Add "<<precProjRefiner.getCumCG()<<" PCG iters."<<std::endl;
        cumCG += precProjRefiner.getCumCG();
#endif
#endif
#if PRINT_TR
        std::cout<<"TRUST-REGION: kkt satisf.= "<<kkt<<std::endl;
#if COUNT_CG
        std::cout<<"TRUST-REGION: cumul. CG iters= "<<cumCG<<std::endl;
        std::cout<<"TRUST-REGION: cumul. CG restarts= "<<cumRestart<<std::endl;
#endif
#endif
        /* Update negative difference z-zS, compute its infinity norm (and euclidean norm) */   
        dz = negstp;
        zz = z;
        pzs = zS;
        g = gAL;
        norIstp = -1.;
        gts = 0.;
        while (zz != z+dZ){
            aux = *zz++-*pzs++;
            *dz++ = aux;
            gts += *g++*aux;
            norIstp = MAX(norIstp,ABS(aux));
        }

/*        if (norIstp<=diffTol){
#if PRINT_TR
            std::cout<<"Primal progress is too low, norIstp= "<<norIstp<<
                       " after "<<it<<" iters, break with kkt= "<<kkt<<", nec= "<<nec<<std::endl;
            std::cout<<"Max. primal iters: "<<maxPit<<std::endl;
#endif
            break;
        } */

        /* Evaluate model at candidate */
#if REFIN_PCG_RES
        vmod = val+qmod.evalModelCandidate();        
#endif
#if REFIN_MOR_TOR
        vmod = val+precProjRefiner.getValMod();
#endif       
        /* Evaluate AL and separated gradients at candidate */
#ifdef PEROCP
//        std::cout<<"Periodic shooting from candidate"<<std::endl;
        vcan = shooter.shootFromCandidatePeriodic(halfRho);
#else
        vcan = shooter.shootFromCandidate(halfRho);
#endif
        /* Compute shooting AL gradient difference gALshoot(zS)-gALshoot(z) */
        dgs = dffGshoot;
        pgs = gALshoot; 
        pgsS = gALshootS;
        while (pgs != gALshoot+dGshoot)
            *dgs++ = *pgsS++-*pgs++;
#if DPC
        /* Compute path-slacks AL gradient difference gALpath(zS)-gALpath(z) */
        dgp = dffGpath;
        pgp = gALpath;
        pgpS = gALpathS;
        while (pgp != gALpath+dGpath)
        	*dgp++ = *pgpS++-*pgp++;
#endif
        /** 
        ** SR1 updates on separated blocks 
        **/
        phs = hALshoot;    
        dgs = dffGshoot;
        dzs = negstp;
#if DPC
        php = hALpath;
        dgp = dffGpath;
        dzp = dzs+dZshoot;
#endif
        /* Loop over shooting blocks */
        for (n=0;n<Ns_;++n){ 
            /* On shooting AL hessian */
            /* dgs<-dgs+phs*dzs and srtio<-(dgs+phs*dzs)'*dzs */
            ddgs = dgs;
            pphs = phs;
            dddzs = dzs;
            srtio = 0.;
            ndz2 = 0.;
            ndg2 = 0.;
            for (i=0;i<dxux;++i){
                aux = 0.;
                ddzs = dzs;
                for (j=0;j<dxux;++j){
                    aux += *pphs**ddzs;
                    ++pphs;
                    ++ddzs;
                }
                *ddgs += aux;
                auxz = *dddzs;
                auxg = *ddgs;
                srtio += auxg*auxz;
                ndz2 += auxz*auxz;
                ndg2 += auxg*auxg;
                ++ddgs;
                ++dddzs;
            } 
            srtio2 = srtio*srtio;
            /* Check denominator and if zero skip update on current shooting-block */
            /* !!!! Change SR1 update criterion !!!! */
            //if (ABS(srtio)>=POSZER){
            if (srtio2>sr_skip2*ndz2*ndg2){
                /* SR1 update and increment phs and dgs pointers by dxux */
                ddgs = dgs; 
                for (i=0;i<dxux;++i){
                    aux = *dgs;
                    dddgs = ddgs;
                    for (j=0;j<dxux;++j){
                        *phs -= *dddgs*aux/srtio;
                        ++phs;
                        ++dddgs;
                    }
                    ++dgs; 
                }
                dzs += dxu;
            }
#if DPC
            /* On path AL hessian */
            /* dgp <- dgp + php*dzp and srtio <- (dgp + php*dzp)'*dzp */           
            ddgp = dgp;        
            pphp = php;   
            dddzp = dzp;
            srtio = 0.;
            ndz2 = 0.;
            ndg2 = 0.;
            for (i=0;i<dpc;++i){
                aux = 0.;
                ddzp = dzp;
                for (j=0;j<dpc;++j){
                    aux += *pphp**ddzp;
                    ++pphp;
                    ++ddzp;
                }
                *ddgp += aux;
                auxz = *dddzp;
                auxg = *ddgp;
                srtio += auxg*auxz;
                ndz2 += auxz*auxz;
                ndg2 += auxg*auxg; 
                ++ddgp;
                ++dddzp; 
            } 
            srtio2 = srtio*srtio;
            /* Check denominator and if zero, skip update on current slack-block */
            /* !!!! Change SR1 update criterion !!!! */
        //    if (ABS(srtio)>=POSZER){
            if (srtio2>sr_skip2*ndz2*ndg2){
                /* SR1 update and increment php and dgp pointers by dpc */ 
                ddgp = dgp;
                for (i=0;i<dpc;++i){
                    aux = *dgp;
                    dddgp = ddgp;
                    for (j=0;j<dpc;++j){
                        *php -= *dddgp*aux/srtio;
                        ++php;
                        ++dddgp;
                    }
                    ++dgp;
                }
                dzp += dpc;
            }
#endif
        }
        /**
        ** Trust-region management 
        **/
        /* Adjust trust region radius on first iteration */
        if (it==0)
            trRad = MIN(trRad,norIstp);
        /* Coefficient for trust region update */
        if (vcan-val+gts<=0.){
            alpha = sig3;
        }
        else {
            alpha = MAX(sig1,0.5*gts/(vcan-val+gts));
        }
        /* Actual objective reduction */
        actRed = val-vcan;
        /* Predicted objective reduction */
        preRed = val-vmod;
        /* Update trust-region radius */
    //    testDecrease = 1*(actRed<=eta1*preRed)+2*(actRed>eta1*preRed)*(actRed<eta2*preRed)+3*(actRed>=eta2*preRed);
        testDecrease = 1*(actRed<eta0*preRed)+2*(actRed>=eta0*preRed)*(actRed<eta1*preRed);
        testDecrease += 3*(actRed>=eta1*preRed)*(actRed<eta2*preRed)+4*(actRed>=eta2*preRed);  
        switch (testDecrease){
            case 1:
                /* Decrease trust-region radius */
            //    trRad = 0.5*(sig2*trRad+sig1*MIN(norIstp,trRad));
            //    std::cout<<"nor2stp="<<nor2stp<<std::endl;
            //    std::cout<<"alpha="<<alpha<<std::endl;
                trRad = MIN(MAX(alpha,sig1)*norIstp,sig2*trRad);
                break;
            case 2:
                /* Update trust region radius */
            //    trRad *= 0.5*(sig3+sig1);
            //    std::cout<<"nor2stp="<<nor2stp<<std::endl;
            //    std::cout<<"alpha="<<alpha<<std::endl;
                trRad = MAX(sig1*trRad,MIN(alpha*norIstp,sig2*trRad));
                break;
            case 3:
                /* Update trust region radius */
            //    trRad *= 0.5*(1.0+sig3);
            //    std::cout<<"nor2stp="<<nor2stp<<std::endl;
            //    std::cout<<"alpha="<<alpha<<std::endl;
                trRad = MAX(sig1*trRad,MIN(alpha*norIstp,sig3*trRad));
                break;
            case 4:
                /* Increase trust region radius */
            //    std::cout<<"nor2stp="<<nor2stp<<std::endl;
            //    std::cout<<"alpha="<<alpha<<std::endl;
                trRad = MAX(trRad,MIN(alpha*norIstp,sig3*trRad));
                break;
            default:  
#if PRINT_TR    
                std::cerr<<"PROBLEM IN TRUST RADIUS UPDATE, ABORT WITH KKT="<<kkt<<std::endl;
                std::cout<<"actRed= "<<actRed<<std::endl;
                std::cout<<"preRed= "<<preRed<<std::endl;
                std::cout<<"testDecrease="<<testDecrease<<std::endl;
#endif
                break;
        }     
        /** 
        ** Iteration successful: update iterate to candidate point 
        **/
        testIterate = 0*(actRed<=eta0*preRed)+1*(actRed>eta0*preRed);
        if (testIterate){
#if PRINT_TR
            std::cout<<">>> SUCCESSFUL ITERATION <<<"<<std::endl;
            std::cout<<"actRed= "<<actRed<<std::endl;
            std::cout<<"preRed= "<<preRed<<std::endl;
#endif
            /* Update iterate to candidate point */
            memcpy(z,zS,dZBL);
            /* Update separated shooting Al gradients */
            memcpy(gALshoot,gALshootS,dGshootBL);
#if DPC
            /* Update slack-path of full AL gradient */
            memcpy(gALpath,gALpathS,dGpathBL);
#endif
            /* Build shooting part of full AL gradient at updated iterate */
            memset(gAL,0,dZshootBL);
            g = gAL;
            pgs = gALshoot;
            for (n=0;n<Ns_;++n){
                for (i=0;i<dxux;++i)
                    *(g+i) += *pgs++;
                g += dxu;
            }
            /* Update kkt */
            evalKKTprim();
            /*  Update CG tolerance:
              - sqrt(kkt2) for local superlinear convergence 
              - kkt2 for quadratic local convergence */
#if QUAD_CV
            cgTol2 = kkt2*MIN(0.01,kkt2);    
#endif
#if SUPL_CV
            cgTol2 = kkt2*MIN(0.01,kkt);
#endif
            /* Update vAL value */
            val = vcan;
        }
        /* Go to next iteration */
		++it;
	}
    numPrimIt = it;
}

/**
** Real-time primal loop, truncated trust-region SR1
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::primalRealTime(const double kktTol){

	int it,i,j,n,testIterate,testDecrease;
    double vmod,vcan,val,actRed,preRed;
    double aux,ssz,alpha;
    double auxz,auxg;
    double srtio,srtio2,ndz,ndz2,ndg,ndg2;
    double n2r,n2g,gts;
    double *zz;
	double *zm,*zM;
	double *pz,*dz,*pzs;
	double *dzs,*dzp,*ddzs,*ddzp,*dddzs,*dddzp;
	double *dgp,*pgp,*pgpS,*ddgp,*dddgp;
	double *dgs,*pgs,*pgsS,*ddgs,*dddgs;
    double *g;
	double *phs,*pphs,*php,*pphp;

    /**
    ** Initialization phase  
    **/
	/* Compute augmented Lagrangian and its full gradient */
    val = shooter.shootFromIterate(halfRho);

	/* Get KKT satisfaction */
    evalKKTprim();

    /* Initialize squared tolerance for refinement phase kkt2*min(0.01,...)
    	- kkt for local superlinear convergence 
    	- kkt2 for quadratic local convergence */
#if QUAD_CV
	cgTol2 = kkt2*MIN(0.01,kkt2);
#endif
#if SUPL_CV
    cgTol2 = kkt2*MIN(0.01,kkt);
#endif

	/* Initialize trust-region radius: coef*kkt */
	trRad = 0.1*kkt;

    /**
    ** Trust-region iterations 
    **/
	it = 0;
#if PRINT_TR
    std::cout<<"val initial= "<<std::setprecision(10)<<val<<std::endl;
    std::cout<<"kkt initial= "<<std::setprecision(10)<<kkt<<std::endl;
#endif
#if COUNT_CG
    cumCG = 0;
    cumRestart = 0;
#endif 
	/* Set initial step-size */
    ssz = sini;
    while ((it<maxPit)&&(trRad>trRadTol)&&(kkt>kktTol)&&(cumCG<maxCumCG)){
#if PRINT_TR
        std::cout<<"T-R ITERATION "<<it+1<<std::endl;
        std::cout<<"T-R radius= "<<trRad<<std::endl;
        std::cout<<std::setprecision(5)<<"val= "<<val<<std::endl;
        std::cout<<std::setprecision(5)<<"kkt= "<<kkt<<", kktTol= "<<kktTol<<std::endl;
//        std::cout<<"Running step-size= "<<ssz<<std::endl;
#endif 
		/**
        ** Generate candidate for next iterate inside trust-region 
        **/
        /** 1. Cauchy phase */
#if PRINT_TR
        std::cout<<"Start Cauchy phase..."<<std::endl;
#endif
#if CCH_PRJ_SRCH
        /* Use backtracking projected search */
//        std::cout<<"!!!!! Initial step-size= !!!!"<<ssz<<std::endl;
        activDetector.findProjectedSearch(&ssz,&vmod,&norIstp,trRad);
#if DETEC_CAUCH_ERR
        if (activDetector.getError()){
            errPrim = 1;
            std::cout<<"ERROR IN CAUCHY, interruption at TR iter "<<it<<" with kkt="<<kkt<<std::endl;
            return;
        }
#endif
//        ssz = sini;
#endif
#if CCH_INTER_EXTRA
        /* Use interpolation-extrapolation strategy (see TRON) */
        activDetector.findInterpExtrap(&ssz,&vmod,&norIstp,trRad);
#if DETEC_CAUCH_ERR
        if (activDetector.getError()){
            errPrim = 1;
            std::cout<<"ERROR IN CAUCHY, interruption at TR iter "<<it<<" with kkt="<<kkt<<std::endl;
            return;
        }
#endif
//        ssz = sini;
#endif

        /** 2. Refinement phase */
#if PRINT_TR            
        std::cout<<"Start refinement phase..."<<std::endl;
#endif
#if REFIN_PCG_RES
        /* PCG on active face at Cauchy point with restart if problem bound hit during search */
    //    precSubRefiner.findCandidateRestart(cgTol2,trRad);
        precSubRefiner.findCandidateRestartRT(cgTol2,trRad,cumCG);
#if COUNT_CG
//        std::cout<<"Add "<<precSubRefiner.getCumCG()<<" PCG iters."<<std::endl;
        cumCG += precSubRefiner.getCumCG();
        cumRestart += precSubRefiner.getCumRestart();
#endif
#endif
#if REFIN_MOR_TOR
        /* Projected searches along directions generated by PCG */
        precProjRefiner.findCandidate(vmod,cgTol2,trRad);
#if COUNT_CG
//        std::cout<<"Add "<<precProjRefiner.getCumCG()<<" PCG iters."<<std::endl;
        cumCG += precProjRefiner.getCumCG();
#endif
#endif
#if PRINT_TR
        std::cout<<"TRUST-REGION: kkt satisf.= "<<kkt<<std::endl;
#if COUNT_CG
        std::cout<<"TRUST-REGION: cumul. CG iters= "<<cumCG<<std::endl;
        std::cout<<"TRUST-REGION: cumul. CG restarts= "<<cumRestart<<std::endl;
#endif
#endif
        /* Update negative difference z-zS, compute its infinity norm (and euclidean norm) */   
        dz = negstp;
        zz = z;
        pzs = zS;
        g = gAL;
        norIstp = -1.;
        gts = 0.;
        while (zz != z+dZ){
            aux = *zz++-*pzs++;
            *dz++ = aux;
            gts += *g++*aux;
            norIstp = MAX(norIstp,ABS(aux));
        }

        /* Evaluate model at candidate */
#if REFIN_PCG_RES
        vmod = val+qmod.evalModelCandidate();        
#endif
#if REFIN_MOR_TOR
        vmod = val+precProjRefiner.getValMod();
#endif       
        /* Evaluate AL and separated gradients at candidate */
        vcan = shooter.shootFromCandidate(halfRho);

        /* Compute shooting AL gradient difference gALshoot(zS)-gALshoot(z) */
        dgs = dffGshoot;
        pgs = gALshoot; 
        pgsS = gALshootS;
        while (pgs != gALshoot+dGshoot)
            *dgs++ = *pgsS++-*pgs++;
#if DPC
        /* Compute path-slacks AL gradient difference gALpath(zS)-gALpath(z) */
        dgp = dffGpath;
        pgp = gALpath;
        pgpS = gALpathS;
        while (pgp != gALpath+dGpath)
        	*dgp++ = *pgpS++-*pgp++;
#endif
        /** 
        ** SR1 updates on separated blocks 
        **/
        phs = hALshoot;    
        dgs = dffGshoot;
        dzs = negstp;
#if DPC
        php = hALpath;
        dgp = dffGpath;
        dzp = dzs+dZshoot;
#endif
        /* Loop over shooting blocks */
        for (n=0;n<Ns_;++n){ 
            /* On shooting AL hessian */
            /* dgs<-dgs+phs*dzs and srtio<-(dgs+phs*dzs)'*dzs */
            ddgs = dgs;
            pphs = phs;
            dddzs = dzs;
            srtio = 0.;
            ndz2 = 0.;
            ndg2 = 0.;
            for (i=0;i<dxux;++i){
                aux = 0.;
                ddzs = dzs;
                for (j=0;j<dxux;++j){
                    aux += *pphs**ddzs;
                    ++pphs;
                    ++ddzs;
                }
                *ddgs += aux;
                auxz = *dddzs;
                auxg = *ddgs;
                srtio += auxg*auxz;
                ndz2 += auxz*auxz;
                ndg2 += auxg*auxg;
                ++ddgs;
                ++dddzs;
            } 
            srtio2 = srtio*srtio;
            /* Check denominator and if zero skip update on current shooting-block */
            /* !!!! Change SR1 update criterion !!!! */
            //if (ABS(srtio)>=POSZER){
            if (srtio2>sr_skip2*ndz2*ndg2){
                /* SR1 update and increment phs and dgs pointers by dxux */
                ddgs = dgs; 
                for (i=0;i<dxux;++i){
                    aux = *dgs;
                    dddgs = ddgs;
                    for (j=0;j<dxux;++j){
                        *phs -= *dddgs*aux/srtio;
                        ++phs;
                        ++dddgs;
                    }
                    ++dgs; 
                }
                dzs += dxu;
            }
#if DPC
            /* On path AL hessian */
            /* dgp <- dgp + php*dzp and srtio <- (dgp + php*dzp)'*dzp */           
            ddgp = dgp;        
            pphp = php;   
            dddzp = dzp;
            srtio = 0.;
            ndz2 = 0.;
            ndg2 = 0.;
            for (i=0;i<dpc;++i){
                aux = 0.;
                ddzp = dzp;
                for (j=0;j<dpc;++j){
                    aux += *pphp**ddzp;
                    ++pphp;
                    ++ddzp;
                }
                *ddgp += aux;
                auxz = *dddzp;
                auxg = *ddgp;
                srtio += auxg*auxz;
                ndz2 += auxz*auxz;
                ndg2 += auxg*auxg; 
                ++ddgp;
                ++dddzp; 
            } 
            srtio2 = srtio*srtio;
            /* Check denominator and if zero, skip update on current slack-block */
            /* !!!! Change SR1 update criterion !!!! */
        //    if (ABS(srtio)>=POSZER){
            if (srtio2>sr_skip2*ndz2*ndg2){
                /* SR1 update and increment php and dgp pointers by dpc */ 
                ddgp = dgp;
                for (i=0;i<dpc;++i){
                    aux = *dgp;
                    dddgp = ddgp;
                    for (j=0;j<dpc;++j){
                        *php -= *dddgp*aux/srtio;
                        ++php;
                        ++dddgp;
                    }
                    ++dgp;
                }
                dzp += dpc;
            }
#endif
        }
        /**
        ** Trust-region management 
        **/
        /* Adjust trust region radius on first iteration */
        if (it==0)
            trRad = MIN(trRad,norIstp);
        /* Coefficient for trust region update */
        if (vcan-val+gts<=0.){
            alpha = sig3;
        }
        else {
            alpha = MAX(sig1,0.5*gts/(vcan-val+gts));
        }
        /* Actual objective reduction */
        actRed = val-vcan;
        /* Predicted objective reduction */
        preRed = val-vmod;
        /* Update trust-region radius */
    //    testDecrease = 1*(actRed<=eta1*preRed)+2*(actRed>eta1*preRed)*(actRed<eta2*preRed)+3*(actRed>=eta2*preRed);
        testDecrease = 1*(actRed<eta0*preRed)+2*(actRed>=eta0*preRed)*(actRed<eta1*preRed);
        testDecrease += 3*(actRed>=eta1*preRed)*(actRed<eta2*preRed)+4*(actRed>=eta2*preRed);  
        switch (testDecrease){
            case 1:
                /* Decrease trust-region radius */
            //    trRad = 0.5*(sig2*trRad+sig1*MIN(norIstp,trRad));
            //    std::cout<<"nor2stp="<<nor2stp<<std::endl;
            //    std::cout<<"alpha="<<alpha<<std::endl;
                trRad = MIN(MAX(alpha,sig1)*norIstp,sig2*trRad);
                break;
            case 2:
                /* Update trust region radius */
            //    trRad *= 0.5*(sig3+sig1);
            //    std::cout<<"nor2stp="<<nor2stp<<std::endl;
            //    std::cout<<"alpha="<<alpha<<std::endl;
                trRad = MAX(sig1*trRad,MIN(alpha*norIstp,sig2*trRad));
                break;
            case 3:
                /* Update trust region radius */
            //    trRad *= 0.5*(1.0+sig3);
            //    std::cout<<"nor2stp="<<nor2stp<<std::endl;
            //    std::cout<<"alpha="<<alpha<<std::endl;
                trRad = MAX(sig1*trRad,MIN(alpha*norIstp,sig3*trRad));
                break;
            case 4:
                /* Increase trust region radius */
            //    std::cout<<"nor2stp="<<nor2stp<<std::endl;
            //    std::cout<<"alpha="<<alpha<<std::endl;
                trRad = MAX(trRad,MIN(alpha*norIstp,sig3*trRad));
                break;
            default:  
#if PRINT_TR    
                std::cerr<<"PROBLEM IN TRUST RADIUS UPDATE, ABORT WITH KKT="<<kkt<<std::endl;
                std::cout<<"actRed= "<<actRed<<std::endl;
                std::cout<<"preRed= "<<preRed<<std::endl;
                std::cout<<"testDecrease="<<testDecrease<<std::endl;
#endif
                break;
        }     
        /** 
        ** Iteration successful: update iterate to candidate point 
        **/
        testIterate = 0*(actRed<=eta0*preRed)+1*(actRed>eta0*preRed);
        if (testIterate){
#if PRINT_TR
            std::cout<<">>> SUCCESSFUL ITERATION <<<"<<std::endl;
            std::cout<<"actRed= "<<actRed<<std::endl;
            std::cout<<"preRed= "<<preRed<<std::endl;
#endif
            /* Update iterate to candidate point */
            memcpy(z,zS,dZBL);
            /* Update separated shooting Al gradients */
            memcpy(gALshoot,gALshootS,dGshootBL);
#if DPC
            /* Update slack-path of full AL gradient */
            memcpy(gALpath,gALpathS,dGpathBL);
#endif
            /* Build shooting part of full AL gradient at updated iterate */
            memset(gAL,0,dZshootBL);
            g = gAL;
            pgs = gALshoot;
            for (n=0;n<Ns_;++n){
                for (i=0;i<dxux;++i)
                    *(g+i) += *pgs++;
                g += dxu;
            }
            /* Update kkt */
            evalKKTprim();
            /*  Update CG tolerance:
              - sqrt(kkt2) for local superlinear convergence 
              - kkt2 for quadratic local convergence */
#if QUAD_CV
            cgTol2 = kkt2*MIN(0.01,kkt2);    
#endif
#if SUPL_CV
            cgTol2 = kkt2*MIN(0.01,kkt);
#endif
            /* Update vAL value */
            val = vcan;
        }
        /* Go to next iteration */
		++it;
	}
#if PRINT_DUA
    numPrimIt = it;
#endif
#if PRINT_TR
    std::cout<<"Aborted TR iterations, max. primal iters: "<<maxPit<<std::endl;
    std::cout<<"Break after "<<it<<" iterations with kkt: "<<kkt<<std::endl;
    std::cout<<"val final: "<<std::setprecision(10)<<val<<std::endl;
    std::cout<<"trRad="<<trRad<<"<= trRadTol="<<trRadTol<<std::endl;
#endif
}

/**************************************** Dual loops *******************************************/
/**
** Real-time dual loop (only one dual update)
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::solveRealTime(const double kktTol){

    double *mm,*ee;

#if JACO_PRE
#if REFIN_PCG_RES
    precSubRefiner.setRegCoefPrec(1.);
#endif
#endif
    /* Fixed number of primal iterations */
    primalRealTime(kktTol);
    /* Evaluate squared 2-norm of equality constraints */
    evalEconNorm();

#if PRINT_DUA
    std::cout<<"kkt= "<<kkt<<std::endl;
    std::cout<<"nec= "<<nec<<std::endl;
    std::cout<<numPrimIt<<" primal iterations."<<std::endl;
#if COUNT_CG
    std::cout<<cumCG<<" CG iters."<<std::endl;
    std::cout<<cumRestart<<" CG restarts."<<std::endl;
#endif
#endif
    /* One dual update */
    mm = mu;
    ee = econ;
    while (mm != mu+dL)
        *mm++ += 2.*halfRho*(*ee++);
    /* Extract optimal control input */
    memcpy(uopt,z+dx,dubl);
}

/**
** LANCELOT outer loop: adapt tolerance on primal optimality depending on constraints satisfaction
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::solveLancelot(){

    int dit;
    double *mm,*ee;
    double ecTol,kktTol;
    double rho;

    halfRho = halfRhoIni;
    /**
    * Outer loop
    */
    ecTol = ecTolIni;
    kktTol = kktTolIni;
    globCG = 0;
    globRestart = 0;
    dit = 0;
    while (1){
        /* Stop if max. # iterations reached */
        if (dit>=maxDit)
            break;
#if JACO_PRE
#if REFIN_PCG_RES
        precSubRefiner.setRegCoefPrec(1.);
#endif
#endif
        /* Set initial SR1 hessian approximation */
        initHessianApprox();
        /* Primal iterations stopped when kkt satisfaction below kktTol */
        primalFull(kktTol);
        /* Evaluate squared 2-norm of equality constraints */
        evalEconNorm();

        globCG += cumCG;
        globRestart += cumRestart;

#if PRINT_DUA
        std::cout<<"Outer iter "<<dit+1<<": kkt="<<kkt<<std::endl;
        std::cout<<"Outer iter "<<dit+1<<": nec="<<nec<<std::endl;
        std::cout<<"Outer iter "<<dit+1<<": "<<numPrimIt<<" primal iterations."<<std::endl;
#if COUNT_CG
        std::cout<<"Outer iter "<<dit+1<<": "<<cumCG<<" CG iters."<<std::endl;
        std::cout<<"Outer iter "<<dit+1<<": "<<cumRestart<<" CG restarts."<<std::endl;
#endif
        std::cout<<"Outer iter "<<dit+1<<", ec tolerance: "<<ecTol<<std::endl;
        std::cout<<"Outer iter "<<dit+1<<", kkt tolerance: "<<kktTol<<std::endl;
#endif

        /* If nan, abort */
        if ((isnan(kkt))||(kkt>=1E5)||(isnan(nec))||(nec>=1E5)){
            std::cout<<"norIstp= "<<norIstp<<std::endl;
            std::cerr<<"Abort, kkt="<<kkt<<", nec= "<<nec<<std::endl;
            displayIterate();
            return;
        }

        /* Abort if error in primal loop */
        if (errPrim) 
            return;

        /* Stop if absolute feasibility tolerance reached */
        if (nec<=necTolAbs)
            break;

        /* Lancelot-like update */
        if (nec<=ecTol){
            rho = 2.*halfRho;
#if PRINT_DUA
            std::cout<<"Sufficiently feasible, update dual variable"<<std::endl;
#endif
            /* First-order dual update, freeze penalty */
            mm = mu;
            ee = econ;
            while (mm != mu+dL)
                *mm++ += rho*(*ee++);
            
            /* Update tolerances */
            ecTol /= pow(rho,0.9);
            kktTol /= rho;
        }
        else {
#if PRINT_DUA
            std::cout<<"Not sufficiently feasible, increase penalty parameter"<<std::endl;
#endif
            /* Increase penalty, freeze dual */
            halfRho *= mulPen;
            rho = 2.*halfRho;
            /* Update tolerances */
            ecTol = 0.1/pow(rho,0.1);
            kktTol = 1./rho; 
        }
        /* Go to next outer iteration */
        ++dit;
    }
    /* Extract optimal control input */
    memcpy(uopt,z+dx,dubl);
}

/**
** Full dual loop (classic augmented Lagrangian outer iterations)
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::solveFull(){

    int dit;
    double *mm,*ee;

//    std::cout<<"sta solveFull"<<std::endl;

    /* Set initial penalty parameter */
    halfRho = halfRhoIni;

    /* Set initial guess on primal optimizer */
/*    xes = shooter.getXest();
    zz = z;
    for (n=0;n<Ns_;++n){
        memcpy(zz,xes,dxbl);
        zz += dxu;
    } */

    /**
    ** Outer (dual) loop
    **/
    globCG = 0;
    globRestart = 0;
    dit = 0;
    while (1){
        /* Stop if max. # dual steps reached */
        if (dit>=maxDit)
            break;

#if JACO_PRE
#if REFIN_PCG_RES
        precSubRefiner.setRegCoefPrec(1.);
#endif
#endif
        /* Set initial SR1 approx. */
        initHessianApprox();
        /* Primal iterations */
        primalFull(kktTolAbs);
        /* Evaluate squared 2-norm of equality constraints */
        evalEconNorm();

        globCG += cumCG;
        globRestart += cumRestart;
#if PRINT_DUA
        std::cout<<"Dual iter "<<dit+1<<": kkt="<<kkt<<std::endl;
        std::cout<<"Dual iter "<<dit+1<<": nec="<<nec<<std::endl;
        std::cout<<"Dual iter "<<dit+1<<": "<<numPrimIt<<" primal iterations."<<std::endl;
#if COUNT_CG
        std::cout<<"Dual iter "<<dit+1<<": "<<cumCG<<" CG iters."<<std::endl;
        std::cout<<"Dual iter "<<dit+1<<": "<<cumRestart<<" CG restarts."<<std::endl;
#endif
#endif
        /* Stop if feasibility tolerance reached */
        if (nec<=necTolAbs)
            break;
        /* If nan, abort */
        if ((isnan(kkt))||(kkt>=1E5)||(isnan(nec))||(nec>=1E5)){
            std::cout<<"norIstp= "<<norIstp<<std::endl;
            std::cerr<<"Abort, kkt="<<kkt<<", nec= "<<nec<<std::endl;
            displayIterate();
            return;
        }
        /* Dual update */
        mm = mu;
        ee = econ;
        while (mm != mu+dL)
            *mm++ += 2.*halfRho*(*ee++);
        /* Increase penalty parameter */
        halfRho *= mulPen;

        ++dit;
    }
    /* Extract optimal control input */
    memcpy(uopt,z+dx,dubl);
}

/******************************************************
******************** For debugging ******************
********************************************************/
/**
** Display members: primal optimizer z, dual optimizer mu,...
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::displayMembers(const int acc){

    int i;

    std::cout<<"Trust-region solver: z="<<std::endl;
    for (i=0;i<dZ;++i)
        std::cout<<std::setprecision(acc)<<z[i]<<std::endl;

/*    std::cout<<"Trust-region solver: zS="<<std::endl;
    for (i=0;i<dZ;++i)
        std::cout<<zS[i]<<std::endl; */

    std::cout<<"Trust-region solver: econ="<<std::endl;
    for (i=0;i<dL;++i)
        std::cout<<std::setprecision(acc)<<econ[i]<<std::endl;

    std::cout<<"Trust-region solver: mu="<<std::endl;
    for (i=0;i<dL;++i)
        std::cout<<mu[i]<<std::endl; 

/*    std::cout<<"Trust-region solver: gAL="<<std::endl;
    for (i=0;i<dZ;++i)
        std::cout<<gAL[i]<<std::endl; */
}

/**
** Display problem bounds
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::displayBounds(){

    int i;

    for (i=0;i<dZ;++i)
        std::cout<<"zmin["<<i<<"]="<<zmin[i]<<", zmax["<<i<<"]="<<zmax[i]<<std::endl;

}

/**
** Displays candidate point zS 
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::displayCandidate(){

    int i;

    std::cout<<"Candidate zS= "<<std::endl;
    for (i=0;i<dZ;++i)
        std::cout<<zS[i]<<std::endl;
}

/** 
** Displays primal iterate z
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::displayIterate(){

    int i;

    std::cout<<"Iterate z= "<<std::endl;
    for (i=0;i<dZ;++i)
        std::cout<<z[i]<<std::endl;
}

/**
** Displays dual 
**/
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::displayDual(){
    
    double *mm;

    std::cout<<"Dual variable mu="<<std::endl;
    mm = mu;
    while (mm != mu+dL)
        std::cout<<*mm++<<std::endl;
}

/* Display shooting part of AL hessian */
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::displayHalShoot(){

    int n,i;
    double *hal;

    std::cout<<"hALshoot: "<<std::endl;    
    hal = hALshoot;
    for (n=0;n<Ns_;++n){
        std::cout<<"Block "<<n<<std::endl;
        for (i=0;i<dhS;++i){  
            if (i%dxux == dxux-1){
                std::cout<<*hal++<<std::endl;
            }
            else {
                std::cout<<*hal++<<" ";
            }
        }
    }

}

/* Read data from MAT files and store in z & mu */
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::readMATfiles(){

    // TO BE IMPLEMENTED

}

/* Display initial state xest */
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::displayXest(){

    double *xe,*xest;

    std::cout<<"OCP solver, initial state xest="<<std::endl;
    xest = shooter.getXest();
    xe = xest;
    while (xe != xest+dx)
        std::cout<<*xe++<<std::endl;

}
#if SCAL_VAR
/* Display state transformation */
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::displayStateTransfo(){
    shooter.displayStateTransfo();
}

/* Display input transfo */
template <class dynT,class costT,class pconT,class mayerT> void NonlinOCPsolver<dynT,costT,pconT,mayerT>::displayInputTransfo(){
    shooter.displayInputTransfo();
}
#endif