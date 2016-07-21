#include <string.h>
#include <iomanip>

#include "SpgOCPsolver.h"
#include "macros.h"


/**
* Default constructor 
*/
template <class dynT,class costT,class pconT,class mayerT> SpgOCPsolver<dynT,costT,pconT,mayerT>::SpgOCPsolver(){

	z = NULL;
    znex = NULL;
	dffz = NULL;
	uopt = NULL;

	zmin = NULL;
	zmax = NULL;

	mu = NULL;
	econ = NULL;

	gAL = NULL;
    dffGal = NULL;
	gALshoot = NULL;
	gALpath = NULL;
	dffGpath = NULL;

    vALhist = NULL;
} 

/**
* Default destructor
*/
template <class dynT,class costT,class pconT,class mayerT> SpgOCPsolver<dynT,costT,pconT,mayerT>::~SpgOCPsolver(){
	

	if (z != NULL)
		delete [] z;

    if (znex != NULL)
        delete [] znex;

	if (dffz != NULL)
		delete [] dffz;

	if (uopt != NULL)
		delete [] uopt;


	if (zmin != NULL)
		delete [] zmin;

	if (zmax != NULL)
		delete [] zmax;


	if (mu != NULL)
		delete [] mu;

	if (econ != NULL)
		delete [] econ;


	if (gAL != NULL)
		delete [] gAL;

    if (dffGal != NULL)
        delete [] dffGal;

	if (gALshoot != NULL)
		delete [] gALshoot;

	if (dffGpath != NULL)
		delete [] dffGpath;

    if (vALhist != NULL)
        delete [] vALhist;
}

/**
* Copy constructor
*/
template <class dynT,class costT,class pconT,class mayerT> SpgOCPsolver<dynT,costT,pconT,mayerT>::SpgOCPsolver(const SpgOCPsolver<dynT,costT,pconT,mayerT>& obj){
	

	shooter = obj.shooter;


    Ns_ = obj.Ns_;

	dZshoot = obj.dZshoot;
   	dZpath = obj.dZpath;
    dZ = obj.dZ;
    dZBL = obj.dZBL;
   
    dLshoot = obj.dLshoot;
    dLpath = obj.dLpath;
    dL = obj.dL;
    dLBL = obj.dLBL;
    
    du = obj.du;
    dubl = obj.dubl;
   	dx = obj.dx;
	dxu = obj.dxu;
    dxux = obj.dxux;
    dpc = obj.dpc;
    
    dGshoot = obj.dGshoot;
    dGshootBL = obj.dGshootBL;
   	dGpath = obj.dGpath;
    dGpathBL = obj.dGpathBL;
    
    dhS = obj.dhS;
    dhP = obj.dhP;
   	dHshoot = obj.dHshoot; 
    dHpath = obj.dHpath;

    /* Primal variables */
    z = new double[dZ];
    std::copy(obj.z,obj.z+dZ,z);

    znex = new double[dZ];
    std::copy(obj.znex,obj.znex+dZ,znex);

    dffz = new double[dZ];
    std::copy(obj.dffz,obj.dffz+dZ,dffz);

    uopt = new double[du];
    std::copy(obj.uopt,obj.uopt+du,uopt);


    /* Bounds on primal optimizer */
    zmin = new double[dZ];
    std::copy(obj.zmin,obj.zmin+dZ,zmin);

    zmax = new double[dZ];
    std::copy(obj.zmax,obj.zmax+dZ,zmax);


    /* Duals and equality constraints */
    mu = new double[dL];
    std::copy(obj.mu,obj.mu+dL,mu);

    econ = new double[dL];
    std::copy(obj.econ,obj.econ+dL,econ);


    /* Gradients */
    gAL = new double[dZ];
    std::copy(obj.gAL,obj.gAL+dZ,gAL);

    dffGal = new double[dZ];
    std::copy(obj.dffGal,obj.dffGal+dZ,dffGal);

    gALshoot = new double[dGshoot];
    std::copy(obj.gALshoot,obj.gALshoot+dGshoot,gALshoot);
    
    gALpath = obj.gALpath;
    
  	dffGpath = new double[dGpath];
  	std::copy(obj.dffGpath,obj.dffGpath+dGpath,dffGpath);

  	
  	/* Primal solver options */
	maxPit = obj.maxPit;

    gamSpg = obj.gamSpg;

    histLenSpg = obj.histLenSpg;

    sigOneSpg = obj.sigOneSpg;
    sigTwoSpg = obj.sigTwoSpg;

    aminSpg = obj.aminSpg;
    amaxSpg = obj.amaxSpg;

    ainiSpg = obj.ainiSpg;


    /* Control parameters */
    halfRho = obj.halfRho;
    halfRhoIni = obj.halfRhoIni;

    kkt2 = obj.kkt2;
    kkt = obj.kkt;
    nec2 = obj.nec2;
    nec = obj.nec;
    kktTol2 = obj.kktTol2;
    necTol2 = obj.necTol2;
    diffTol = obj.diffTol;


    /* History for non-monotone line-search */
    vALhist = new double[histLenSpg];
    std::copy(obj.vALhist,obj.vALhist+histLenSpg,vALhist);
}


/**
* Assignment operator
*/
template <class dynT,class costT,class pconT,class mayerT> SpgOCPsolver<dynT,costT,pconT,mayerT>& SpgOCPsolver<dynT,costT,pconT,mayerT>::operator= (const SpgOCPsolver<dynT,costT,pconT,mayerT>& obj){
	
	double *z_,*znex_,*dffz_,*uopt_;
    double *zmin_,*zmax_;
	double *mu_,*econ_;
	double *gAL_,*dffGal_;
	double *gALshoot_,*dffGshoot_;
	double *dffGpath_;
    double *vALhist_;


	if (this != &obj){

		shooter = obj.shooter;

        /* Dimensions */
        Ns_ = obj.Ns_;

		dZshoot = obj.dZshoot;
   		dZpath = obj.dZpath;
    	dZ = obj.dZ;
    	dZBL = obj.dZBL;
   
    	dLshoot = obj.dLshoot;
    	dLpath = obj.dLpath;
    	dL = obj.dL;
    	dLBL = obj.dLBL;
    
    	du = obj.du;
    	dubl = obj.dubl;
   		dx = obj.dx;
		dxu = obj.dxu;
    	dxux = obj.dxux;
    	dpc = obj.dpc;
    
    	dGshoot = obj.dGshoot;
    	dGshootBL = obj.dGshootBL;
   		dGpath = obj.dGpath;
    	dGpathBL = obj.dGpathBL;
    
    	dhS = obj.dhS;
    	dhP = obj.dhP;
   		dHshoot = obj.dHshoot; 
    	dHpath = obj.dHpath;

    	/* Options */
		maxPit = obj.maxPit;

    	gamSpg = obj.gamSpg;

    	histLenSpg = obj.histLenSpg;

    	sigOneSpg = obj.sigOneSpg;
    	sigTwoSpg = obj.sigTwoSpg;

    	aminSpg = obj.aminSpg;
    	amaxSpg = obj.amaxSpg;
    	ainiSpg = obj.ainiSpg;

        /* Control parameters */
        halfRhoIni = obj.halfRhoIni;
        halfRho = obj.halfRho;

        kkt2 = obj.kkt2;
        kkt = obj.kkt;

        nec2 = obj.nec2;
        nec = obj.nec;

        kktTol2 = obj.kktTol2;
        necTol2 = obj.necTol2;

        diffTol = obj.diffTol;

    	/* */
    	if (z != NULL)
    		delete [] z;
    	if (obj.z != NULL){
    		z_ = new double[dZ]; 
    		std::copy(obj.z,obj.z+dZ,z_);
    		z = z_;
    	}
    	else {
    		z = NULL;
    	}

        /* */
        if (znex != NULL)
            delete [] znex;
        if (obj.znex != NULL){
            znex_ = new double[dZ];
            std::copy(obj.znex,obj.znex+dZ,znex_);
            znex = znex_;
        }
        else {
            znex = NULL;
        }

    	/* */
    	if (dffz != NULL)
    		delete [] dffz;
    	if (obj.dffz != NULL){
    		dffz_ = new double[dZ];
    		std::copy(obj.dffz,obj.dffz+dZ,dffz_);
    		dffz = dffz_;
    	}
    	else {
    		dffz = NULL;
    	}

        /* */
        if (uopt != NULL)
            delete [] uopt;
        if (obj.uopt != NULL){
            uopt_ = new double[du];
            std::copy(obj.uopt,obj.uopt+du,uopt_);
            uopt = uopt_;
        }
        else {
            uopt = NULL;
        }

    	/* */
    	if (zmin != NULL)
    		delete [] zmin;
    	if (obj.zmin != NULL){
    		zmin_ = new double[dZ];
    		std::copy(obj.zmin,obj.zmin+dZ,zmin_);
    		zmin = zmin_;
    	}
    	else {
    		zmin = NULL;
    	}

    	/* */
    	if (zmax != NULL)
    		delete [] zmax;
    	if (obj.zmax != NULL){
    		zmax_ = new double[dZ];
    		std::copy(obj.zmax,obj.zmax+dZ,zmax_);
    		zmax = zmax_;
    	}
    	else {
    		zmax = NULL;
    	}

    	/* */
    	if (mu != NULL)
    		delete [] mu;
    	if (obj.mu != NULL){
    		mu_ = new double[dL];
    		std::copy(obj.mu,obj.mu+dL,mu_);
    		mu = mu_;
    	}
    	else {
    		mu = NULL;
    	}

    	/* */
    	if (econ != NULL)
    		delete [] econ;
    	if (obj.econ != NULL){
    		econ_ = new double[dL];
    		std::copy(obj.econ,obj.econ+dL,econ_);
    		econ = econ_;
    	}
    	else {
    		econ = NULL;
    	}

    	/* */
    	if (gAL != NULL)
            delete [] gAL;
        if (obj.gAL != NULL){
            gAL_ = new double[dZ];
            std::copy(obj.gAL,obj.gAL+dZ,gAL_);
            gAL = gAL_;
        }
        else {
            gAL = NULL;
        }

        /* */
        if (dffGal != NULL)
            delete [] dffGal;
        if (obj.dffGal != NULL){
            dffGal_ = new double[dZ];
            std::copy(obj.dffGal,obj.dffGal+dZ,dffGal_);
            dffGal = dffGal_;
        }
        else {
            dffGal = NULL;
        }

        /* */
        if (gALshoot != NULL)
            delete [] gALshoot;
        if (obj.gALshoot != NULL){
            gALshoot_ = new double[dGshoot];
            std::copy(obj.gALshoot,obj.gALshoot+dGshoot,gALshoot_);
            gALshoot = gALshoot_;
        }
        else {
            gALshoot = NULL; 
        }

        /* */
        gALpath = obj.gALpath;

        /* */
        if (dffGpath != NULL)
            delete [] dffGpath;
        if (obj.dffGpath != NULL){
            dffGpath_ = new double[dGpath];
            std::copy(obj.dffGpath,obj.dffGpath+dGpath,dffGpath_);
            dffGpath = dffGpath_;
        }
        else {
            dffGpath = NULL;
        }

        /* */
        if (vALhist != NULL)
            delete [] vALhist;
        if (obj.vALhist != NULL){
            vALhist_ = new double[histLenSpg];
            std::copy(obj.vALhist,obj.vALhist+histLenSpg,vALhist_);
            vALhist = vALhist_;
        }
        else {
            vALhist = NULL;
        }
	}
	return *this;     
} 

/**
* Initialization method
*/
template <class dynT,class costT,class pconT,class mayerT> void SpgOCPsolver<dynT,costT,pconT,mayerT>::init(){

	int i,n,ibgx,indx,ibgu,indu;
    double *xmin,*xmax,*umin,*umax,*zl,*zu;

    /* Initialize shooter */
    shooter.init();
    shooter.displayShootingParameters();

    /* Initialize dimensions */
    Ns_ = shooter.getNs_();
    dZshoot = shooter.getDprimShoot();
    dZpath = shooter.getDprimPath();
    dZ = shooter.getDprim();
    dZBL = dZ*SZDBL;
    dL = shooter.getDdua();
    dLBL = dL*SZDBL;
    du = shooter.getDu();
    dubl = du*SZDBL;
    dx = shooter.getDx();
    dxu = dx+du;
    dxux = dxu+dx;
    dpc = shooter.getDpc();
    dGshoot = Ns_*dxux;
    dGshootBL = dGshoot*SZDBL;
    dGpath = Ns_*dpc;
    dGpathBL = dGpath*SZDBL;
    dhS = dxux*dxux;
    dhP = dpc*dpc;
    dHshoot = Ns_*dhS;
    dHpath = Ns_*dhP;

    /* Allocate primal optimizer */
    if (z != NULL)
        delete [] z;
    z = new double[dZ];
    memset(z,0,dZBL);

    /* Allocate next primal optimizer */
    if (znex != NULL)
        delete [] znex;
    znex = new double[dZ];
    memset(znex,0,dZBL);

    /* Allocate difference between iterates */
    if (dffz != NULL)
        delete [] dffz;
    dffz = new double[dZ];
    memset(dffz,0,dZBL);

    /* Allocate optimal input */
    if (uopt != NULL)
        delete [] uopt;
    uopt = new double[du];

    /* Allocate dual optimizer */
    if (mu != NULL)
        delete [] mu;
    mu = new double[dL];
    
    /* Allocate equality constraints */
    if (econ != NULL)
        delete [] econ;
    econ = new double[dL];
    
    /* Allocate primal bounds (fixed and varying) */
    if (zmin != NULL)
        delete [] zmin;
    zmin = new double[dZ];

    if (zmax != NULL)
        delete [] zmax;
    zmax = new double[dZ];

    
    /* Initialize bounds on slacks for path-constraints */    
    memset(zmin+dZpath,0,dZpath*SZDBL);
    memset(zmax+dZpath,POSINF,dZpath*SZDBL);

    /* Initialize bounds zmin & zmax */
    umin = shooter.getUmin();   
    umax = shooter.getUmax();
    xmin = shooter.getXmin();
    xmax = shooter.getXmax();
        
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
    std::copy(xmax,xmax+dx,zu);    

    /* Initialize AL gradients */
    if (gAL != NULL)
        delete [] gAL;
    gAL = new double[dZ];
    
    if (dffGal != NULL)
        delete [] dffGal;
    dffGal = new double[dZ];

    gALpath = gAL+dZshoot;
    
    if (gALshoot != NULL)
        delete [] gALshoot;
    gALshoot = new double[dGshoot];
    
    /* Initialize differences of AL gradients */
    if (dffGpath != NULL)
        delete [] dffGpath;
    dffGpath = new double[dGpath];

    /* Initialize vAL history for non-monotone line-search */
    if (vALhist != NULL)
        delete [] vALhist;

    if (histLenSpg<0){
        std::cerr<<"ERROR<SpgOCPsolver>: history length nonpositive, abort."<<std::endl;
        return;
    }

    vALhist = new double[histLenSpg];
    
    /* Shooter initialization (members must have been allocated in Simulator) */
    shooter.setPrimal(z);
    shooter.setCandidate(znex);
    shooter.setDual(mu);
    shooter.setEcon(econ);
    shooter.setGradient(gAL);
    shooter.setGradientShooting(gALshoot);
    shooter.setGradientPath(gALpath);
    shooter.setDiffGradientPath(dffGpath);

    /* Set penalty parameter to initial penalty parameter in options */
    halfRho = halfRhoIni;
}

/**
* Display members: primal optimizer z, dual optimizer mu,...
*/
template <class dynT,class costT,class pconT,class mayerT> void SpgOCPsolver<dynT,costT,pconT,mayerT>::displayMembers(){

    int i;

    std::cout<<"SPG solver: z="<<std::endl;
    for (i=0;i<dZ;++i)
        std::cout<<z[i]<<std::endl;

    std::cout<<"SPG solver: mu="<<std::endl;
    for (i=0;i<dL;++i)
        std::cout<<mu[i]<<std::endl;

    std::cout<<"SPG solver: econ="<<std::endl;
    for (i=0;i<dL;++i)
        std::cout<<econ[i]<<std::endl;

    std::cout<<"SPG solver: gAL="<<std::endl;
    for (i=0;i<dZ;++i)
        std::cout<<gAL[i]<<std::endl;
}

/** 
* Initialize primal optimizer 
*/
template <class dynT,class costT,class pconT,class mayerT> void SpgOCPsolver<dynT,costT,pconT,mayerT>::setPrimal(double* z_){
    assert(dZ>=1);
    std::copy(z_,z_+dZ,z);
}

/**
* Initialize dual optimizer 
*/
template <class dynT,class costT,class pconT,class mayerT> void SpgOCPsolver<dynT,costT,pconT,mayerT>::setDual(double* mu_){
    assert(dL>=1);
    std::copy(mu_,mu_+dL,mu);
}

/**
* Cold-start primal-dual optimizer
*/
template <class dynT,class costT,class pconT,class mayerT> void SpgOCPsolver<dynT,costT,pconT,mayerT>::coldStartPrimalDual(){
    
    std::cout<<"Cold-start primal optimizer."<<std::endl;
    /* Set primal to 0 */
    memset(z,0,dZBL);

    std::cout<<"Cold-start dual optimizer."<<std::endl;
    /* Set dual to 0 */
    memset(mu,0,dLBL);
}

/**
* Shift primal & dual variables in an NMPC-like fashion
*/
template <class dynT,class costT,class pconT,class mayerT> void SpgOCPsolver<dynT,costT,pconT,mayerT>::shiftPrimalDual(){

    std::cout<<"Shift primal optimizer."<<std::endl;  
    /* Shift shooting-part of primal sequence */
    std::copy(z+dxu,z+dZshoot,z); 
    /* Shift slacks-part of primal sequence */  
    std::copy(z+dZshoot+dpc,z+dZ,z+dZshoot);
  
    std::cout<<"Shift dual optimizer."<<std::endl;  
    /* Shift shooting-part of dual sequence */
    std::copy(mu+dx,mu+dx+dLshoot,mu);
    /* Shift slacks-part of dual sequence */
    std::copy(mu+dx+dLshoot+dpc,mu+dL,mu+dx+dLshoot);
}

/**
* Evaluate squared primal KKT satisfaction \|P(x-g(x))-x\|_2^2
*/
template <class dynT,class costT,class pconT,class mayerT> void SpgOCPsolver<dynT,costT,pconT,mayerT>::evalKKTprim(){

    int k;
    double dd,az,al,au,ag;

    kkt2 = 0.;
    for (k=0;k<dZ;++k){

        az = *(z+k);
        ag = *(gAL+k);
        al = *(zmin+k);
        au = *(zmax+k);

        dd = PROJ(az-ag,al,au)-az;

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

        kkt2 += dd*dd;
    }

    kkt = sqrt(kkt2);
}

/** 
* Evaluate squared 2-norm of equality constraints
*/
template <class dynT,class costT,class pconT,class mayerT> void SpgOCPsolver<dynT,costT,pconT,mayerT>::evalEconNorm(){

    double aux,*e,*end;

    end = econ+dL;

    e = econ;
    nec2 = 0.;
    while (e != end){
        aux = *e++; 
        nec2 += aux*aux;
    }

    nec = sqrt(nec2);
}    

/**
* "Real-time" primal loop stopped after fixed amount of iterations
*/
template <class dynT,class costT,class pconT,class mayerT> void SpgOCPsolver<dynT,costT,pconT,mayerT>::primalRealTime(){

    int i,it,bck;
    double ssz,maxpg,lmb;
    double az,dd; 
    double val,valNex,aux,val_max;
    double rhs,lhs,specNum,specDen;

    /* Initialize history array */
    for (i=0;i<histLenSpg;++i)
        *(vALhist+i) = NEGINF;

    /* Evaluate vAL */
    val = shooter.shootStatesFromIterate(halfRhoIni); 

    /* Store initial AL value into history array */
    *vALhist = val;

    /* Evaluate gAL */
    shooter.shootAdjointsFromIterate();

    /* Compute initial step-size */
    maxpg = -1.;
    for (i=0;i<dZ;++i){
        az = *(z+i);
        dd = ABS(PROJ(az-*(gAL+i),*(zmin+i),*(zmax+i))-az);
        maxpg = MAX(dd,maxpg);
    }

    if (maxpg<=POSZER){
        std::cout<<"Reached critical point, no need to proceed in SPG."<<std::endl;
        return;
    }

    memcpy(dffz,z,dZBL);
    memcpy(dffGal,gAL,dZBL);

    /**
    * SPG loop
    */ 
    ssz = 1./maxpg;
    it = 0;
    while (1){
//        std::cout<<"SPG iter "<<it<<", val= "<<val<<std::endl
        if (it>=maxPit)
            break;

        /* Find max. vAL over history */
        val_max = *vALhist;
        for (i=1;i<histLenSpg;++i)
            val_max = MAX(val_max,*(vALhist+i));

        /* Projected non-monotone line-search */
        lmb = ssz;
        bck = 0;
        while (1){

            /* Compute candidate */
            for (i=0;i<dZ;++i)
                *(znex+i) = PROJ(*(z+i)-*(gAL+i)*lmb,*(zmin+i),*(zmax+i));                

            /* Evaluate left handside*/
            lhs = shooter.shootStatesFromCandidate(halfRhoIni);

            /* Evaluate right handside */
            aux = 0.;
            for (i=0;i<dZ;++i)
                aux += *(gAL+i)*(*(znex+i)-*(z+i));
            aux *= gamSpg;
            rhs = val_max+aux;

            /* Check sufficient decrease */
            if (lhs<=rhs){
                val = lhs;
                /* Update primal iterate */
                memcpy(z,znex,dZBL);
                /* Update difference between primal iterates */
                for (i=0;i<dZ;++i)
                    *(dffz+i) -= *(z+i);
                break;
            }

            /* Check max. bck iters */
            if (bck>=maxBckIt){
                std::cerr<<"ERROR<SpgOCPsolver>: max. # iterations reached in non-monotone line-search."<<std::endl;
                return;
            }

            lmb *= 0.5*(sigTwoSpg+sigOneSpg);
            ++bck;
        }

         /* Update gradient */
        shooter.shootAdjointsFromIterate();
        for (i=0;i<dZ;++i)  
            *(dffGal+i) -= *(gAL+i);

        /* Spectral step-size computation */
        specDen = 0.;
        for (i=0;i<dZ;++i)
            specDen += *(dffz+i)**(dffGal+i);

        switch(specDen>POSZER){

            case 0:     
                ssz = amaxSpg;
                break;

            case 1:
                specNum = 0.;
                for (i=0;i<dZ;++i){
                    aux = *(dffz+i);
                    specNum += aux*aux;
                }
            //    std::cout<<"Spectral step-size."<<std::endl;
                ssz = MIN(amaxSpg,MAX(aminSpg,specNum/specDen));
                break;

            default:
                std::cerr<<"ERROR<SpgOCPsolver>: problem in spectral step-size computation."<<std::endl;
                return;
        }

        /* Shift history */
        memcpy(vALhist+1,vALhist,(histLenSpg-1)*SZDBL);
        *vALhist = val;


        /* Go to next iteration */
        memcpy(dffz,z,dZBL);
        memcpy(dffGal,gAL,dZBL);
        ++it;
    }
}

/**
* "Real-time" dual loop 
*/
template <class dynT,class costT,class pconT,class mayerT> void SpgOCPsolver<dynT,costT,pconT,mayerT>::solveRealTime(){

    int i;
    double *mm,*ee;
   
//    std::cout<<">>>> Real-time AL loop <<<<"<<std::endl;

    /* Fixed number of primal SPG iterations */
    primalRealTime();

    /* Evaluate 2-norm of equality constraints */
    evalEconNorm();

    /* Evaluate KKT satisfaction of augmented Lagrangian problem */
    evalKKTprim();

    std::cout<<"kkt="<<kkt<<std::endl;
    std::cout<<"nec="<<nec<<std::endl;

    /* One dual update */
    mm = mu;
    ee = econ;
    for (i=0;i<dL;++i)
        *mm++ += 2*halfRhoIni*(*ee++);

    /* Extract optimal control input */
    memcpy(uopt,z+dx,dubl);

}

/**
* Homotopy loop 
*/
template <class dynT,class costT,class pconT,class mayerT> void SpgOCPsolver<dynT,costT,pconT,mayerT>::solveHomotopy(){
    
    
    


    
    
    
    
    
    
}







