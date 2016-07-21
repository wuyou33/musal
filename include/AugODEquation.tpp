//
//  AugODEquation.cpp
//  MusAL
//
//  Created by Jean on 3/12/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//
#include "macros.h"
#include "AugODEquation.h"

/**
** Default constructor 
**/
template <class dynT,class cosT> AugODEquation<dynT,cosT>::AugODEquation(){

	/* Set dimensions */
    ds = ode.getDx();
    dy = cost.getDy();
    du = ode.getDu();
    dx = ds+1;
    da = dx+du;
    
    /* Set pointers to jacobians transpose */
    dynJxT = ode.getJacXt();
    dynJuT = ode.getJacUt();
    
    /* Set pointers to gradients */
    cosGx = cost.getGx();
    cosGu = cost.getGu();

    Dq = NULL;
}

/**
** Copy constructor
**/
template <class dynT,class cosT> AugODEquation<dynT,cosT>::AugODEquation(const AugODEquation<dynT,cosT>& obj){

    
    ode = obj.ode;
    cost = obj.cost;

    dx = obj.dx;
    da = obj.da;
    du = obj.du;
    ds = obj.ds;
    dy = obj.dy;

    yref = obj.yref;
    uref = obj.uref;

    dynJxT = ode.getJacXt();
    cosGx = cost.getGx();

    dynJuT = ode.getJacUt();
    cosGu = cost.getGu();

    Dq = obj.Dq;
}

/** 
** Assignment operator 
**/
template <class dynT,class cosT> AugODEquation<dynT,cosT>& AugODEquation<dynT,cosT>::operator= (const AugODEquation<dynT,cosT>& obj){
    

    if (this != &obj){

        ode = obj.ode;
        cost = obj.cost;

        dx = obj.dx;
        da = obj.da;
        du = obj.du;
        ds = obj.ds;
        dy = obj.dy;

        yref = obj.yref;
        uref = obj.uref;

        dynJxT = ode.getJacXt();
        cosGx = cost.getGx();

        dynJuT = ode.getJacUt();
        cosGu = cost.getGu();

        Dq = obj.Dq;
    }
    
    return *this;

}

/** 
** Destructor 
**/
template <class dynT,class cosT> AugODEquation<dynT,cosT>::~AugODEquation(){
}

/**
** Evaluate augmented dynamics rhs 
**/
template <class dynT,class cosT> void AugODEquation<dynT,cosT>::rhs(double* xrhs,double* x,double* u){

    /* Dynamics */
    ode.rhs(xrhs,x,u);

    /* Stage-cost */ 
    cost.evalStage(xrhs+ds,x,u,yref,uref); 
}

/**
** Implement augmented adjoint dynamics:
**	- arhs: right handside 
** 	- x: state
**	- u: input
**	- adj: augmented adjoint
**/
template <class dynT,class cosT> void AugODEquation<dynT,cosT>::adjRhs(double* arhs,double* x,double* u,double* adj){
    
    int i,j,jj;
    double aux,adjc;
    double *prd;

    /* Eval. state-jacobian transpose */
    ode.jxRhsT(x,u);
    /* Eval. input-jacobian transpose */
    ode.juRhsT(x,u);
    /* Eval. state-gradient */
    cost.evalGradXstage(x,u,yref,uref);
    /* Eval. input-gradient */
    cost.evalGradUstage(x,u,yref,uref);
    
    prd = arhs;

    /* Cost component of adjoint */
    adjc = *(adj+ds);

    /**
    * State-adjoint 
    */
    jj = 0;
    for (i=0;i<ds;++i){
        /* State-jacobian transpose */
        aux = 0.;
        for (j=0;j<ds;++j){
            aux -= *(dynJxT+jj)**(adj+j);
            ++jj;
        }
        /* Stage-cost state gradient */
        aux -= *(cosGx+i)*adjc;
        /* Gather product */
        *prd++ = aux;
    }    
    /* Cost part of output adjoint vector */
    *prd++ = 0.;

    /**
    * Input-adjoint 
    */
    jj = 0;
    for (i=0;i<du;++i){
        /* Input-jacobian transpose */
        aux = 0.;
        for (j=0;j<ds;++j){
            aux -= *(dynJuT+jj)**(adj+j);
            ++jj;
        }
        /* Stage-cost input gradient */
        aux -= *(cosGu+i)*adjc;       
        /* Gather product */
#if UNSCAL_VAR
        *prd++ = aux;
#endif
#if SCAL_VAR
        *prd++ = *(Dq+i)*aux;
#endif
    }
}
