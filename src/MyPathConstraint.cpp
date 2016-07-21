//
//  MyPathConstraint.cpp
//  MusAL
//
//  Created by Jean on 3/14/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#include <string.h>
#include <algorithm>
#include <iostream>

#include "MyPathConstraint.h"

/**
* Default constructor 
*/
MyPathConstraint::MyPathConstraint(){

    memset(jacXt,0,dx*dpc*SZDBL);
    memset(jacUt,0,du*dpc*SZDBL);
}


/** 
* Destructor 
*/
MyPathConstraint::~MyPathConstraint(){
	
}


/**
* Compute path-constraint g(x,u)+r (x:state, u:input, r:slack)
*/
void MyPathConstraint::eval(double* con,double* x,double* u,double* r){
    
    int i;
    
    /* USER-DEFINED */
    
    
    /* Add slack variable */
//    for (i=0;i<dpc;++i)
//    	*(con+i) += *(r+i);
    
}


/**
* Transpose of state-jacobian of path-constraint (ROW MAJOR)
*/
void MyPathConstraint::evalJacXt(double* x,double* u){
    
    /* User defined */
}


/**
* Tranpose of input-jacobian of path-constraint (ROW MAJOR)
*/
void MyPathConstraint::evalJacUt(double* x,double* u){
    
    /* User defined */
}


/**
* Evaluate product of vector against state-jacobian transpose 
*/
void MyPathConstraint::prdJacXt(double* vout,double* vin){

	int i,j,jj;
	double rw;

	jj = 0;
	for (i=0;i<dx;++i){
		rw = 0.;
		for (j=0;j<dpc;++j){
			rw += *(jacXt+jj)**(vin+j);
			++jj;
		}
		*(vout+i) += rw;
	}
	
}
 
    
/**
* Evaluate product of vector against input-jacobian transpose 
*/
void MyPathConstraint::prdJacUt(double* vout,double* vin){

	int i,j,jj;
	double rw;

	jj = 0;
	for (i=0;i<du;++i){
		rw = 0.;	
		for (j=0;j<dpc;++j){
			rw += *(jacUt+jj)**(vin+j);
			++jj;
		}
		*(vout+i) += rw;
	}

}
