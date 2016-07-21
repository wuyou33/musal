//
//  MyTrajectory.cpp
//  MusAL
//
//  Created by Jean on 4/23/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//
#include <string.h>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <assert.h>

#include "MyTrajectory.h"
#include "macros.h"


/**
 * Default constructor
 */
MyTrajectory::MyTrajectory(){
    
    yRef = new double[DY*LREFMAX];
    memset(yRef,-1,DY*LREFMAX*SZDBL);

    uRef = new double[DU*LREFMAX];
    memset(uRef,-1,DU*LREFMAX*SZDBL);

    ti = -1.;
    tf = -1.;
    ts = -1.;
    Ns = -1;
}


/**
* Destructor
*/
MyTrajectory::~MyTrajectory(){

    if (yRef != NULL)
    	delete [] yRef;
    
    if (uRef != NULL)
    	delete [] uRef;
}


/**
 * Copy constructor
 */
MyTrajectory::MyTrajectory(const MyTrajectory& obj){

	yRef = new double[DY*LREFMAX];
	std::copy(obj.yRef,obj.yRef+DY*LREFMAX,yRef);

	uRef = new double[DU*LREFMAX];
	std::copy(obj.uRef,obj.uRef+DU*LREFMAX,uRef);

	ti = obj.ti;
	tf = obj.tf;
	ts = obj.ts;	
	Ns = obj.Ns;
}


/**
 * Assignment operator
 */
MyTrajectory& MyTrajectory::operator= (const MyTrajectory& obj){

	double *yRef_,*uRef_;

	if (this != &obj){

		if (yRef != NULL)
			delete [] yRef;

		yRef_ = new double[DY*LREFMAX];
		std::copy(obj.yRef,obj.yRef+DY*LREFMAX,yRef_);
		yRef = yRef_;

		if (uRef != NULL)
			delete [] uRef;

		uRef_ = new double[DU*LREFMAX];
		std::copy(obj.uRef,obj.uRef+DU*LREFMAX,uRef_);
		uRef = uRef_;

		ti = obj.ti;
		tf = obj.tf;
		ts = obj.ts;
		Ns = obj.Ns;
	}

	return *this;
}


/**
* USER-DEFINED: compute output reference point at time instant t
*/
void MyTrajectory::computeOutputPoint(double t,double* y){

	/* USER-DEFINED METHOD */

	/* Walid's crazyflie */
/*	*y = 0.04*t*cos(M_PI*t);
	*(y+1) = 0.04*t*sin(M_PI*t);
	*(y+2) = 0.1*t;
	*(y+3) = 0.;
	*(y+4) = 0.;
	*(y+5) = 0.;
	*(y+6) = 0.;
	*(y+7) = 0.;
	*(y+8) = 0.;
	*(y+9) = 0.;
	*(y+10) = 0.;
	*(y+11) = 0.; */

	/* Inverted pendulum */
/*	*y = 0.;
	*(y+1) = 0.;
	*(y+2) = M_PI;
	*(y+3) = 0.; */

	/* Inverted pendulum with polynomial dynamics */
/*	*y = 0.;
	*(y+1) = 0.;
	*(y+2) = -1.;
	*(y+3) = 0.;
	*(y+4) = 0.; */

	/* Crane */
/*	*y = 0.1;
	*(y+1) = 0.;
	*(y+2) = 0.;
	*(y+3) = 0.; */

	/* DC motor */
/*	*y = 2.0; */

	/* Unicycle */
//	*y = 0.7*cos(0.5*t);
//	*(y+1) = 0.7*sin(0.5*t);
//	*(y+2) = 0.5*M_PI+0.5*t;
}


/**
* USER-DEFINED: compute input reference point at time instant t
*/
void MyTrajectory::computeInputPoint(double t,double* u){

	/* USER-DEFINED METHOD */

	/* Crazyflie */
/*	*u = 35728.;
	*(u+1) = 0.;
	*(u+2) = 0.;
	*(u+3) = 0.; */

	/* Inverted pendulum */
//	*u = 0.;
}


/**
* Generate output reference in [ti,tf] sampled at ts
*/
void MyTrajectory::computeOutputReference(){

	assert(ti<tf);
	assert(Ns>=2);

	int isInt,n;
	double t,aux,auxf;
	double *yy;

	if (Ns>LREFMAX){
		std::cerr<<"ERROR<computeOutputReference>: Too many reference points, consider increasing LREFMAX in macros.h ."<<std::endl;
		return;
	}
	ts = (tf-ti)/(Ns-1);
	n = 0;
	t = ti;
	yy = yRef;
	while (1){
		if (n>=Ns)
			break;
		computeOutputPoint(t,yy);
		yy += dy;
		t += ts;
		++n;
	}
}

/**
* Generate input reference in [ti,tf] sampled at ts
*/
void MyTrajectory::computeInputReference(){

	assert(ti<tf);
	assert(Ns>=2);

	int isInt,n;
	double t;
	double *uu;

	if (Ns>LREFMAX){
		std::cerr<<"ERROR<computeInputReference>: Too many reference points, consider increasing LREFMAX in macros.h"<<std::endl;
		return;
	}
	ts = (tf-ti)/(Ns-1);
	n = 0;
	t = ti;
	uu = uRef;
	while (1){
		if (n>=Ns)
			break;
		computeInputPoint(t,uu);
		uu += du;
		t += ts;
		++n;	
	}
}

/** 
* Display references
*/
void MyTrajectory::displayReference(){

	int i,n;
	double *yy,*uu;

	/* Output reference */
	std::cout<<"<======================>"<<std::endl;
	std::cout<<"Output reference:"<<std::endl;
	n=0;
	yy = yRef;
	while (1){
		if (n>Ns)
			break;
		std::cout<<"n="<<n<<std::endl;
		for (i=0;i<dy;++i)
			std::cout<<*yy++<<std::endl;
		++n;
	}
	/* Input reference */
	std::cout<<"Input reference:"<<std::endl;
	n = 0;
	uu = uRef;
	while (1){
		if (n>Ns)
			break;
		std::cout<<"n="<<std::endl;
		for (i=0;i<du;++i)
			std::cout<<*uu++<<std::endl;
		++n;
	}
	std::cout<<"<=======================>"<<std::endl;
}