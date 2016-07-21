//
//  MyMayer.cpp
//  MusAL
//
//  Created by Jean on 3/14/15.
//  Copyright (c) 2015 Jean. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <string.h>

#include "MyMayer.h"
#include "macros.h"

/**
** Default constructor 
**/
MyMayer::MyMayer(){

	memset(gX,0,dx*SZDBL);
}

/**
** Destructor 
**/
MyMayer::~MyMayer(){}

/**
** Compute Mayer term 
**/
double MyMayer::eval(double* x){
    
    double may;
    
    /* USER-DEFINED */
#if SCAL_VAR
    /* State scaling */
    int i;
    for (i=0;i<dx;++i)
    	*(xsc+i) = *(cs+i)+*(Ds+i)**(x+i);
#endif


    may = 0.;
    
    return may;
}

/**
** Compute gradient of Mayer term
**/
void MyMayer::grad(double* x){
	/* USER-DEFINED */
	int i;

	for (i=0;i<dx;++i){
#if SCAL_VAR
		/* State scaling */
		*(xsc+i) = *(cs+i)+*(Ds+i)**(x+i);
#endif
		*(gX+i) = 0.;
	}
}