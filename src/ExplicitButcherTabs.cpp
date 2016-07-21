#include <string.h>
#include <iostream>

#include "macros.h"
#include "ExplicitButcherTabs.h"


/**
 * Default constructor
 */
ExplicitButcherTabs::ExplicitButcherTabs(){

	/* Set data for explicit trapezoidal */
	memset(A_erk_trap,0,4*SZDBL);
	A_erk_trap[2] = 1.;

	memset(b_erk_trap,0,2*SZDBL);
	b_erk_trap[0] = 0.5;
	b_erk_trap[1] = 0.5;

	memset(c_erk_trap,0,2*SZDBL);
	c_erk_trap[1] = 1.; 

	/* Set data for Runge */
    memset(A_erk_rung,0,4*SZDBL);
    A_erk_rung[2] = 0.5;

    memset(b_erk_rung,0,2*SZDBL);
    b_erk_rung[1] = 1.;

    memset(c_erk_rung,0,2*SZDBL);
    c_erk_rung[1] = 0.5;
    
    /* Set data for ERK4 */
	memset(A_erk_4,0,16*SZDBL);
	A_erk_4[4] = 0.5;
	A_erk_4[9] = 0.5;
	A_erk_4[14] = 1.;
    
	memset(b_erk_4,0,4*SZDBL);
	b_erk_4[0] = 1./6.;
	b_erk_4[1] = 2./6.;
	b_erk_4[2] = 2./6.;
	b_erk_4[3] = 1./6.;

	memset(c_erk_4,0,4*SZDBL);
	c_erk_4[1] = 0.5;
	c_erk_4[2] = 0.5;
	c_erk_4[3] = 1.;
}



/**
* Destructor
*/
ExplicitButcherTabs::~ExplicitButcherTabs(){

}