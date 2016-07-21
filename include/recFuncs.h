#ifndef RECFUNCS_H
#define RECFUNCS_H

/**
** Some recursive methods for compile-time evaluation (C++11) 
**/
/* constexpr recursive computation of x^j */
constexpr double recPow(double x,int j){

    return (j==0) ? 1. : x*recPow(x,j-1);
}

/* Meta-template power */
template <int N> inline double Power(double x){
	
	return x*Power<N-1>(x);
}

template <> inline double Power<0>(double x){

	return 1.0;
}

/* constexpr recursive computation of n! */
constexpr double recFact(int n){

    return (n==0) ? 1. : n*recFact(n-1);
}

/* Meta-template factorial */
template <int N> inline double Factorial(){

	return N*Factorial<N-1>();
}	

template <> inline double Factorial<0>(){

	return 1.0;
}

/* constexpr recursive computation of cos(x) (j-terms of the power series)*/
constexpr double recCos(double x,int j){

    return (j==0) ? 1 : recPow(-1.,j)*recPow(x,2*j)/recFact(2*j)+recCos(x,j-1);
}

/* Meta-template cosine */
template <int N> inline double Cosine(double x){
	return Cosine<N-1>(x)+(N%2 ? -1. : 1.)*Power<2*N>(x)/Factorial<2*N>();
}

template <> inline double Cosine<0>(double x){
	
	return 1.0;
}

#endif










