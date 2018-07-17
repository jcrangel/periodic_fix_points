#ifndef ITIKBANKS_H
#define ITIKBANKS_H


#include "all.h"
//Itikbanks model from Wei2014
struct itikBanks{
	//The parameters
	double N;

	double a12,a21,a13,a31;//1,1.5,2.5,0.2,0.6,4.5,0.5,1
	double r2,r3,d3,k3;

itikBanks(double a12_,double a21_,double a13_,double a31_,
		      double r2_,double r3_,double d3_,double k3_){
	a12=a12_;a21=a21_;a13=a13_;a31=a31_;
	r2=r2_;r3=r3_;d3=d3_;k3=k3_;	
}	
void operator() (const stateType &x, stateType &dxdt,double t);
    {	
		dxdt[0] = x[0]*(1-x[0])-a12*x[0]*x[1]-a13*x[0]*x[2];
	    dxdt[1] = r2*x[1]*(1-x[1])-a21*x[0]*x[1];
	    dxdt[2] = (r3*x[0]*x[2])/(k3+x[0]) - a31*x[0]*x[2]-d3*x[2];
	}
//Returns the jacobian of the system evaluated on the actual values of the solution
matrix<double > Dfu();

};


#endif