#ifndef ITIKBANKS_H
#define ITIKBANKS_H


#include "all.h"
//Itikbanks model from Wei2014
struct itikBanks{
	//The parameters
	double N;

	double a12,a21,a13,a31;//1,1.5,2.5,0.2,0.6,4.5,0.5,1
	double r2,r3,d3,k3;

itikBanks(std::vector<double> parameters){
	a12=parameters[0];a21=parameters[1];a13=parameters[2];a31=parameters[3];
	r2=parameters[4];r3=parameters[5];d3=parameters[6];k3=parameters[7];	
}	
void operator() (const stateType &x, stateType &dxdt,double t)
    {	
		dxdt[0] = x[0]*(1-x[0])-a12*x[0]*x[1]-a13*x[0]*x[2];
	    dxdt[1] = r2*x[1]*(1-x[1])-a21*x[0]*x[1];
	    dxdt[2] = (r3*x[0]*x[2])/(k3+x[0]) - a31*x[0]*x[2]-d3*x[2];
	}
};

struct itikBanksJacobian
{

	double a12,a21,a13,a31;//1,1.5,2.5,0.2,0.6,4.5,0.5,1
	double r2,r3,d3,k3;
	std::vector<stateType> state;
    std::vector<double> times;
	boost::math::barycentric_rational<double> T;
    boost::math::barycentric_rational<double> H;
    boost::math::barycentric_rational<double> E;

    itikBanksJacobian(std::vector<double> param,std::vector<stateType> xx, std::vector<double> tt):
    state(xx), times(tt), T(times.data(), state[0].data(), times.size()),
    H(times.data(), state[1].data(), times.size()),
    E(times.data(), state[2].data(), times.size()) {};

    void operator() (const stateType &x, stateType &dxdt,double t){

	    dxdt[0]= 1 - E(t)*a13 - H(t)*a12 - 2*T(t);
	    dxdt[1]= -T(t)*a12;
	    dxdt[2]= -T(t)*a13;
	    dxdt[3]= -H(t)*a21;
	    dxdt[4]= - T(t)*a21 - H(t)*r2 - r2*(H(t) - 1);
	    dxdt[5]= 0;
	    dxdt[6]= (E(t)*r3)/(T(t) + k3) - E(t)*a31 - pow((x[2]*T(t)*r3)/(T(t) + k3),2);
	    dxdt[7]= 0;
	    dxdt[8]= (T(t)*r3)/(T(t) + k3) - T(t)*a31;
 
    }


};

#endif