/*
In this file are the methods & classes for finding fixed points on space 
That is the algorithms 2.1 and 2.2 of Wei 2014
*/

#ifndef FIXPOINTSINSPACE_H
#define FIXPOINTSINSPACE_H

#include "all.h"

typedef std::vector<fixPoint> fixPoints;

//Return a set(vector) of fixed points
template <class T>
fixPoints al21(double xmax,double ymax, double zmax, double M,double N, double L,T &functionName,std::vector<double> parameters,double tau,double d)
{
fixPoints S;
double xstep=xmax/M;
double ystep=ymax/N;
double zstep=zmax/L;
double xmin=0;
double ymin=0;
double zmin=0;

fixPoints O;

for(int i = xmin ;i <= xmax ; i += xstep)
	for(int j = ymin ;j <= ymax ; j += ystep)
		for(int k = zmin ;k <= zmax ; k += zstep){
			O =al22(i,i+xstep,j,j+ystep,k,k+zstep,O,
				functionName,parameters,tau,d,6);	
		}
S = O;

return S;

}

//%Evaluate F(x,y,z) 
template <class T>
Vector3d evalFunInLast(T &functionName,
			   std::vector<double> initialCondition,double tau,double d){

	std::vector<double > xp={0, tau};
	initialCondition[CONTROL_POS] = initialCondition[CONTROL_POS] + d ;
	std::vector<stateType> V;
	std::vector<double> t;

    steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6), 
            functionName, initialCondition, xp[0], xp[1], 0.001);
//initialCondition has the last
Vector3d v(initialCondition.data());

return v;
}

struct pointxyz{
	double x;	
	double y;	
	double z;
	pointxyz(double x_,double y_,double z_):
	x(x_),y(y_),z(z_){}	
}


template <class T>
fixPoints al22(double xi,double xf,double yi,double yf,double zi,double zf,
			    fixPoints S,T &functionName,double tau,double d,int deepness)
{

if (deepness <= 0)
    return S;

//MatrixXd Fx = MatrixXd::Zero(8,3); // #_#
Vector3d Fx= evalFunInLast(functionName,[xi,yi,zi],tau,d);
pointxyz first(Fx[0]-xi,Fx[1]-yi,Fx[2]-zi);
bool xgood=false;
bool ygood=false;
bool zgood=false;
int n=1;

for(double i: std::vector<double> {xi,xf})
    for(double j: std::vector<double> {yi,yf})
        for (double k: std::vector<double> {zi,zf})
        {
        	//The first element is already computed
            if(n==1){
                n=n+1;
                continue;
            }
            //Step1

        }


}
#endif





















