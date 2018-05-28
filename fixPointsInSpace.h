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

    for(int i = xmin ; i <= xmax ; i += xstep)
        for(int j = ymin ; j <= ymax ; j += ystep)
            for(int k = zmin ; k <= zmax ; k += zstep)
            {
                O =al22(i,i+xstep,j,j+ystep,k,k+zstep,O,
                        functionName,parameters,tau,d,6);
            }
    S = O;

    return S;

}

//%Evaluate F(x,y,z)
template <class T>
Vector3d evalFunInLast(T &functionName,
                       std::vector<double> initialCondition,double tau,double d)
{
    typedef runge_kutta_cash_karp54<stateType> error_stepper_type;
    std::vector<double > xp = {0, tau};
    initialCondition[CONTROL_POS] = initialCondition[CONTROL_POS] + d ;
    std::vector<stateType> V;
    std::vector<double> t;

    size_t steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6),
                               functionName, initialCondition, xp[0], xp[1], 0.001);
//initialCondition has the last
    Vector3d v(initialCondition.data());

    return v;
}

struct pointxyz
{
    double x;
    double y;
    double z;
    pointxyz(double x_,double y_,double z_):
        x(x_),y(y_),z(z_) {}
};


template <class T>
fixPoints al22(double xi,double xf,double yi,double yf,double zi,double zf,
               fixPoints S,T &functionName,double tau,double d,int deepness)
{

    if (deepness <= 0)
        return S;

//MatrixXd Fx = MatrixXd::Zero(8,3); // #_#
    std::vector<Vector3d> Fx;  //Vector of solutions , each solution is a Vector3d
//First solution
    Fx[0] = evalFunInLast(functionName,std::vector<double> {xi,yi,zi},tau,d);
    pointxyz first(Fx[0][0]-xi,Fx[0][1]-yi,Fx[0][2]-zi);
    bool xgood=false;
    bool ygood=false;
    bool zgood=false;
    int n=0;
//pointxyz fminus;
    double mulx;
    double muly;
    double mulz;
//%this is wrong we need to evalute this on the middle points see algorithm
//%2.2 step1
    for(double i: std::vector<double> {xi,xf})
        for(double j: std::vector<double> {yi,yf})
            for (double k: std::vector<double> {zi,zf})
            {
                //The first element is already computed
                if(n==0)
                {
                    n=n+1;
                    continue;
                }
                
                if(DEBUG){
                std::cout << "Step 1"<< std::endl ;
                }

                //Step1
                Fx[n] = evalFunInLast(functionName,std::vector<double> {xi,yi,zi},tau,d);
                //Not effice to create every loop?   //#_#
                pointxyz fminus(Fx[n][0]-xi,Fx[n][1]-yi,Fx[n][2]-zi);

                if(DEBUG){
                std::cout << "Step 2"<< std::endl ;
                }
//            %Step 2 Check if first.x  differs at sign with fminus.x
//            %the evaluation. Compares the first one with everything els
                mulx=fminus.x*first.x;
                muly=fminus.y*first.y;
                mulz=fminus.z*first.z;

            if(mulx <= 0) xgood=true;
            if(muly <= 0) ygood=true;
            if(mulz <= 0) zgood=true;


            //%if actually changes from 0, we take it as a sing change.
            if(equal(first.x,0) && !equal(fminus.x,0)) xgood=true;
            if(equal(first.y,0) && !equal(fminus.y,0)) ygood=true;
            if(equal(first.z,0) && !equal(fminus.z,0)) zgood=true; 
            
            n=n+1;
            }

//%If every element of Fx(xi,yi,zi)-xi has the same sign returns.(The same for %y,z)
if ( !xgood || !ygood || !zgood ) 
    return S;

if(DEBUG){
	std::cout << "Step 3"<< std::endl ;
}	
bool signChange=false;
	
	
Matrix3d I3=  Matrix3d::Identity(3, 3);																
double firstDet = (DFitikBanks( std::vector<double> {xi,yi,zi}, Fx[0],tau,d) - I3 ).determinant();
n=0;
double Det,mult;
    for(double i: std::vector<double> {xi,xf})
        for(double j: std::vector<double> {yi,yf})
            for (double k: std::vector<double> {zi,zf})
            {
                //The first element is already computed
                if(n==0)
                {
                    n=n+1;
                    continue;
                }
//            Det = det(DF(functionName,parameters,[i,j,k],Fx(n,:),tau,d)-[1,0,0;0,1,0;0,0,1]);
            //Step 4
            if(DEBUG)
            std::cout << "step 4" <<std::endl;
            
            mult = Det * firstDet;

            }


}
#endif





















