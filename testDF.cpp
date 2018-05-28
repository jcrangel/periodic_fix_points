#include <iostream>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>

#include "all.h"
#include "itikBanks.h"
#include "newtonPoincare.h"
#include "fixPointsInSpace.h"

using namespace Eigen;
using namespace std;
using namespace boost::numeric::odeint;

int main()
{
std::vector<Vector3d> Fx;  //Vector of solutions , each solution is a Vector3d
Matrix3d I3=  Matrix3d::Identity(3, 3);							
double xi=0.1,yi=0.1,zi=0.1; 									

Vector3d a;
a[0]=1;
a[1]=2;
a[2]=3; 
Fx.push_back(a);

double tau=10,d=0.5;

double firstDet = (DFitikBanks( std::vector<double> {xi,yi,zi}, Fx[0], tau, d ) - I3 ).determinant();

cout<< firstDet;
}