#include <iostream>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>

#include "itikBanks.h"
#include "newtonPoincare.h"
#include "fixPointsInSpace.h"



int main()
{
	//      Vector3d Fx(3);  //Vector of solutions , each solution is a Vector3d
	//      Fx << 0.8, 0.02, 0.03;
	//      Matrix3d I3=  Matrix3d::Identity(3, 3);
	//
	//      double tau=10,d=0.5;
	//      Matrix3d T = DFitikBanks(toStdVectorD(Fx), tau);
	//      T = T - I3;
	//      double firstDet = T.determinant();
	//
	//      cout<< firstDet;
	itikBanks fun(PARAMETERS);
	Matrix3d res = DFode_aprox(fun, stateType{ 0,1,0.06 }, 10, 0.5);

	std::cout << res;
	std::cin.get();
}