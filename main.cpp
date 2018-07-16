/*
To do:
1. Define StateTye as Vector3d so is congruent with all the program.
    -- Create a function that integrate inside use std::vector outside return Vector3d
    --
2. Theres lots code chucks were a ode is solve it, but only the last value is needed. And save all
the state values is unnecesary

Symbols:
        #_# means that part of the should be rewritten to support systems of several size
        !@_@ incongruence with Vector and std::vector
        *_* unnecesary state saving
Note:
 The algorithm finds the same points several times. How this can be avoided?
 -Why to disthinguis beetween std::vector and boost::vector how about just use boost vector
  for everything

  Why is not working?
  -- DFode and DF, DF calculates the jacobian only using the last data. I can't find any reference were
		this should be done like that
	--%algorithm 2.2 line 88 step1




*/

//#include "all.h"
#include "itikBanks.h"
#include "newtonPoincare.h"
#include "fixPointsInSpace.h"
#include <fstream>
#include <iomanip>
//typedef std::vector<fixPoint> fixPoints;
extern std::vector<double> PARAMETERS;
int main(){

		double xmin; double xmax;
		double ymin; double ymax;
		double zmin; double zmax;
		double M;
		std::cin >> xmin >> xmax;
		std::cin >> ymin >> ymax;
		std::cin >> zmin >> zmax;
		std::cin >> M;

		// { 1,2.5,0.5,1.5,4.5,1,0.2,0.5 }
		std::vector<double> PARAMETERS = readParameters();
		double tau, d;
		std::cin >> tau>> d;
		itikBanks fun(PARAMETERS);
		std::vector<fixPoint> S;


		S = al21(xmin,xmax,ymin,ymax,zmin,zmax,M,M,M, fun, tau, d);
		std::ofstream points("points.txt",std::ios::out | std::ios::trunc);
		points << std::fixed << std::setprecision(6);
	for (fixPoint i : S) {
		points << i.solution[0] << "\t" <<
			i.solution[1] << "\t" << i.solution[2] << "\t" << i.stability<<std::endl;
	}


	//Test for function pointIsInSet
	//std::vector<fixPoint> S;
	//S.push_back(fixPoint(true, true, Vector3d(3.26603e-23, 6.85002e-31, 0.00339196) ));
	//S.push_back(fixPoint(true, true, Vector3d(3.26603e-13, 6.85002e-21, 0.00339196) ));
	//S.push_back(fixPoint(true, true, Vector3d(0, 0, 0.00339197) ));
	//S.push_back(fixPoint(true, true, Vector3d(3.26603e-15, 6.85002e-23, 0.00339196) ));
	//S.push_back(fixPoint(true, true, Vector3d(0.328024, 3.37749e-24, 0.268567) ));
	////True Cases
	//std::cout << pointIsInSet(fixPoint(true, true, Vector3d(0, 0, 0.00339197)),S)<<std::endl;
	//std::cout << pointIsInSet(fixPoint(true, true, Vector3d(0.328024, 3.37749e-24, 0.268567)),S)<<std::endl;
	//std::cout << pointIsInSet(fixPoint(true, true, Vector3d(3.26603e-17, 6.85002e-17, 0.00339196)),S)<<std::endl;
	////False Cases
	//std::cout << pointIsInSet(fixPoint(true, true, Vector3d(0, 1, 0.00339197)),S)<<std::endl;
	//std::cout << pointIsInSet(fixPoint(true, true, Vector3d(0.002, 1, 0.00339197)),S)<<std::endl;
	//std::cout << pointIsInSet(fixPoint(true, true, Vector3d(3.26603e-17, 6.85002e-17, 0.00239196)), S) << std::endl;
	//std::cout<<"program finished!"<< std::endl;
	//std::cin.get();//VS windows only

	return 0;
}
