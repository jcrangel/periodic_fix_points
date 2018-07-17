#ifndef AL21MANAGER_H
#define AL21MANAGER_H
/**
This is acts as a interface for algorithm 21 for managing data, etc.


**/
#include "itikbanks.h"
#include "fixpointsinspace.h"
#include <fstream>
#include <iomanip>

void runAl21() {

	double xmin = 0; double xmax = 1;
	double ymin = 0; double ymax = 1;
	double zmin = 0; double zmax = 1;
	double M=10;
	//std::cin >> xmin >> xmax;
	//std::cin >> ymin >> ymax;
	//std::cin >> zmin >> zmax;
	//std::cin >> M;

	// { 1,2.5,0.5,1.5,4.5,1,0.2,0.5 }
	//std::vector<double> PARAMETERS = readParameters();
	std::vector<double> PARAMETERS = { 1,2.5,0.5,1.5,4.5,1,0.2,0.5 };
	double tau=10, d=0.5;
	//std::cin >> tau >> d;
	itikBanks fun(PARAMETERS);
	std::vector<fixPoint> S;


	S = al21(xmin, xmax, ymin, ymax, zmin, zmax, M, M, M, fun, tau, d);
	//std::ofstream points("points.txt", std::ios::out | std::ios::trunc);
	
	LogAndStdout lcout("points_log.txt");
	lcout << "========================================================\n";
	lcout << "Working with\n";
	lcout << "X=[" << xmin << " , " << xmax << " ]\n";
	lcout << "Y=[" << ymin << " , " << ymax << " ]\n";
	lcout << "Z=[" << zmin << " , " << zmax << " ]\n";
	lcout << "M: " << M << "tau: "<< tau << "d: "<<d << "\n";
	lcout << "Parameters:\n ";
	lcout << "(a12,a13,r2,a21,r3,k3,a31,d3)=(";
	for (double i : PARAMETERS) lcout << i << ",";
	lcout << ")\nPoints Founded:\n";
	lcout << std::fixed << std::setprecision(6);
	for (fixPoint i : S) {
		lcout << i.solution[0] << "\t" <<
			i.solution[1] << "\t" << i.solution[2] << "\t" << i.stability << "\n";
	}
	std::cin.get();
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





#endif