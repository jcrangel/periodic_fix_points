#ifndef AL21MANAGER_H
#define AL21MANAGER_H
/**
This is acts as a interface for algorithm 21 for managing data, etc.


**/
#include "itikbanks.h"
#include "fixpointsinspace.h"
#include <fstream>
#include <iomanip>
#include <vector>

void runAl21() {

	Doub xmin = 0; Doub xmax = 1;
	Doub M=10;
	//std::cin >> xmin >> xmax;
	//std::cin >> ymin >> ymax;
	//std::cin >> zmin >> zmax;
	//std::cin >> M;

	// { 1,2.5,0.5,1.5,4.5,1,0.2,0.5 }
	//std::vector<Doub> PARAMETERS = readParameters();
	std::vector<Doub> PARAMETERS = { 1,2.5,0.5,1.5,4.5,1,0.2,0.5 };
	Doub tau=10, d=0.5;
	//std::cin >> tau >> d;
	ItikBanks fun(PARAMETERS);
	std::vector<FixPoint> S;
	std::vector<double> v;



	S = al21(xmin, xmax, M, fun, tau, d);
	//std::ofstream points("points.txt", std::ios::out | std::ios::trunc);

	LogAndStdout lcout("points_log.txt");
	lcout << "========================================================\n";
	lcout << "Working with\n";
	lcout << "X=[" << xmin << " , " << xmax << " ]\n";
	lcout << fun.getSystemSize() <<" intervals \n";
	lcout << "M: " << M << "tau: "<< tau << "d: "<<d << "\n";
	lcout << "Parameters:\n ";
	lcout << "(a12,a13,r2,a21,r3,k3,a31,d3)=(";
	for (Doub i : PARAMETERS) lcout << i << ",";
	lcout << ")\nPoints Founded:\n";
	lcout << std::fixed << std::setprecision(6);
	for (FixPoint i : S) {
		lcout << i.solution[0] << "\t" <<
			i.solution[1] << "\t" << i.solution[2] << "\t" << i.stability << "\n";
	}
	/*std::cin.get();*/
}


#endif
