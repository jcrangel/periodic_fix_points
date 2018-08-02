#ifndef AL21MANAGER_H
#define AL21MANAGER_H
/**
This is acts as a interface for algorithm 21 for managing data, etc.


**/
#include "murinemodelv2.h"
#include "fixpointsinspace.h"
#include <fstream>
#include <iomanip>
#include <vector>

void runAl21() {

	Doub xmin = 0; Doub xmax = 1e12;
	Doub M=10;
	    Doub aH, MuH, rH, KH, aC, MuC, rC, KC, KT, MuD, rI, MuIC, MuI;
    Doub rT, aT, eT, hT, aTBeta, eTBeta, rTBeta, MuBeta, aGammaC;
    Doub MuGamma, gMl, aMlGamma, eMlGamma, MuMl;

	aH = 1e-4;
	MuH = 0.005;
	rH = 10e-2;
	KH = 1;
	aC = 1e-4;
	MuC = 0.01925;
	rC = 0.00004e-2;
	KC = 1;
	KT = 1e12;
	MuD = 0.009625;
	rI = 1e-2;
	MuIC = 1e-7;
	MuI = 1e-2;
	rT = 0.002;
	aT = 0.1136;
	eT = 50;
	hT = 5.2e5;
	aTBeta = 0.69;
	eTBeta = 1e4;
	rTBeta = 5.57e-6;
	MuBeta = 6.93;
	aGammaC = 1.02e-4;
	MuGamma = 0.102;
	gMl = 1.44;
	aMlGamma = 2.89;
	eMlGamma = 3.38e5;
	MuMl = 0.0144;


	std::vector<Doub> PARAMETERS = { aH, MuH, rH, KH, aC, MuC, rC, KC, KT, MuD, rI, MuIC, 
	MuI, rT, aT, eT, hT, aTBeta, eTBeta, rTBeta, MuBeta, aGammaC, MuGamma, gMl, aMlGamma, eMlGamma, MuMl};
	Doub ef =0.05;
	Doub tau=251.4, dose=140000;
	//std::cin >> tau >> d;
	MurineModelv2 fun(PARAMETERS);
	std::vector<FixPoint> S;
	std::vector<double> v;
	


	S = al21(xmin, xmax, M, fun, tau, dose);
	//std::ofstream points("points.txt", std::ios::out | std::ios::trunc);

	LogAndStdout lcout("points_log.txt");
	lcout << "========================================================\n";
	lcout << "Working with\n";
	lcout << "X=[" << xmin << " , " << xmax << " ]\n";
	lcout << fun.getSystemSize() <<" intervals \n";
	lcout << "M: " << M << "tau: "<< tau << "d: "<<dose << "\n";
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
