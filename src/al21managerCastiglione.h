#ifndef AL21MANAGER_H
#define AL21MANAGER_H
/**
This is acts as a interface for algorithm 21 for managing data, etc.


**/
#include "castiglione.h"
#include "fixpointsinspace.h"
#include <fstream>
#include <iomanip>
#include <vector>

void runAl21()
{

	Doub xmin = 0;// 0 zero give division over zero..
	Doub xmax = 1;
	Doub M = 2;
	Doub a0, b0, c0, d0, f0;
	Doub a1, b1, c1, d1, f1;
	Doub d2, f2, e2;
	Doub e3;
	Doub a4, c4, e4;

	a0 = 10e-4;
	b0 = 0.005;
	c0 = 10;
	d0 = 10e-2;
	f0 = 1;
	a1 = 10e-4;
	b1 = 0.005;
	c1 = 10;
	d1 = 10e-2;
	f1 = 1;
	d2 = 0.02;
	e2 = 0.1;
	f2 = 1;
	e3 = 0.1;
	a4 = 10e-2;
	c4 = 10e-7;
	e4 = 10e-2;

	std::vector<Doub> PARAMETERS = { a0, b0, c0, d0, f0, a1, b1, c1, d1, f1,
									 d2, f2, e2, e3, a4, c4, e4 };
	Doub tau = 1000, dose = 1;
	//std::cin >> tau >> d;
	Castiglione fun(PARAMETERS);
	std::vector<FixPoint> S;
	std::vector<double> v;

	S = al21(xmin, xmax, M, fun, tau, dose);
	//std::ofstream points("points.txt", std::ios::out | std::ios::trunc);

	LogAndStdout lcout("points_log.txt");
	lcout << "========================================================\n";
	lcout << "Working with\n";
	lcout << "X=[" << xmin << " , " << xmax << " ]\n";
	lcout << fun.getSystemSize() << " intervals \n";
	lcout << "M: " << M << "tau: " << tau << "d: " << dose << "\n";
	lcout << "Parameters:\n ";
	lcout << "(a0, b0, c0, d0, f0, a1, b1, c1, d1, f1 d2, f2, e2, e3, a4, c4, e4 )=(";
	for (Doub i : PARAMETERS)
		lcout << i << ",";
	lcout << ")\nPoints Founded:\n";
	lcout << std::fixed << std::setprecision(6);
	for (FixPoint i : S)
	{
		lcout << i.solution[0] << "\t" << i.solution[1] << "\t" << i.solution[2] << "\t" << i.stability << "\n";
	}
	/*std::cin.get();*/
}

#endif
