#ifndef AL21MANAGER_H
#define AL21MANAGER_H
/**
This is acts as a interface for algorithm 21 for managing data, etc.


**/
#include "castiglione.h"
#include "parameterswept.h"
#include <fstream>
#include <iomanip>
#include <vector>

void runAl21()
{
	Doub a0, b0, c0, d0, f0;
	Doub a1, b1, c1, d1, f1;
	Doub d2, f2, e2;
	Doub e3;
	Doub a4, c4, e4;

	a0 = 1e-4;
	b0 = 0.005;
	c0 = 10;
	d0 = 1e-2;
	f0 = 1;
	a1 = 1e-4;
	b1 = 0.005;
	c1 = 10;
	d1 = 1e-2;
	f1 = 1;
	d2 = 0.02;
	e2 = 0.1;
	f2 = 1;
	e3 = 0.1;
	a4 = 1e-2;
	c4 = 1e-7;
	e4 = 1e-2;

	std::vector<Doub> PARAMETERS = { a0, b0, c0, d0, f0, a1, b1, c1, d1, f1,
									 d2, f2, e2, e3, a4, c4, e4 };
	//Doub tau = 1000, dose = 0.5;
	//std::cin >> tau >> d;
	Castiglione fun(PARAMETERS);
	StateType x0{1e-4 / 0.005, 1e-4 / 0.005, 0.1, 0, 0};  

	//std::ofstream points("points.txt", std::ios::out | std::ios::trunc);
	Doub taui = 10;
	Doub tauf = 1200;
	Doub di = 0.01;
	Doub df = 1;
	setCriticalCells Cri;
	Doub tumorBound=0.8*f2;// percentage of carrying capacity
	Cri= al23_tau_d(taui, tauf, di, df, 10,10, fun, tumorBound, x0, fun.getTumorIndex());
	LogAndStdout lcout("points_log.txt");

	lcout << "========================================================\n";
	lcout << "Castiglione parameter swepting tau and d\n";
	lcout << "tau=[" << taui << "," << tauf << "]\n";
	lcout << "d=[" << di << "," << df << "]\n";
	lcout << "Parameters:\n ";
	lcout << "(a0, b0, c0, d0, f0, a1, b1, c1, d1, f1 d2, f2, e2, e3, a4, c4, e4 )=(";
	for (Doub i : PARAMETERS)
		lcout << i << ",";

	std::ofstream out("points.txt", std::ios::out );
	for (cell2D cell : Cri) {
		out << cell;
	}
	lcout << "Tumor upper bound: " << tumorBound << "\n";
	lcout << "\nCritical cell point data saved to points.txt\n";
	/*std::cin.get();*/
}

#endif
