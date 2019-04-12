#ifndef AL21MANAGERERANDI_H
#define AL21MANAGERERANDI_H
/**
This is acts as a interface for algorithm 21 for managing data, etc.


**/
#include "erandi_sys.h"
#include "fixpointsinspace.h"
#include <fstream>
#include <iomanip>
#include <vector>

void runAl21()
{
	Doub M = 2;
	Doub xmin = 1;// 0 zero give division over zero..
	Doub xmax = 1e9;
    Doub r,k,aT,aTBeta,eTBeta,muD,ra;
    Doub thetaD,re, muCa,thetaA,muCi,rTBeta;
    Doub muBeta;
    r = 0.001060289;
	k = 6.75e15;
	aT = 6e-11,
	aTBeta = 0.69,
	eTBeta = 1e4,
	muD = 0.009625,
	ra = 61;
    thetaD = 212;
	re = 6500;
	muCa = 0.01925;
	thetaA = 10;
	muCi = 0.007;
	rTBeta = 5.75e10-6;
    muBeta = 2.7 ;


	std::vector<Doub> PARAMETERS = {r,k,aT,aTBeta,eTBeta,muD,ra, thetaD,re, muCa,thetaA,muCi,rTBeta, muBeta };
	Doub tau = 168, dose = 1e7;
	//std::cin >> tau >> d;
	ErandiSys fun(PARAMETERS);
	std::vector<FixPoint> S;
	std::vector<double> v;

	S = al21(xmin, xmax, M, fun, tau, dose);
	//std::ofstream points("points.txt", std::ios::out | std::ios::trunc);

	LogAndStdout lcout("points_log.txt");
	lcout << "========================================================\n";
	lcout << "Erandi model Working with\n";
	lcout << "X=[" << xmin << " , " << xmax << " ]\n";
	lcout << fun.getSystemSize() << " intervals \n";
	lcout << "M: " << M << ", tau: " << tau << ", d: " << dose << "\n";
	lcout << "Parameters:\n ";
	lcout << "(r,k,aT,aTBeta,eTBeta,muD,ra, thetaD,re, muCa,thetaA,muCi,rTBeta, muBeta)=(";
	for (Doub i : PARAMETERS)
		lcout << i << ",";
	lcout << ")\n"<< S.size()<<" Points Founded:\n";
	lcout << std::fixed << std::setprecision(6);
	for (FixPoint i : S)
	{
		lcout << i.solution <<"\n-----------------------\n";
		if (i.stability)
			lcout << "Stable\n";
		else
			lcout << "Unstable\n";
		//lcout << i.solution[0] << "\t" << i.solution[1] << "\t" << i.solution[2] << "\t" << i.stability << "\n";
	}
	/*std::cin.get();*/
}

void runAl21NoDelay()
{
	Doub M = 10;
	Doub xmin = 1;// 0 zero give division over zero..
	Doub xmax = 1e9;
	Doub r1, r2, r3, r4, r5, r6, r7, r8;
	r1 = 0.0002823582562;
	r2 = 0.005;
	r3 = 2.507341315e-8;
	r4 = 0.00001129433025;
	r5 = 0.005435396431;
	r6 = 0.07863677434;
	r7 = 1.976507793;
	r8 = 0.5761520217;

	std::vector<Doub> PARAMETERS = { r1, r2, r3, r4, r5, r6, r7, r8 };
	Doub tau = 48, dose = 1e7;
	//std::cin >> tau >> d;
	ErandiSysNoDelay fun(PARAMETERS);
	std::vector<FixPoint> S;
	std::vector<double> v;

	S = al21(xmin, xmax, M, fun, tau, dose);
	//std::ofstream points("points.txt", std::ios::out | std::ios::trunc);

	LogAndStdout lcout("points_log.txt");
	lcout << "========================================================\n";
	lcout << "Erandi model No delay Working with\n";
	lcout << "X=[" << xmin << " , " << xmax << " ]\n";
	lcout << fun.getSystemSize() << " intervals \n";
	lcout << "M: " << M << ", tau: " << tau << ", d: " << dose << "\n";
	lcout << "Parameters:\n ";
	lcout << "(r1, r2, r3, r4, r5, r6, r7, r8)=(";
	for (Doub i : PARAMETERS)
		lcout << i << ",";
	lcout << ")\n" << S.size() << " Points Founded:\n";
	lcout << std::fixed << std::setprecision(6);
	for (FixPoint i : S)
	{
		lcout << i.solution << "\n-----------------------\n";
		if (i.stability)
			lcout << "Stable\n";
		else
			lcout << "Unstable\n";
		//lcout << i.solution[0] << "\t" << i.solution[1] << "\t" << i.solution[2] << "\t" << i.stability << "\n";
	}
	/*std::cin.get();*/
}

#endif
