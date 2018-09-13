/*
In this file are the methods & classes for finding fixed points on space
That is the algorithms 2.1 and 2.2 of Wei 2014
*/

#ifndef PARAMETERSWEPT_H
///< .
#define PARAMETERSWEPT_H

#include "util.h"
#include "integrationode.h"
#include "newtonpoincare.h"
#include "fixpointsinspace.h"
#include <algorithm>
#include <map>

/*
The parameter swepting over tau and d. For each (tau,d) 
check if the periodic fix point founded is under and upper bound. 
In this way we check if the therapy controls the tumor under certain threhold.
For the tumor variable, we need to provide the index of the were is the bounded variable in func.
*/
template <class T>
setCriticalCells al23_tau_d(Doub tau0, Doub taum, Doub d0, Doub dm,  int M, int N, 
	T &func, Doub upperBound, StateType x0, int indexBoundedVar)
{
	setCriticalCells Cri;
	Doub htau = (taum - tau0) / M;
	Doub hd = (dm - d0) / N;
	//pt = [];
	//%A map : (tau, d)->{Set of fixed points}
	//std::map<std::string, FixPoint*> S;

	//%A map : (tau, d)-># 1 if the periodic fix point under the threshold, 2 if not,
	//0 periodic point 
	// false otherwise
	std::map<std::string, bool> S;
	
	//%A map : (tau, d)->#Fix point founded
	std::map<std::string, bool> Q;
	
	StateType periodicPoint;
	//%Step 2 , we're swepting over r2 and d. In the corner points
	for (int x = 0; x <= M; x++) { //tau = r2
		for (int y = 0; y <= N; y++) { //% d = d
			StateType ic;
			ic = x0;
			Doub i = tau0 + x * htau;
			Doub j = d0 + y * hd;
			//Chanching d using  j and we let the system to run withput therapy for one period
			//of size i
			Doub max = findperiodicPointGetMax(0, i, ic, func, i, j,20);
			//TODO: Save also when the algorithm doesn't found the periodic point
			std::string key = strtokey(i, j);

			if (max == -1) {
				Q.insert(std::pair<std::string, bool>(key, false));
				continue; // No point founded
			}else
				Q.insert(std::pair<std::string, bool>(key, true));
			//Save the first three colums of O.That is all the fixed points in (i, j).
			//zaiss = size(O);


			S.insert(std::pair<std::string, bool>(key, max <= upperBound));

			//The numbers of rows is the numbers fixed points.
		}
	}
	//Step 3
	for (int x = 0; x < M; x++) { 
		for (int y = 0; y < N; y++) { 
			Doub i = tau0 + (x+0.5) * htau;
			Doub j = d0 + (y+0.5) * hd;
			StateType ic;
			ic = x0;
			//Chanching d using  j and we let the system to run withput therapy for one period
			//of size i
			Doub max = findperiodicPointGetMax(0, i, ic, func, i, j, 20);
			std::string key = strtokey(i, j);

			if (max == -1) {
				Q.insert(std::pair<std::string, bool>(key, false));
				continue; // No point founded
			}
			else
				Q.insert(std::pair<std::string, bool>(key, true));
			//Save the first three colums of O.That is all the fixed points in (i, j).
			//zaiss = size(O);

			S.insert(std::pair<std::string, bool>(key, max <= upperBound));
			//The numbers of rows is the numbers fixed points.
		}
	}

//Step 4
#pragma omp declare reduction (merge : setCriticalCells : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp parallel for reduction(merge: Cri)
	for (int x = 0; x < M; x++) {
		for (int y = 0; y < N; y++) {

			Doub i = tau0 + x * htau;
			Doub iplus = tau0 + (x + 1)*htau;
			Doub j = d0 + y * hd;
			Doub jplus = d0 + (y + 1)*hd;

			Doub ihalf = tau0 + (x + 0.5)*htau;
			Doub jhalf = d0 + (y + 0.5)*hd;

			//if some point have a diferent value that cell is critical
			if (!(S[strtokey(i, j)] == S[strtokey(iplus, j)] == S[strtokey(i, jplus)] ==
				S[strtokey(iplus, jplus)] == S[strtokey(ihalf, jhalf)] )) {
				//We normalize the tau with the maximum tau
				//Cri.push_back(cell2D(point(i/taum, j), point(iplus / taum, j), point(i / taum, jplus), point(iplus / taum, jplus)));
				
				//Make this call thread safe. All should share variable Cri
				al24_tau_d(i, iplus, j, jplus, func, upperBound, x0, indexBoundedVar,Cri,2);

			}

			// call al24

		}
	}

	return Cri;
};


template <class T>
void al24_tau_d(Doub tau0, Doub taum, Doub d0, Doub dm,
	T &func, Doub upperBound, StateType x0, int indexBoundedVar, setCriticalCells &Cri,int level)
{
	int M = 2;
	int N = 2;
	Doub htau = (taum - tau0) / M;
	Doub hd = (dm - d0) / N;
	int maxLevel = 5;
	//pt = [];
	//%A map : (tau, d)->{Set of fixed points}
	//std::map<std::string, FixPoint*> S;

	//%A map : (tau, d)-># true if the periodic fix point under the threshold
	// false otherwise
	std::map<std::string, bool> S;

	//%A map : (tau, d)->#Number stable fix points
	std::map<std::string, bool> Q;

	StateType periodicPoint;
	//%Step 2 , we're swepting over r2 and d. In the corner points
	for (int x = 0; x <= M; x++) { //tau = r2
		for (int y = 0; y <= N; y++) { //% d = d
			StateType ic;
			ic = x0;
			Doub i = tau0 + x * htau;
			Doub j = d0 + y * hd;
			//Chanching d using  j and we let the system to run withput therapy for one period
			//of size i
			Doub max = findperiodicPointGetMax(0, i, ic, func, i, j, 20);
			//TODO: Save also when the algorithm doesn't found the periodic point
			std::string key = strtokey(i, j);

			if (max == -1) {
				Q.insert(std::pair<std::string, bool>(key, false));
				continue; // No point founded
			}
			else
				Q.insert(std::pair<std::string, bool>(key, true));
			//Save the first three colums of O.That is all the fixed points in (i, j).
			//zaiss = size(O);

			S.insert(std::pair<std::string, bool>(key, max <= upperBound));

			//The numbers of rows is the numbers fixed points.
		}
	}
	//Step 3
	for (int x = 0; x < M; x++) {
		for (int y = 0; y < N; y++) {
			Doub i = tau0 + (x + 0.5) * htau;
			Doub j = d0 + (y + 0.5) * hd;
			StateType ic;
			ic = x0;
			//Chanching d using  j and we let the system to run withput therapy for one period
			//of size i
			Doub max = findperiodicPointGetMax(0, i, ic, func, i, j, 20);
			std::string key = strtokey(i, j);

			if (max == -1) {
				Q.insert(std::pair<std::string, bool>(key, false));
				continue; // No point founded
			}
			else
				Q.insert(std::pair<std::string, bool>(key, true));
			//Save the first three colums of O.That is all the fixed points in (i, j).
			//zaiss = size(O);

			S.insert(std::pair<std::string, bool>(key, max <= upperBound));
			//The numbers of rows is the numbers fixed points.
		}
	}

	//Step 4
	for (int x = 0; x < M; x++) {
		for (int y = 0; y < N; y++) {

			Doub i = tau0 + x * htau;
			Doub iplus = tau0 + (x + 1)*htau;
			Doub j = d0 + y * hd;
			Doub jplus = d0 + (y + 1)*hd;

			Doub ihalf = tau0 + (x + 0.5)*htau;
			Doub jhalf = d0 + (y + 0.5)*hd;
			//Step 5
			//if some point have a diferent value that cell is critical
			if ( !(S[strtokey(i, j)] == S[strtokey(iplus, j)] == S[strtokey(i, jplus)] ==
				S[strtokey(iplus, jplus)] == S[strtokey(ihalf, jhalf)]) || 
				!(Q[strtokey(i, j)] == Q[strtokey(iplus, j)] == Q[strtokey(i, jplus)] ==
				Q[strtokey(iplus, jplus)] == Q[strtokey(ihalf, jhalf)]) 
				){
				//We normalize the tau with the maximum tau
				if (level == maxLevel)// Insert t					
					Cri.push_back(cell2D(point(i , j), point(iplus , j), point(i , jplus), point(iplus , jplus)));
				else
					al24_tau_d(i, iplus, j, jplus, func, upperBound, x0, indexBoundedVar, Cri, level+1);
			}

			// call al24

		}
	}

};


#endif
