/*
In this file are the methods & classes for finding fixed points on space
That is the algorithms 2.1 and 2.2 of Wei 2014
*/

#ifndef FIXPOINTSINSPACE_H
///< .
#define FIXPOINTSINSPACE_H

#include "util.h"
#include "integrationode.h"
#include "newtonpoincare.h"
#include <algorithm>

/**********************************************************************************************/ /**
 * @fn	template <class T> std::vector<fixPoint> al21(Doub inputDomainMin, Doub inputDomainMax, Doub divisions, T &functionName, Doub tau, Doub d)
 *
 * @brief	Generic version of Algorithm 2.1 of Wei2014. Searching on N-Dimensional space given
 * 			by the size the system functionName. Every dimension is partioned into the same number
 * 			of division
 *
 * @author	Iron
 * @date	7/17/2018
 *
 * @tparam	T	Generic type parameter.
 * @param 		  	inputDomainMin	The input domain minimum, this is for every state x1,x2,x3,...,xn
 * @param 		  	inputDomainMax	The input domain maximum.
 * @param 		  	divisions	  	The number of divisions, the same for all states .
 * @param [in,out]	functionName  	Name of the function that is a class with a functor.
 * @param 		  	tau			  	Param tau.
 * @param 		  	d			  	A Doub to process.
 *
 * @return	A std::vector&lt;fixPoint&gt;
 *
 * ### tparam	T	Generic type parameter for resieving any class represting a ODE system
 * 					    with a functor.
 * ### param 		 	stateMin	.
 * ### param 		 	stateMax	Vector that contains the Xmax,Ymax,Zmax...
 **************************************************************************************************/

template <class T>
std::vector<FixPoint> al21(Doub inputDomainMin, Doub inputDomainMax, Doub divisions, T &functionName, Doub tau, Doub d)
{
	std::vector<FixPoint> S;
	std::vector<Doub> A; // The set of points from which will be create the points space (cartesian product).
	Doub step = (inputDomainMax - inputDomainMin) / divisions;
	//Create the first set of intervals e.g [0 0.2 0.4 0.6 0.8 1]
	// for (Doub i = inputDomainMin; i < inputDomainMax; i += step) {
	// 	A.push_back(i);
	// }
	int N = functionName.getSystemSize();
	// std::vector< std::vector<Doub> > pointSpace;
	// if (N > 1) {
	// 	pointSpace = cartesianProduct(A, A);// Creates A x A
	// 	if (N > 2)
	// 	{
	// 		for (int i = 1; i < N; i++) // Creates A x A x A ...
	// 			pointSpace = cartesianProduct(pointSpace, A);
	// 	}
	// }
	// Iterates over all subdomains
	std::vector<Doub> point(N, inputDomainMin);

	for (int i = 0; i < std::pow(divisions + 1, N); i++)
	{
		if (DEBUG0)
		{
			printVector(point);
			std::cout << "Plus step:" << step << std::endl;
		}

		std::vector<Doub> pointPlusStep(point);
		std::for_each(pointPlusStep.begin(), pointPlusStep.end(), [step](Doub &d) { d += step; }); //add the step to each element
		al22(point, pointPlusStep, S, functionName, tau, d, 1);
		nextPointSubdomain(point, inputDomainMin, inputDomainMax, step);

		if (DEBUG0)
			std::cout << "END searching on subspace" << std::endl;
	}

	return S;
};

/**********************************************************************************************/ /**
 * @fn	void shiftRowleft(std::vector<Doub> &xi, std::vector<Doub> &xm, std::vector<Doub> &xf, int row)
 *
 * @brief	Shift row to the left
 *
 * @author	Iron
 * @date	7/31/2018
 *
 * @param [in,out]	xi 	The column of xi.
 * @param [in,out]	xm 	The column of xm.
 * @param [in,out]	xf 	The column of xf.
 * @param 		  	row	The row index to be shifted.
 **************************************************************************************************/

void shiftRowleft(std::vector<Doub> &xi, std::vector<Doub> &xm, std::vector<Doub> &xf, int row)
{
	Doub temp1;
	temp1 = xi[row];
	xi[row] = xm[row];
	xm[row] = xf[row];
	xf[row] = temp1;
};

void shiftRowRight(std::vector<Doub> &xi, std::vector<Doub> &xm, std::vector<Doub> &xf, int row)
{
	Doub temp1;
	temp1 = xf[row];
	xf[row] = xm[row];
	xm[row] = xi[row];
	xi[row] = temp1;
};

/**********************************************************************************************/ /**
 * @fn	template <class T> void al22(StateType xi, StateType xf, std::vector<fixPoint> &S, T &functionName, Doub tau, Doub d, int deepness)
 *
 * @brief	Algorithm 2.2 of wei 2014, Generic version for N dimensional systems
 *
 * @author	Iron
 * @date	7/17/2018
 *
 * @tparam	T	Generic type parameter.
 * @param 		  	xi				The vector for the initial values of each state.
 * @param 		  	xf				The vector for the final values of each state.
 * @param [in,out]	S				A std::vector&lt;fixPoint&gt; to process.
 * @param [in,out]	functionName	Name of the function.
 * @param 		  	tau				The tau.
 * @param 		  	d				A Doub to process.
 * @param 		  	deepness		The deepness.
 *
 * ### tparam	T	Generic type parameter.
 * ### param 		 	intervals	The intervals.
 **************************************************************************************************/

template <class T>
void al22(std::vector<Doub> xi, std::vector<Doub> xf,
	std::vector<FixPoint> &S, T &functionName, Doub tau, Doub d, int deepness)
{
	if (DEBUG1)
		std::cout << deepness << std::endl;

	if (deepness <= 0)
		return;
	int N = functionName.getSystemSize();

	std::vector<VectorEigen> fx; //Vector of solutions , each solution is a eigen Vector of size N
								 //First solution
	fx.push_back(evalFunInLast(functionName, xi, tau, d));

	StateType first;
	for (int i = 0; i < N; i++)
	{
		//	(Fx[0][0] - xi, Fx[0][1] - yi, Fx[0][2] - zi);
		first.push_back(fx[0][i] - xi[i]);
	}

	std::vector<bool> xgood(N, false);
	int n = 0;
	//pointxyz fminus;

	std::vector<Doub> mulx(N, 0);
	//Doub mulx;
	//Doub muly;
	//Doub mulz;

	std::vector<Doub> xm;
	for (int i = 0; i < N; i++)
	{
		xm.push_back((xf[i] + xi[i]) / 2);
	}
	//Doub xm = (xf + xi) / 2;
	//Doub ym = (yf + yi) / 2;
	//Doub zm = (zf + zi) / 2;

	// Cartesian product of { xi,xm,xf } x { xi,xm,xf } x ... N times
	// std::vector<std::vector<Doub> > pointSpace;
	// if (N > 1)
	// {
	// 	// Create A x A
	// 	pointSpace = cartesianProduct(std::vector<Doub>{xi[0], xm[0], xf[0]},
	// 								  std::vector<Doub>{xi[1], xm[1], xf[1]});
	// 	if (N > 2)
	// 	{
	// 		for (int i = 2; i < N; i++) // Create A x A x A ...
	// 			pointSpace = cartesianProduct(pointSpace, std::vector<Doub>{xi[i], xm[i], xf[i]});
	// 	}
	// }

	StateType point(xi);
	//3, there's always three options
	for (int i = 1; i < std::pow(point.size(), 3); i++)
	{ // Iterates over the all the points is size n [x1,x2,...,xn]
		//The first element is already computed

		nextPointSubdomain(point, xi, xf);
		//We can actually prevent this if the modify nextPointSubdomain to work in reverse
		StateType reversePoint(point.size());
		std::reverse_copy(point.begin(), point.end(), reversePoint.begin());

		if (DEBUG1)
		{
			std::cout << n << " : Step 1" << std::endl;
		}

		//Step1									We have to use the resersed one
		fx.push_back(evalFunInLast(functionName, reversePoint, tau, d));
		//pointxyz fminus(Fx[n][0] - i, Fx[n][1] - j, Fx[n][2] - k);
		StateType fminus;
		for (int i = 0; i < N; i++)
		{
			//	(Fx[0][0] - point, Fx[0][1] - yi, Fx[0][2] - zi);
			fminus.push_back(fx[n][i] - reversePoint[i]);
		}

		if (DEBUG1)
		{
			std::cout << "Step 2" << std::endl;
		}
		//            %Step 2 Check if first.x  differs at sign with fminus.x
		//            %the evaluation. Compares the first one with everything els
		//mulx = fminus.x*first.x;
		//muly = fminus.y*first.y;
		//mulz = fminus.z*first.z;

		for (int i = 0; i < N; i++)
		{
			mulx[i] = fminus[i] * first[i];
		}

		//if (mulx < 0 || equal(mulx, 0)) xgood = true;
		//if (muly < 0 || equal(muly, 0)) ygood = true;
		//if (mulz < 0 || equal(mulz, 0)) zgood = true;
		for (int i = 0; i < N; i++)
		{
			if (mulx[i] < 0 || equal(mulx[i], 0))
				xgood[i] = true;
		}

		//%if actually changes from 0, we take it as a sign change.
		//if (equal(first.x, 0) && !equal(fminus.x, 0)) xgood = true;
		//if (equal(first.y, 0) && !equal(fminus.y, 0)) ygood = true;
		//if (equal(first.z, 0) && !equal(fminus.z, 0)) zgood = true;
		for (int i = 0; i < N; i++)
		{
			if (equal(first[i], 0) && !equal(fminus[i], 0))
				xgood[i] = true;
		}
		n = n + 1;
	}

	//If every element of Fx(xi,yi,zi)-xi has the same sign returns.(The same for %y,z)
	//if (!(xgood && ygood && zgood))  !(xgood && false && zgood) = !(false) = true, if there a false exit
	//	return;
	//
	// If there a false exit
	if (std::find(xgood.begin(), xgood.end(), false) != xgood.end())
		return;

	if (DEBUG1)
	{
		std::cout << "Step 3" << std::endl;
	}
	bool signChange = false;

	MatrixEigen I = MatrixEigen::Identity(N, N);
	//Here we could use the F(xi,yi,zi) calculated on step2
	Doub firstDet = (DFode_aprox(functionName, xi, tau, d) - I).determinant();
	n = 0;
	Doub Det, mult;

	point = xi;
	//3, there's always three options
	for (int i = 1; i < std::pow(point.size(), 3); i++)
	{ // Iterates over the all the points is size n [x1,x2,...,xn]
		//The first element is already computed
		nextPointSubdomain(point, xi, xf);
		StateType reversePoint(point.size());
		std::reverse_copy(point.begin(), point.end(), reversePoint.begin());
		//Det = (DF(toStdVectorD(Fx[n]), tau) - I).determinant();
		Det = (DFode_aprox(functionName, reversePoint, tau, d) - I).determinant();
		//Step 4
		if (DEBUG1)
			std::cout << "step 4" << std::endl;

		mult = Det * firstDet;
		if (mult <= 0)
		{
			// If true, then it was a change in the determinant sing
			signChange = true;
			//Put a break that get out of here to step 5
			break;
		}
		n = n + 1;
	}

	/*
	N - Dimentional case
	For sending to al22 all 2^N subintervals.
	For a 3D dimentional system we have the columns:
	xi		xm		xf
	0.1		0.3		0.5  x<--
	0.1		0.3		0.5  y<--
	0.1		0.3		0.5	 z<--

	and we send to al22(xi , xm ...) with xi and xm as the initial
	and final column of data. First row for the first dependent variable
	e.g x, the second for y, and so on.
	For obtain one subdomain we shift the each row to the left :
	xi		xm		xf
	0.3		0.5  	0.1		<--
	0.1		0.3		0.5
	0.1		0.3		0.5
	To send all possible subdomains we make all possible permutations
	of shifts. In this case are 2^3 = 8 in total. Each shift is controlled
	in a array of shifts <<bitArray>> 0 if there's no shift and 1 if there's
	a shift,[0,0,0,1,0,1,...] we generate all 0 & 1's permutations with
	algorithm in p. 437 Rosen- Discrete Math . We make left shift if there's a 1
	and no shift if there's a 0.

	Probably this can be done using nextPointSubdomain
	*/
	if (signChange)
	{
		if (DEBUG1)
			std::cout << "step 5" << std::endl;

		std::vector<Doub> xiTemp;
		std::vector<Doub> xmTemp;
		std::vector<Doub> xfTemp;
		std::vector<int> bitArray(N, 0);
		xiTemp = xi;
		xmTemp = xm;
		xfTemp = xf;

		for (int i = 0; i < pow(2, N) - 1; i++)
		{
			al22(xiTemp, xmTemp, S, functionName, tau, d, deepness - 1);
			int j = 0;
			while (bitArray[j] == 1)
			{
				bitArray[j] = 0;
				shiftRowRight(xiTemp, xmTemp, xfTemp, j);
				j++;
			}
			bitArray[j] = 1;
			shiftRowleft(xiTemp, xmTemp, xfTemp, j);
		}
		al22(xiTemp, xmTemp, S, functionName, tau, d, deepness - 1); //TODO: Necesary??
	}

	//%else
	//%newton
	//%disp('step 6')

	point = xi;
	//3, there's always three options
	for (int i = 1; i < std::pow(point.size(), 3); i++)
	{ // Iterates over the all the points is size n [x1,x2,...,xn]
		//The first element is already computed
		nextPointSubdomain(point, xi, xf);
		StateType reversePoint(point.size());
		std::reverse_copy(point.begin(), point.end(), reversePoint.begin());

		FixPoint fixPt = newtonPoincare(functionName, reversePoint, tau, d);
		if (fixPt.convergent)
		{
			//If have negatives dont insert it , sometime wild -0 appear
			// is always true than -0 < 0 is false?
			if (pointHaveNegatives(fixPt))
				goto END;

			if (!pointIsInSet(fixPt, S)) //
				S.push_back(fixPt);

			goto END; //Yes, a goto
		}
	}
END:
	return;
};

#endif
