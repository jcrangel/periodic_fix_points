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

/**********************************************************************************************//**
 * @fn	template <class T> std::vector<fixPoint> al21(double inputDomainMin, double inputDomainMax, double divisions, T &functionName, double tau, double d)
 *
 * @brief	Generic version of Algorithm 2.1 of Wei2014. Searching on N-Dimensional space given
 * 			by the size the system functionName. Every dimension is partioned into the same number
 * 			of division
 *
 * @author	Iron
 * @date	7/17/2018
 *
 * @tparam	T	Generic type parameter.
 * @param 		  	inputDomainMin	The input domain minimum.
 * @param 		  	inputDomainMax	The input domain maximum.
 * @param 		  	divisions	  	The number of divisions, the same for all dimensions .
 * @param [in,out]	functionName  	Name of the function that is a class with a functor.
 * @param 		  	tau			  	Param tau.
 * @param 		  	d			  	A double to process.
 *
 * @return	A std::vector&lt;fixPoint&gt;
 *
 * ### tparam	T	Generic type parameter for resieving any class represting a ODE system
 * 					    with a functor.
 * ### param 		 	stateMin	.
 * ### param 		 	stateMax	Vector that contains the Xmax,Ymax,Zmax...
 **************************************************************************************************/

template <class T>
std::vector<fixPoint> al21(double inputDomainMin, double inputDomainMax, double divisions, T &functionName, double tau, double d)
{
	std::vector<fixPoint> S;
	std::vector<double> A; // The set of points from which will be create the points space (cartesian product). 
	double step = (inputDomainMax - inputDomainMin) / divisions;

	for (double i = inputDomainMin; i < inputDomainMax; i += step) {
		A.push_back(i);
	}
	int N = functionName.getSystemSize();
	std::vector< std::vector<double> > pointSpace;
	if (N > 1) {
		pointSpace = cartesianProduct(A, A);// Creates A x A
		if (N > 2)
		{
			for (int i = 1; i < N; i++) // Creates A x A x A ...
				pointSpace = cartesianProduct(pointSpace, A);
		}
	}

	for (std::vector<double> i : pointSpace) {
		if (DEBUG0) {
			printVector(i);
			std::cout << "Plus step:" << step << std::endl;
		}
		std::vector<double> iplusStep(i);
		std::for_each(iplusStep.begin(), iplusStep.end(), [step](double& d) {d += step; }); //add the step to each element
		al22(i, iplusStep, S, functionName, tau, d, 1);

		if (DEBUG0) std::cout << "END searching on subspace" << std::endl;
	}

	return S;
};


void shiftRowleft(std::vector<double> &xi, std::vector<double> &xm, std::vector<double> &xf, int row) {
	double temp1;
	temp1 = xi[row];
	xi[row] = xm[row];
	xm[row] = xf[row];
	xf[row] = temp1;

};

void shiftRowRight(std::vector<double> &xi, std::vector<double> &xm, std::vector<double> &xf, int row) {
	double temp1;
	temp1 = xf[row];
	xf[row] = xm[row];
	xm[row] = xi[row];
	xi[row] = temp1;

};


//For saving a point on the space [x,y,z] 

/**********************************************************************************************//**
 * @struct	pointxyz
 *
 * @brief	A pointxyz.
 *
 * @author	Iron
 * @date	7/17/2018
 **************************************************************************************************/

struct pointxyz
{
	/** @brief	The x coordinate */
	double x;
	/** @brief	The y coordinate */
	double y;
	/** @brief	The z coordinate */
	double z;
	pointxyz(double x_, double y_, double z_) :

		/**********************************************************************************************//**
		 * @fn	x(x_), y(y_), z(z_)
		 *
		 * @brief	Constructor
		 *
		 * @author	Iron
		 * @date	7/17/2018
		 *
		 * @param	parameter1	The first parameter.
		 **************************************************************************************************/

		x(x_), y(y_), z(z_) {}


};

/**********************************************************************************************//**
 * @fn	template <class T> void al22(double xi, double xf, double yi, double yf, double zi, double zf, std::vector<fixPoint> &S, T &functionName, double tau, double d, int deepness)
 *
 * @brief	Al 22
 *
 * @author	Iron
 * @date	7/17/2018
 *
 * @tparam	T	Generic type parameter.
 * @param 		  	xi				The xi.
 * @param 		  	xf				The xf.
 * @param 		  	yi				The yi.
 * @param 		  	yf				The yf.
 * @param 		  	zi				The zi.
 * @param 		  	zf				The zf.
 * @param [in,out]	S				A std::vector&lt;fixPoint&gt; to process.
 * @param [in,out]	functionName	Name of the function.
 * @param 		  	tau				The tau.
 * @param 		  	d				A double to process.
 * @param 		  	deepness		The deepness.
 **************************************************************************************************/

template <class T>
void al22(double xi, double xf, double yi, double yf, double zi, double zf,
	std::vector<fixPoint> &S, T &functionName, double tau, double d, int deepness)
{
	if (DEBUG1)
		std::cout << deepness << std::endl;

	if (deepness <= 0)
		return;

	//MatrixXd Fx = MatrixXd::Zero(8,3); // #_#
	std::vector<Vector3d> Fx;  //Vector of solutions , each solution is a Vector3d
//First solution
	Fx.push_back(evalFunInLast(functionName, std::vector<double> {xi, yi, zi}, tau, d));

	pointxyz first(Fx[0][0] - xi, Fx[0][1] - yi, Fx[0][2] - zi);
	bool xgood = false;
	bool ygood = false;
	bool zgood = false;
	int n = 0;
	//pointxyz fminus;
	double mulx;
	double muly;
	double mulz;
	double xm = (xf + xi) / 2;
	double ym = (yf + yi) / 2;
	double zm = (zf + zi) / 2;
	//%this is wrong we need to evalute this on the middle points see algorithm
	//%2.2 step1
	for (double i : std::vector<double>{ xi,xm,xf })
		for (double j : std::vector<double>{ yi,xm,yf })
			for (double k : std::vector<double>{ zi,xm,zf })
			{
				//The first element is already computed
				if (n == 0)
				{
					n = n + 1;
					continue;
				}

				if (DEBUG1) {
					std::cout << "Step 1" << std::endl;
				}

				//Step1
				Fx.push_back(evalFunInLast(functionName, std::vector<double> {i, j, k}, tau, d));
				//Not effice to create every loop?   //#_#
				pointxyz fminus(Fx[n][0] - i, Fx[n][1] - j, Fx[n][2] - k);

				if (DEBUG1) {
					std::cout << "Step 2" << std::endl;
				}
				//            %Step 2 Check if first.x  differs at sign with fminus.x
				//            %the evaluation. Compares the first one with everything els
				mulx = fminus.x*first.x;
				muly = fminus.y*first.y;
				mulz = fminus.z*first.z;

				if (mulx < 0 || equal(mulx, 0)) xgood = true;
				if (muly < 0 || equal(muly, 0)) ygood = true;
				if (mulz < 0 || equal(mulz, 0)) zgood = true;


				//%if actually changes from 0, we take it as a sign change.
				if (equal(first.x, 0) && !equal(fminus.x, 0)) xgood = true;
				if (equal(first.y, 0) && !equal(fminus.y, 0)) ygood = true;
				if (equal(first.z, 0) && !equal(fminus.z, 0)) zgood = true;

				n = n + 1;
			}

	//If every element of Fx(xi,yi,zi)-xi has the same sign returns.(The same for %y,z)
	if (!(xgood && ygood && zgood))
		return;

	if (DEBUG1) {
		std::cout << "Step 3" << std::endl;
	}
	bool signChange = false;

	Matrix3d I3 = Matrix3d::Identity(3, 3);
	//Here we could use the F(xi,yi,zi) calculated on step2
	double firstDet = (DFode_aprox(functionName, stateType{ xi,yi,zi }, tau, d) - I3).determinant();
	n = 0;
	double Det, mult;
	for (double i : std::vector<double>{ xi,xm,xf }) {
		for (double j : std::vector<double>{ yi,xm,yf }) {
			for (double k : std::vector<double>{ zi,xm,zf })
			{
				//The first element is already computed
				if (n == 0)
				{
					n = n + 1;
					continue;
				}
				//Det = (DF(toStdVectorD(Fx[n]), tau) - I3).determinant();
				Det = (DFode_aprox(functionName, stateType{ i,j,k }, tau, d) - I3).determinant();
				//Step 4
				if (DEBUG1)
					std::cout << "step 4" << std::endl;

				mult = Det * firstDet;
				if (mult <= 0) {
					// If true, then it was a change in the determinant sing
					signChange = true;
					//Put a break that get out of here to step 5
					break;
				}
				n = n + 1;


			}
			if (mult <= 0)
				break;
		}
		if (mult <= 0)
			break;
	}




	if (signChange) {
		if (DEBUG1)
			std::cout << "step 5" << std::endl;
		//% cut at half

		std::vector<double> X = { xi, xm, xf };
		std::vector<double> Y = { yi, ym, yf };
		std::vector<double> Z = { zi, zm, zf };

		std::vector<fixPoint> temp;
		for (int i = 0; i < STATE_SIZE - 1; ++i)
			for (int j = 0; j < STATE_SIZE - 1; ++j)
				for (int k = 0; k < STATE_SIZE - 1; ++k) {
					al22(X[i], X[i + 1], Y[j], Y[j + 1], Z[k], Z[k + 1], S,
						functionName, tau, d, deepness - 1);
					//S.insert(S.end(),temp.begin(),temp.end());
				}
	}



	//%else
	//%newton
	//%disp('step 6')

	for (double i : std::vector<double>{ xi, xm, xf })
		for (double j : std::vector<double>{ yi, ym, yf })
			for (double k : std::vector<double>{ zi, zm, zf }) {

				fixPoint fixPt = newtonPoincare(functionName, stateType{ i,j,k }, tau, d);
				if (fixPt.convergent) {
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

/**********************************************************************************************//**
 * @fn	template <class T> void al22(std::vector<double> intervals, std::vector<fixPoint> &S, T &functionName, double tau, double d, int deepness)
 *
 * @brief	Algorithm 2.2 of wei 2014, Generic version for N dimensional systems
 *
 * @author	Iron
 * @date	7/17/2018
 *
 * @tparam	T	Generic type parameter.
 * @param 		  	intervals   	The intervals.
 * @param [in,out]	S				A std::vector&lt;fixPoint&gt; to process.
 * @param [in,out]	functionName	Name of the function.
 * @param 		  	tau				The tau.
 * @param 		  	d				A double to process.
 * @param 		  	deepness		The deepness.
 *
 * ### tparam	T	Generic type parameter.
 **************************************************************************************************/

template <class T>
void al22(stateType xi, stateType xf,
	std::vector<fixPoint> &S, T &functionName, double tau, double d, int deepness)
{
	if (DEBUG1)
		std::cout << deepness << std::endl;

	if (deepness <= 0)
		return;
	int N = functionName.getSystemSize();

	//MatrixXd Fx = MatrixXd::Zero(8,3); // #_#
	std::vector<VectorXd> fx;  //Vector of solutions , each solution is a eigen Vector of size N 
							   //First solution
	fx.push_back(evalFunInLast(functionName, xi, tau, d));

	stateType first;
	for (int i = 0; i < N; i++) {
		//	(Fx[0][0] - xi, Fx[0][1] - yi, Fx[0][2] - zi);
		first.push_back(fx[0][i] - xi[i]);
	}

	std::vector<bool> xgood(N, false);
	int n = 0;
	//pointxyz fminus;

	std::vector<double> mulx(N, 0);
	//double mulx;
	//double muly;
	//double mulz;

	std::vector<double> xm;
	for (int i = 0; i < N; i++) {
		xm.push_back( (xf[i] + xi[i]) / 2) ;
	}
	//double xm = (xf + xi) / 2;
	//double ym = (yf + yi) / 2;
	//double zm = (zf + zi) / 2;

	// Cartesian product of { xi,xm,xf } x { xi,xm,xf } x ... N times
	std::vector< std::vector <double> > pointSpace;
	if (N > 1) {
		// Create A x A
		pointSpace = cartesianProduct(std::vector<double>{xi[0], xm[0], xf[0]},
			std::vector<double>{xi[1], xm[1], xf[1]});
		if (N > 2)
		{
			for (int i = 2; i < N; i++) // Create A x A x A ...
				pointSpace = cartesianProduct(pointSpace, std::vector<double>{xi[i], xm[i], xf[i]});
		}
	}

	for (stateType xi : pointSpace) { // Iterates over the all the points is size n [x1,x2,...,xn]
				//The first element is already computed
		if (n == 0)
		{
			n = n + 1;
			continue;
		}

		if (DEBUG1) {
			std::cout << "Step 1" << std::endl;
		}

		//Step1				
		fx.push_back(evalFunInLast(functionName, xi, tau, d));
		//pointxyz fminus(Fx[n][0] - i, Fx[n][1] - j, Fx[n][2] - k);
		stateType fminus;
		for (int i = 0; i < N; i++) {
			//	(Fx[0][0] - xi, Fx[0][1] - yi, Fx[0][2] - zi);
			fminus.push_back(fx[n][i] - xi[i]);
		}

		if (DEBUG1) {
			std::cout << "Step 2" << std::endl;
		}
		//            %Step 2 Check if first.x  differs at sign with fminus.x
		//            %the evaluation. Compares the first one with everything els
		//mulx = fminus.x*first.x;
		//muly = fminus.y*first.y;
		//mulz = fminus.z*first.z;

		for (int i = 0; i < N; i++) {
			mulx[i] = fminus[i] * first[i];
		}

		//if (mulx < 0 || equal(mulx, 0)) xgood = true;
		//if (muly < 0 || equal(muly, 0)) ygood = true;
		//if (mulz < 0 || equal(mulz, 0)) zgood = true;
		for (int i = 0; i < N; i++) {
			if (mulx[i] < 0 || equal(mulx[i], 0)) xgood[i] = true;
		}

		//%if actually changes from 0, we take it as a sign change.
		//if (equal(first.x, 0) && !equal(fminus.x, 0)) xgood = true;
		//if (equal(first.y, 0) && !equal(fminus.y, 0)) ygood = true;
		//if (equal(first.z, 0) && !equal(fminus.z, 0)) zgood = true;
		for (int i = 0; i < N; i++) {
			if (equal(first[i], 0) && !equal(fminus[i], 0)) xgood[i] = true;
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

	if (DEBUG1) {
		std::cout << "Step 3" << std::endl;
	}
	bool signChange = false;

	MatrixXd I = MatrixXd::Identity(N, N);
	//Here we could use the F(xi,yi,zi) calculated on step2
	double firstDet = (DFode_aprox(functionName, xi, tau, d) - I).determinant();
	n = 0;
	double Det, mult;

	for (stateType xi : pointSpace) {

		//The first element is already computed
		if (n == 0)
		{
			n = n + 1;
			continue;
		}
		//Det = (DF(toStdVectorD(Fx[n]), tau) - I).determinant();
		Det = (DFode_aprox(functionName, xi, tau, d) - I).determinant();
		//Step 4
		if (DEBUG1)
			std::cout << "step 4" << std::endl;

		mult = Det * firstDet;
		if (mult <= 0) {
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
	xiTemp will be
	*/
	if (signChange) {
		if (DEBUG1)
			std::cout << "step 5" << std::endl;

		std::vector<double> xiTemp;
		std::vector<double> xmTemp;
		std::vector<double> xfTemp;
		std::vector<int>  bitArray(N, 0);
		xiTemp = xi;
		xmTemp = xm;
		xfTemp = xf;

		for (int i = 0; i < pow(2, N) - 1; i++) {
			al22(xiTemp,xmTemp,S,functionName, tau, d, deepness - 1);
			int j = 0;
			while (bitArray[j] == 1) {
				bitArray[j] = 0;
				shiftRowRight(xiTemp, xmTemp, xfTemp, j);
				j++;
			}
			bitArray[j] = 1;
			shiftRowleft(xiTemp, xmTemp, xfTemp, j);
		}
		al22(xiTemp, xmTemp, S, functionName, tau, d, deepness - 1);

	}



	//%else
	//%newton
	//%disp('step 6')

	for (stateType xi : pointSpace)
	{

		fixPoint fixPt = newtonPoincare(functionName, xi, tau, d);
		if (fixPt.convergent) {
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










