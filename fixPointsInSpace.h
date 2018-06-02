/*
In this file are the methods & classes for finding fixed points on space
That is the algorithms 2.1 and 2.2 of Wei 2014
*/

#ifndef FIXPOINTSINSPACE_H
#define FIXPOINTSINSPACE_H

#include "all.h"


//Return a set(vector) of fixed points
template <class T>
std::vector<fixPoint> al21(double xmax, double ymax, double zmax, double M, double N, double L, T &functionName, double tau, double d)
{
	std::vector<fixPoint> S;
	double xstep = xmax / M;
	double ystep = ymax / N;
	double zstep = zmax / L;
	double xmin = 0;
	double ymin = 0;
	double zmin = 0;

	std::vector<fixPoint> temp;

//concat
//vector1.insert( vector1.end(), vector2.begin(), vector2.end() );

	double i, j, k;
	for (i = xmin; i <= xmax; i += xstep)
		for (j = ymin; j <= ymax; j += ystep)
			for ( k = zmin; k <= zmax; k += zstep)
			{
				//Since where pasing S as a reference, each al22 modify S. If we want parallelism S must be
				//part of a reduction or something
				al22(i, i + xstep, j, j + ystep, k, k + zstep, S,functionName, tau, d, 6) ;
                //S.insert(S.end(),temp.begin(),temp.end());
			}
	return S;
}

//%Evaluate F(x,y,z)
template <class T>
Vector3d evalFunInLast(T &functionName, stateType initialCondition, double tau, double d)
{
	vectorBoost x(3);
	std::copy(initialCondition.begin(), initialCondition.end(), x.begin());
	itikBanks_jacobi_stiff J(PARAMETERS);
	itikBanks_stiff fun(PARAMETERS);
	stateType xp = { 0, tau };
	initialCondition[CONTROL_POS] = initialCondition[CONTROL_POS] + d;

	size_t num_of_steps = integrate_const(make_dense_output< rosenbrock4< double > >(1.0e-6, 1.0e-6),
		std::make_pair(fun, J),
		x, 0.0, 50.0, 0.01);

		std::vector<double> res;
		std::copy(x.begin(), x.end(), res.begin());

		Vector3d v(res.data());
		return v;
	//initialCondition has the last
}

struct pointxyz
{
	double x;
	double y;
	double z;
	pointxyz(double x_, double y_, double z_) :
		x(x_), y(y_), z(z_) {}
};


template <class T>
void al22(double xi, double xf, double yi, double yf, double zi, double zf,
	std::vector<fixPoint> &S, T &functionName, double tau, double d, int deepness)
{

	if (deepness <= 0)
		return ;

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
	//%this is wrong we need to evalute this on the middle points see algorithm
	//%2.2 step1
	for (double i : std::vector<double>{ xi,xf })
		for (double j : std::vector<double>{ yi,yf })
			for (double k : std::vector<double>{ zi,zf })
			{
				//The first element is already computed
				if (n == 0)
				{
					n = n + 1;
					continue;
				}

				if (DEBUG) {
					std::cout << "Step 1" << std::endl;
				}

				//Step1
				Fx.push_back(evalFunInLast(functionName, std::vector<double> {xi, yi, zi}, tau, d));
				//Not effice to create every loop?   //#_#
				pointxyz fminus(Fx[n][0] - xi, Fx[n][1] - yi, Fx[n][2] - zi);

				if (DEBUG) {
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
		return ;

	if (DEBUG) {
		std::cout << "Step 3" << std::endl;
	}
	bool signChange = false;

	Matrix3d I3 = Matrix3d::Identity(3, 3);
	double firstDet = (DFitikBanks(toStdVectorD(Fx[0]), tau) - I3).determinant();
	n = 0;
	double Det, mult;
	for (double i : std::vector<double>{ xi,xf }) {
		for (double j : std::vector<double>{ yi,yf }) {
			for (double k : std::vector<double>{ zi,zf })
			{
				//The first element is already computed
				if (n == 0)
				{
					n = n + 1;
					continue;
				}
				Det = (DFitikBanks(toStdVectorD(Fx[n]), tau) - I3).determinant();
				//Step 4
				if (DEBUG)
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

	double xm = (xf + xi) / 2;
	double ym = (yf + yi) / 2;
	double zm = (zf + zi) / 2;


	if (signChange) {
		//%Do step 5;
		//%disp('step 5')
		//% cut at half

		std::vector<double> X = { xi, xm, xf };
		std::vector<double> Y = { yi, ym, yf };
		std::vector<double> Z = { zi, zm, zf };

        std::vector<fixPoint> temp;
		for (int i = 0; i < STATE_SIZE; ++i)
			for (int j = 0; j < STATE_SIZE; ++j)
				for (int k = 0; k < STATE_SIZE; ++k) {
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
					S.push_back(fixPt);
					goto END; //Yes, a goto
				}
			}
END:

	return ;
}

#endif










