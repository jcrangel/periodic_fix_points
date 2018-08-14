/*
Procedures to obtain the jacobian following the procedure in eq (6)&(7) of the
//paper Wei 2014
*/

#ifndef NEWTONPOINCARE_H
#define NEWTONPOINCARE_H

#include "util.h"
#include "integrationode.h"

//Gets the jacobian following the procedure in eq (6) & (7) of the
//paper Wei 2014
template <class T>
FixPoint newtonPoincare(T &fun, StateType initialCondition, Doub tau, Doub d) {
	int n = 0;            //%initialize iteration counter
	Doub eps = 1;          //%initialize error
	StateType xp{ 0, tau }; //% define the span of the computational domain
	Doub tol = 1e-6;
	int maxStep = 100;
	int N = fun.getSystemSize();
	int controlIndex = fun.getControlIndex();
	VectorEigen Xk = toEigenVector(initialCondition);
	VectorEigen Pxk;
	MatrixEigen df;
	MatrixEigen A;
	MatrixEigen I = MatrixEigen::Identity(N, N);
	VectorEigen y;
	VectorEigen Xplus;

	//StateType initCond;

	while (n < maxStep) {
		df = DFode_aprox(fun, toStateType(Xk), tau, d);  //Jacobian
		initialCondition = toStateType(Xk);
		initialCondition[controlIndex] = initialCondition[controlIndex] + d;

		//  STIFF INTEGRATION simple------------------
		StateType u(N);
		integrateSystem(fun, initialCondition, u, xp[0], xp[1]);
		Pxk = toEigenVector(u);
		//  STIFF INTEGRATION END ------------------

		A = df - I;
		//Solve Ay = (Pxk-Xk)
		y = A.colPivHouseholderQr().solve(Pxk - Xk);
		Xplus = Xk - y;

		eps = y.norm();
		if (eps < tol)
			break;

		Xk = Xplus;  //%update x
		n++;
	}
	// 0 means unstable
	Doub stability = 0;
	bool convergence = false;
	if (n < maxStep) {
		convergence = true;

		if (isStable(df))
			stability = 1;
	}
	FixPoint ss(convergence, stability, Xk);
	return ss;
}

//Calculate the jacobian of the poincare map given by equations (6) & (7)
// N-Dim system
template <class T>
MatrixEigen DFode_aprox(T &fun, StateType initialCondition,
	Doub tau, Doub d)
{
	const Doub EPS = 1.0e-8;
	//int n = initialCondition.size();
	int n = fun.getSystemSize();
	MatrixEigen df(n, n);
	StateType xh = initialCondition;
	VectorEigen fvec = evalFunInLast(fun, initialCondition, tau, d);

	for (int j = 0; j < n; j++)
	{
		Doub temp = xh[j];
		Doub h = EPS * abs(temp);
		if (h == 0.0) h = EPS;
		xh[j] = temp + h;
		h = xh[j] - temp;
		//StateType f=func(xh);
		VectorEigen f = evalFunInLast(fun, xh, tau, d);
		xh[j] = temp;
		for (int i = 0; i < n; i++)// TODO: This is ok? check if this is good in N-Dim
			df(i, j) = (f[i] - fvec[i]) / h;
	}
	return df;
};

#endif