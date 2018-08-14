/**********************************************************************************************/ /**
 * Copyright 2018 Julio C. Rangel
 * @file	castiglione.h
 *
 * @brief From paper. Optimal Control in a Model of Dendritic Cell Transfection Cancer Immunotherapy
 **************************************************************************************************/

#ifndef MURINEMODEL_H
#define MURINEMODEL_H

#include "util.h"
#include <cmath>
class Castiglione
{
private:
	std::vector<Doub> parameters;
	const int systemSize = 5; // number of states
	int controlIndex = 3;	 // The index of the state where the control happens
	const int numParameters = 17;
	Doub a0, b0, c0, d0, f0;
	Doub a1, b1, c1, d1, f1;
	Doub d2, f2, e2;
	Doub e3;
	Doub a4, c4, e4;

	void loadParameters()
	{
		a0 = parameters[0];
		b0 = parameters[1];
		c0 = parameters[2];
		d0 = parameters[3];
		f0 = parameters[4];
		a1 = parameters[5];
		b1 = parameters[6];
		c1 = parameters[7];
		d1 = parameters[8];
		f1 = parameters[9];
		d2 = parameters[10];
		f2 = parameters[11];
		e2 = parameters[12];
		e3 = parameters[13];
		a4 = parameters[14];
		c4 = parameters[15];
		e4 = parameters[16];
	}

	void changeParameters(std::vector<Doub> newParameters)
	{
		parameters = newParameters;
		loadParameters();
	}
public:
	Castiglione(std::vector<Doub> parameters) : parameters(parameters)
	{
		loadParameters();
	}
	template <class T>
	void operator()(const T &x, T &dxdt, Doub /*t*/)
	{
		Doub H = x[0];
		Doub C = x[1];
		Doub M = x[2];
		Doub D = x[3];
		Doub I = x[4];

		dxdt[0] = a0 - b0 * H + c0 * D*(d0 * H*(1 - H / f0));
		dxdt[1] = a1 - b1 * C + c1 * I*(M + D)*(d1 * C*(1 - C / f1));
		dxdt[2] = (d2 * M*(1 - M / f2)) - e2 * M*C;
		dxdt[3] = -e3 * D*C;
		dxdt[4] = a4 * H*D - c4 * C*I - e4 * I;
	}
	/**
 * This is needed for the Stiff integrator
 *
*/
//void operator()(const VectorBoost &x, VectorBoost &dxdt, Doub /*t*/)
//{
//}

/**
 * The jacobian of the system
*/
	void operator()(const VectorBoost &x, MatrixBoost &M, Doub /*t*/, VectorBoost &dfdt) {
	}

	int getSystemSize() {
		return systemSize;
	}
	int getControlIndex()
	{
		return controlIndex;
	}
};

#endif
