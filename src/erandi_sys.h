/**********************************************************************************************/ /**
 * Copyright 2018 Julio C. Rangel
 * @file	erandi_sys.h
 *
 * @brief. 
 * The Erandi's ODE system with no delay 
 * 
 **************************************************************************************************/

#ifndef ERANDI_SYS_H
#define ERANDI_SYS_H

#include "util.h"
#include <cmath>
class ErandiSys
{
private:
	std::vector<Doub> parameters;
	const int systemSize = 6; // number of states
	int controlIndex = 3;	 // The index of the state where the control happens
	const int numParameters = 14;
    Doub r,k,aT,aTBeta,eTBeta,muD,ra;
    Doub thetaD,re, muCa,thetaA,muCi,rTBeta;
    Doub muBeta;

	void loadParameters()
	{
		r = parameters[0];
		k = parameters[1];
		aT = parameters[2];
		aTBeta = parameters[3];
		eTBeta = parameters[4];
		muD = parameters[5];
		ra = parameters[6];
		thetaD = parameters[7];
		re = parameters[8];
		muCa = parameters[9];
		thetaA = parameters[10];
		muCi = parameters[11];
		rTBeta = parameters[12];
		muBeta = parameters[13];
	}

	void changeParameters(std::vector<Doub> newParameters)
	{
		parameters = newParameters;
		loadParameters();
	}
public:
	ErandiSys(std::vector<Doub> parameters) : parameters(parameters)
	{
		loadParameters();
	}

	void operator()(const StateType &x, StateType &dxdt, Doub /*t*/)
	{
		Doub T = x[0];
		Doub Ca = x[1];
		Doub Ci = x[2];
		Doub Den = x[3];
		Doub D = x[4];
		Doub F = x[5];

		dxdt[0] = r * T * std::log(k / T) - (aT * T * Ca) * (( (1-aTBeta)*eTBeta / (eTBeta + F) ) + aTBeta);
		dxdt[1] = ra*Ci*(Den /(Den +thetaD)) + re*Den*(Ca /(Ca + thetaA)) - muCa*Ca; 
		dxdt[2] = -ra*Ci*(Den /(Den +thetaD)) - muCi*Ci; 
		dxdt[3] = -muD * Den;
		dxdt[4] = aTBeta * T - muBeta * F;
	}

	int getSystemSize() {
		return systemSize;
	}
	int getControlIndex()
	{
		return controlIndex;
	}
};


class ErandiSysNoDelay
{
private:
	std::vector<Doub> parameters;
	const int systemSize = 4; // number of states
	int controlIndex = 3;	 // The index of the state where the control happens
	const int numParameters = 8;
	Doub r1, r2, r3, r4, r5, r6, r7, r8;

	void loadParameters()
	{
		r1 = parameters[0];
		r2 = parameters[1];
		r3 = parameters[2];
		r4 = parameters[3];
		r5 = parameters[4];
		r6 = parameters[5];
		r7 = parameters[6];
		r8 = parameters[7];
	}

	void changeParameters(std::vector<Doub> newParameters)
	{
		parameters = newParameters;
		loadParameters();
	}
public:
	ErandiSysNoDelay(std::vector<Doub> parameters) : parameters(parameters)
	{
		loadParameters();
	}

	void operator()(const StateType &x, StateType &dxdt, Doub /*t*/)
	{
		Doub T = x[0];
		Doub C = x[1];
		Doub F = x[2];
		Doub D = x[3];

		dxdt[0] = r1 * T * (1 -r2*T) - r3 * C * (T /(T + 1)) * ((1 + 0.69*F) / (1 + F));
		dxdt[1] = C*(D / (D + 1)) + r4 - r5*C;
		dxdt[2] = r6 * T - r7 * F;
		dxdt[3] = -r8 * D;
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
