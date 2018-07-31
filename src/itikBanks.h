/*
//Itikbanks model from Wei2014
*/



#ifndef ITIKBANKS_H
#define ITIKBANKS_H


#include "util.h"

/**********************************************************************************************//**
 * @struct	itikBanks
 *
 * @brief	The class that defines the ODE system itikBanks .
 *
 * @author	Iron
 * @date	7/17/2018
 **************************************************************************************************/

class ItikBanks
{
private:
	std::vector<Doub> parameters;//TODO: Check for the correct size
	const int systemSize = 3; // number of ODEs of the system
	int controlIndex = 2; // The index of the state where the control happens

	//These are for expressing more easy the system formulas
	Doub a12, a21, a13, a31;
	Doub r2, r3, d3, k3;

	void loadParameters() {
		a12 = parameters[0];
		a13 = parameters[1];
		r2 = parameters[2];
		a21 = parameters[3];
		r3 = parameters[4];
		k3 = parameters[5];
		a31 = parameters[6];
		d3 = parameters[7];
	}

	void changeParameters(std::vector<Doub> newParameters) {
		parameters = newParameters;
		loadParameters();
	}

    /**********************************************************************************************//**
     * @fn	itikBanks(std::vector<Doub> parameters)
     *
     * @brief	Constructor
     *
     * @author	Iron
     * @date	7/17/2018
     *
     * @param	parameters	Options for controlling the operation.
     **************************************************************************************************/
public:
    ItikBanks(std::vector<Doub> parameters): parameters(parameters){
		loadParameters();
	}

    void operator()(const StateType &x, StateType &dxdt, Doub /*t*/)
    {
        dxdt[0] = x[0] * (1 - x[0]) - a12 * x[0] * x[1] - a13 * x[0] * x[2];
        dxdt[1] = r2 * x[1] * (1 - x[1]) - a21 * x[0] * x[1];
        dxdt[2] = (r3 * x[0] * x[2]) / (k3 + x[0]) - a31 * x[0] * x[2] - d3 * x[2];
    }

    void operator()(const VectorBoost &x, VectorBoost &dxdt, Doub /*t*/)
    {
        dxdt[0] = x[0] * (1 - x[0]) - a12 * x[0] * x[1] - a13 * x[0] * x[2];
        dxdt[1] = r2 * x[1] * (1 - x[1]) - a21 * x[0] * x[1];
        dxdt[2] = (r3 * x[0] * x[2]) / (k3 + x[0]) - a31 * x[0] * x[2] - d3 * x[2];
    }


    /**********************************************************************************************//**
     * @fn	void operator()(const VectorBoost &x, MatrixBoost &M, Doub , VectorBoost &dfdt)
     *
     * @brief	The jacobian  needed for the stiff solver
     *
     * @author	Iron
     * @date	7/17/2018
     *
     * @param 		  	x		  	A VectorBoost to process.
     * @param [in,out]	M		  	A MatrixBoost to process.
     * @param 		  	parameter3	The third parameter.
     * @param [in,out]	dfdt	  	The dfdt.
     **************************************************************************************************/

    void operator()(const VectorBoost &x, MatrixBoost &M, Doub /*t*/, VectorBoost &dfdt)
	{

		Doub T = x[0];
		Doub H = x[1];
		Doub E = x[2];
		//The jacobian
			M(0, 0) = 1 - E * a13 - H * a12 - 2 * T;
			M(0, 1) = -T * a12;
			M(0, 2) = -T * a13;

			M(1, 0) = -H * a21;
			M(1, 1) = -T * a21 - H * r2 - r2 * (H - 1);
			M(1, 2) = 0;

			M(2, 0) = (E * r3) / (T + k3) - E * a31 - (E * T * r3) / pow((T + k3), 2);
			M(2, 1) = 0;
			M(2, 2) = (T * r3) / (T + k3) - T * a31 - d3;
			//since dx/dt = f, and f is not time dependent
			dfdt[0] = 0;
			dfdt[1] = 0;
			dfdt[2] = 0;

	}

	int getSystemSize() {
		return systemSize;
	}
	int getControlIndex() {
		return controlIndex;
	}
};

struct ItikBanksJacobian
{
	/**********************************************************************************************//**
	 * @property	Doub a12, a21, a13, a31
	 *
	 * @brief	Gets the 31
	 *
	 * @return	a 31.
	 **************************************************************************************************/

	Doub a12, a21, a13, a31;

	/**********************************************************************************************//**
	 * @property	Doub r2, r3, d3, k3
	 *
	 * @brief	Gets the k 3
	 *
	 * @return	The k 3.
	 **************************************************************************************************/

	Doub r2, r3, d3, k3;
	/** @brief	A Doub to process */
	Doub T;
	/** @brief	The height */
	Doub H;
	/** @brief	A Doub to process */
	Doub E;

	ItikBanksJacobian(std::vector<Doub> parameters,Doub T_,Doub H_, Doub E_):
		T(T_),
		H(H_),

		/**********************************************************************************************//**
		 * @fn	E(E_)
		 *
		 * @brief	Constructor
		 *
		 * @author	Iron
		 * @date	7/17/2018
		 *
		 * @param	parameter1	The first parameter.
		 **************************************************************************************************/

		E(E_)
	{
		a12 = parameters[0];
		a13 = parameters[1];
		r2 = parameters[2];
		a21 = parameters[3];
		r3 = parameters[4];
		k3 = parameters[5];
		a31 = parameters[6];
		d3 = parameters[7];
	}

	/**********************************************************************************************//**
	 * @fn	void operator()(const StateType &x, StateType &dxdt, Doub t)
	 *
	 * @brief	Function call operator
	 *
	 * @author	Iron
	 * @date	7/17/2018
	 *
	 * @param 		  	x   	A StateType to process.
	 * @param [in,out]	dxdt	The dxdt.
	 * @param 		  	t   	A Doub to process.
	 **************************************************************************************************/

	void operator()(const StateType &x, StateType &dxdt, Doub t)
	{
		dxdt[0] = (1 - E * a13 - H * a12 - 2 * T) * x[0] + (-T * a12)* x[3] + (-T * a13)* x[6];
		dxdt[1] = (1 - E * a13 - H * a12 - 2 * T) * x[1] + (-T * a12)* x[4] + (-T * a13)* x[7];
		dxdt[2] = (1 - E * a13 - H * a12 - 2 * T) * x[2] + (-T * a12)* x[5] + (-T * a13)* x[8];

		dxdt[3] = (-H * a21) * x[0] + (-T * a21 - H * r2 - r2 * (H - 1)) * x[3];
		dxdt[4] = (-H * a21) * x[1] + (-T * a21 - H * r2 - r2 * (H - 1)) * x[4];
		dxdt[5] = (-H * a21) * x[2] + (-T * a21 - H * r2 - r2 * (H - 1)) * x[5];

		dxdt[6] = ((E * r3) / (T + k3) - E * a31 - (E * T * r3) / pow((T + k3), 2))* x[0] + ((T * r3) / (T + k3) - T * a31 - d3)* x[6];
		dxdt[7] = ((E * r3) / (T + k3) - E * a31 - (E * T * r3) / pow((T + k3), 2))* x[1] + ((T * r3) / (T + k3) - T * a31 - d3)* x[7];
		dxdt[8] = ((E * r3) / (T + k3) - E * a31 - (E * T * r3) / pow((T + k3), 2))* x[2] + ((T * r3) / (T + k3) - T * a31 - d3)* x[8];
	}

};

// %Gets the jacobian following the procedure in eq (6)&(7) with state values at Fx (only one)
// Receive the function in row mayor order, and then return it has a matrix


//How to read paramters for this functon 

/**********************************************************************************************//**
 * @fn	std::vector<Doub> readParameters()
 *
 * @brief	Reads the parameters
 *
 * @author	Iron
 * @date	7/17/2018
 *
 * @return	The parameters.
 **************************************************************************************************/

std::vector<Doub> readParameters() {
	int numParam = 8;
	std::vector<Doub> param(numParam);
	for (Doub &i : param)
		std::cin >> i;
	return param;
}



#endif
