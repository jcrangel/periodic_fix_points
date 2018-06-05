#ifndef ITIKBANKS_H
#define ITIKBANKS_H


#include "all.h"


//Itikbanks model from Wei2014

struct itikBanks
{
    //The parameters
    double N;

	double a12, a21, a13, a31;
    double r2, r3, d3, k3;

    itikBanks(std::vector<double> parameters)
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
    void operator()(const stateType &x, stateType &dxdt, double /*t*/)
    {
        dxdt[0] = x[0] * (1 - x[0]) - a12 * x[0] * x[1] - a13 * x[0] * x[2];
        dxdt[1] = r2 * x[1] * (1 - x[1]) - a21 * x[0] * x[1];
        dxdt[2] = (r3 * x[0] * x[2]) / (k3 + x[0]) - a31 * x[0] * x[2] - d3 * x[2];
    }
};
struct itikBanks_stiff
{
    //The parameters
    double N;

	double a12, a21, a13, a31;
    double r2, r3, d3, k3;

    itikBanks_stiff(std::vector<double> parameters)
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

    void operator()(const stateType &x, stateType &dxdt, double /*t*/)
    {
        dxdt[0] = x[0] * (1 - x[0]) - a12 * x[0] * x[1] - a13 * x[0] * x[2];
        dxdt[1] = r2 * x[1] * (1 - x[1]) - a21 * x[0] * x[1];
        dxdt[2] = (r3 * x[0] * x[2]) / (k3 + x[0]) - a31 * x[0] * x[2] - d3 * x[2];
    }

    void operator()(const vectorBoost &x, vectorBoost &dxdt, double /*t*/)
    {
        dxdt[0] = x[0] * (1 - x[0]) - a12 * x[0] * x[1] - a13 * x[0] * x[2];
        dxdt[1] = r2 * x[1] * (1 - x[1]) - a21 * x[0] * x[1];
        dxdt[2] = (r3 * x[0] * x[2]) / (k3 + x[0]) - a31 * x[0] * x[2] - d3 * x[2];
    }

    //The jacobian
    void operator()(const vectorBoost &x, matrixBoost &M, double /*t*/, vectorBoost &dfdt)
	{
		double T = x[0];
		double H = x[1];
		double E = x[2];
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

			dfdt[0] = 0;
			dfdt[1] = 0;
			dfdt[2] = 0;

	}

};

//This jacobian is calculated over all one solution X(t)

struct itikBanksJacobianUt
{

    double a12, a21, a13, a31;
    double r2, r3, d3, k3;
    std::vector<stateType> state;
	stateType times;
    boost::math::barycentric_rational<double> T;
    boost::math::barycentric_rational<double> H;
    boost::math::barycentric_rational<double> E;

    itikBanksJacobianUt(std::vector<double> parameters, std::vector< stateType > xx, stateType tt)
                    :state(xx),
                    times(tt),
                    T(times.data(), state[0].data(), times.size()),
                    H(times.data(), state[1].data(), times.size()),
                    E(times.data(), state[2].data(), times.size())
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
	//Equation (7) Wei 2014
	// Row ordr
	// a b
	// c d
	// is returned as  [a b c d ] ^ T

	//
	//{ {D1, D2, D3},    { { x1, x2, x3 },
	//{ D4, D5, D6 }, X  { x4, x5, x6 },
	//{ D7, D8, D9 }}    { x7, x8, x9 }}
	//=
	//(D1 x1 + D2 x4 + D3 x7 | D1 x2 + D2 x5 + D3 x8 | D1 x3 + D2 x6 + D3 x9
	//	D4 x1 + D5 x4 + D6 x7 | D4 x2 + D5 x5 + D6 x8 | D4 x3 + D5 x6 + D6 x9
	//	D7 x1 + D8 x4 + D9 x7 | D7 x2 + D8 x5 + D9 x8 | D7 x3 + D8 x6 + D9 x9)
    void operator()(const stateType &x, stateType &dxdt, double t)
    {
		dxdt[0] = (1 - E(t) * a13 - H(t) * a12 - 2 * T(t)) * x[0] + (-T(t) * a12)* x[3] + (-T(t) * a13)* x[6];
		dxdt[1] = (1 - E(t) * a13 - H(t) * a12 - 2 * T(t)) * x[1] + (-T(t) * a12)* x[4] + (-T(t) * a13)* x[7];
		dxdt[2] = (1 - E(t) * a13 - H(t) * a12 - 2 * T(t)) * x[2] + (-T(t) * a12)* x[5] + (-T(t) * a13)* x[8];

		dxdt[3] = (-H(t) * a21) * x[0] + (-T(t) * a21 - H(t) * r2 - r2 * (H(t) - 1)) * x[3];
		dxdt[4] = (-H(t) * a21) * x[1] + (-T(t) * a21 - H(t) * r2 - r2 * (H(t) - 1)) * x[4];
		dxdt[5] = (-H(t) * a21) * x[2] + (-T(t) * a21 - H(t) * r2 - r2 * (H(t) - 1)) * x[5];

		dxdt[6] = ((E(t) * r3) / (T(t) + k3) - E(t) * a31 - (E(t) * T(t) * r3) / pow((T(t) + k3), 2))* x[0] + ((T(t) * r3) / (T(t) + k3) - T(t) * a31 - d3)* x[6];
		dxdt[7] = ((E(t) * r3) / (T(t) + k3) - E(t) * a31 - (E(t) * T(t) * r3) / pow((T(t) + k3), 2))* x[1] + ((T(t) * r3) / (T(t) + k3) - T(t) * a31 - d3)* x[7];
		dxdt[8] = ((E(t) * r3) / (T(t) + k3) - E(t) * a31 - (E(t) * T(t) * r3) / pow((T(t) + k3), 2))* x[2] + ((T(t) * r3) / (T(t) + k3) - T(t) * a31 - d3)* x[8];
    }

	//The jacobian
	//matrix<double> Dfu(double t) {
	//	matrix<double> M(STATE_SIZE, STATE_SIZE);
	//	M(0, 0) = 1 - E(t) * a13 - H(t) * a12 - 2 * T(t);
	//	M(0, 1) = -T(t) * a12;
	//	M(0, 2) = -T(t) * a13;

	//	M(1, 0) = -H(t) * a21;
	//	M(1, 1) = -T(t) * a21 - H(t) * r2 - r2 * (H(t) - 1);
	//	M(1, 2) = 0;

	//	M(2, 0) = (E(t) * r3) / (T(t) + k3) - E(t) * a31 - (E(t) * T(t) * r3) / pow((T(t) + k3), 2);
	//	M(2, 1) = 0;
	//	M(2, 2) = (T(t) * r3) / (T(t) + k3) - T(t) * a31 - d3;
	//}
};

//This ja

struct itikBanksJacobian
{

	double a12, a21, a13, a31;
	double r2, r3, d3, k3;
	double T;
	double H;
	double E;

	itikBanksJacobian(std::vector<double> parameters,double T_,double H_, double E_):
		T(T_),
		H(H_),
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
	void operator()(const stateType &x, stateType &dxdt, double t)
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

Matrix3d DFitikBanks(stateType Fx,double tau){

    //std::vector<stateType> V;
    //std::vector<double> t;
	stateType xp={0,tau};

    itikBanksJacobian Dfu(PARAMETERS,Fx[0],Fx[1],Fx[2]);
    typedef runge_kutta_dopri5<stateType> error_stepper_type;


	stateType identityMatrixVector = {1,0,0,0,1,0,0,0,1};
	size_t steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6),
		Dfu, identityMatrixVector, xp[0], xp[1], 0.001);
            //push_back_state_and_time(V, t));
    // Returns the last matrix since
    // DF(x) = v(tau).See page 12 research notebook
   // return reshapeVectorToMatrix(V.back());
	return reshapeVectorToMatrix(identityMatrixVector);
}
//Version where a stiff solver is used
// Matrix3d DFitikBanks_stiff(stateType Fx,double tau){

//     //std::vector<stateType> V;
//     //std::vector<double> t;
// 	stateType xp={0,tau};

//     itikBanksJacobian Dfu(PARAMETERS,Fx[0],Fx[1],Fx[2]);
//     typedef runge_kutta_dopri5<stateType> error_stepper_type;


// 	vectorBoost identityMatrixVector(8);
// 	identityMatrixVector[0] =1;
// 	identityMatrixVector[1] =0;
// 	identityMatrixVector[2] =0;
// 	identityMatrixVector[3] =0;
// 	identityMatrixVector[4] =1;
// 	identityMatrixVector[5] =0;
// 	identityMatrixVector[6] =0;
// 	identityMatrixVector[7] =1;

// 	size_t num_of_steps = integrate_const(make_dense_output< rosenbrock4< double > >(1.0e-6, 1.0e-6),
// 						std::make_pair(functionName, functionName),identityMatrixVector, 0.0, tau, 0.01);

// 	size_t steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6),
// 		Dfu, identityMatrixVector, xp[0], xp[1], 0.001);
//             //push_back_state_and_time(V, t));
//     // Returns the last matrix since
//     // DF(x) = v(tau).See page 12 research notebook
//    // return reshapeVectorToMatrix(V.back());
// 	return reshapeVectorToMatrix(identityMatrixVector);
// }



#endif
