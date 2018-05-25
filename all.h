#ifndef ALL_H
#define ALL_H

#include <iostream>
#include <cmath>
//Boosst lib
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>

#define STATE_SIZE 3

using namespace Eigen;
using namespace boost::numeric::odeint;

//An array [x1(t),x2(t),x3(t)]	
//typedef std::array< double, STATE_SIZE > stateType;
typedef std::vector< double > stateType;

//[ integrate_observer
	struct push_back_state_and_time
	{
	    std::vector< stateType >& m_states;
	    std::vector< double >& m_times;

	    push_back_state_and_time( std::vector< stateType > &states , std::vector< double > &times )
	    : m_states( states ) , m_times( times ) { }

	    void operator()( const stateType &x , double t )
	    {
	        m_states.push_back( x );
	        m_times.push_back( t );
	    }
	};

//Transpose
void transpose(const std::vector < stateType> u,
					 std::vector < std::vector<double> > & state ){
    //Copy the data as transpose
	//for some stupid reason in windows we need i < u.size in ubuntu ins i <= u.size
    for(int i = 0; i < u.size() ; i++)  
    	for (int j = 0; j < STATE_SIZE ; j++ )
    		state[j][i] = u[i][j];

}



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
    void operator()(const stateType &x, stateType &dxdt, double t)
    {
        dxdt[0] = x[0] * (1 - x[0]) - a12 * x[0] * x[1] - a13 * x[0] * x[2];
        dxdt[1] = r2 * x[1] * (1 - x[1]) - a21 * x[0] * x[1];
        dxdt[2] = (r3 * x[0] * x[2]) / (k3 + x[0]) - a31 * x[0] * x[2] - d3 * x[2];
    }
};

//Reshape a vector Size N^2 into a Matrix NxN, 
//vector is assumed to come in row order
Matrix3d reshapeVectorToMatrix(const stateType x) {

	int N = sqrt(x.size());
	Matrix3d A(N,N);
	
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A(i, j) = x[i*N + j];
		}
	}

	return A;
}


struct itikBanksJacobian
{

    double a12, a21, a13, a31; 
    double r2, r3, d3, k3;
    std::vector<std::vector<double>> state;
    std::vector<double> times;
    boost::math::barycentric_rational<double> T;
    boost::math::barycentric_rational<double> H;
    boost::math::barycentric_rational<double> E;

    itikBanksJacobian(std::vector<double> parameters, std::vector< std::vector<double> > xx, std::vector<double> tt) 
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

//Gets the jacobian following the procedure in eq (6)&(7) of the
//paper Wei 2014
template <class T>
Matrix3d DFode(T &fun,stateType initialCondition,double tau,double d)
{
// xp=[0, tau];
	stateType xp{0,tau};

	initialCondition[2] = initialCondition[2] + d;
	std::vector<stateType> u;
    std::vector<double> tt;
    typedef runge_kutta_cash_karp54<stateType> error_stepper_type;
    size_t steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6),
                   fun, initialCondition, xp[0], xp[1], 0.001,
                   push_back_state_and_time(u, tt));
    // u save the data column wise, the firts colum has the first state data. It would
    //be nice to have in row 
    //This create a transpose of u, with dimentions: STATE_SIZE * steps (3 x steps)
    // with elements equal to zero
    std::vector < std::vector < double>> state(STATE_SIZE, std::vector<double>(u.size()));
    transpose(u,state);
     //Check if its correct
    //  for( size_t i=0; i <= steps; i++ )
    //  {
    //     std::cout << tt[i] << '\t' << state[0][i] << '\t'
    //     	     << state[1][i] << '\t' << state[2][i] << '\n';
    //  }
    std::vector<stateType> V;
    std::vector<double> t;
	//!!!! this param shoul exist only on one place not to be repated
    std::vector<double> param = { 1,2.5,0.5,1.5,4.5,1,0.2,0.5 };
    std::vector<double> identityMatrixVector = {1,0,0,0,1,0,0,0,1};
    itikBanksJacobian Dfu(param,state,tt);
    steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6), 
            Dfu, identityMatrixVector, xp[0], xp[1], 0.001,
            push_back_state_and_time(V, t));
 //   for( size_t i=0; i <= steps; i++ )
	//{
 //       std::cout << V[i][0] << '\t' << V[i][1] << '\t'
 //                 << V[i][2] << '\t' << V[i][3] << '\t' << V[i][4] << '\t'
 //                 << V[i][5] << '\t' << V[i][6] << '\t' << V[i][7] << '\t'
 //                 << V[i][8] << '\n';
 //   }

	// Returns the last matrix since 
	// DF(x) = v(tau).See page 12 research notebook
	return reshapeVectorToMatrix(V.back());
};

class fixPoint {
public:
	bool convergent;
	bool stability;
	Vector3d solution;

	fixPoint(bool convergent_ ,bool stability_,Vector3d solution_) :
		convergent(convergent_),
		stability(stability_),
		solution(solution_) {}

};
//Calculate all eigenvalues of matrix A, if all of them are less tahn 1
// the matrix is stable()true.
bool isStable(MatrixXd A) {
	VectorXcd eivals = A.eigenvalues();
	for (int i = 0; i < eivals.size(); ++i)
		if (std::abs(eivals[i]) >= 1)
			return false;

	return true;
};

template <class T>
fixPoint newtonPoincare(T &fun, stateType initialCondition, double tau, double d) {
	int n = 0;            //%initialize iteration counter
	double eps = 1;          //%initialize error
	std::vector<double> xp{ 0, tau }; //% define the span of the computational domain
	double tol = 1e-6;
	int maxStep = 100;

	Vector3d Xk(initialCondition.data());
	Vector3d Pxk;
	Matrix3d df;
	std::vector<stateType> u;
	std::vector<double> t;
	Matrix3d A;
	Matrix3d I = Matrix3d::Identity();
	Vector3d y;
	Vector3d Xplus;
	size_t steps;
	typedef runge_kutta_cash_karp54<stateType> error_stepper_type;

	//stateType initCond;

	while (n < maxStep) {					//!@_@
		df = DFode(fun, std::vector<double> {Xk[0], Xk[1], Xk[2]}, tau, d);  //Jacobian
		initialCondition[0] = Xk[0];
		initialCondition[1] = Xk[1];
		initialCondition[2] = Xk[2]+d;

	    steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6),
			fun,				//!@_@							
			initialCondition,
			xp[0], xp[1], 0.001, push_back_state_and_time(u, t));
		//[t, u] = ode15s(functionName, xp, [Xk(1), Xk(2), Xk(3) + d], [], [parameters d]);
		
		//!@This should be done in a bette wway
		Pxk[0] = u.back()[0];
		Pxk[1] = u.back()[1];
		Pxk[2] = u.back()[2];

		A = df - I;
		//Solve Ay = (Pxk-Xk)
		y = A.colPivHouseholderQr().solve(Pxk-Xk);
		Xplus = Xk - y;
		
		eps = y.norm();
		if (eps < tol)
			break;
	
		Xk = Xplus;  //%update x
		n++;
		
	}
	// 0 means unstable
	double stability = 0;
	bool convergence = false;
	if (n < maxStep){
		convergence = true;
		
		if (isStable(df))
			stability = 1;

	}
	fixPoint ss(convergence, stability, Xk);
	return ss;
}


#endif