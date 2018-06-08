#ifndef ALL_H
#define ALL_H

#include <iostream>
#include <cmath>
//Boosst lib
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>

#define STATE_SIZE 3
#define CONTROL_POS 2
#define DEBUG false

const std::vector<double> PARAMETERS = { 1,2.5,0.2,1.5,4.5,1,0.2,0.5 };

using namespace Eigen;
using namespace boost::numeric::odeint;

//An array [x1(t),x2(t),x3(t)]
//typedef std::array< double, STATE_SIZE > stateType;
typedef std::vector< double > stateType;
typedef boost::numeric::ublas::vector< double > vectorBoost;
typedef boost::numeric::ublas::matrix< double > matrixBoost;


//The class represents eigenvector
class fixPoint
{
public:
    bool convergent;
    bool stability;
    Vector3d solution;

    fixPoint(bool convergent_,bool stability_,Vector3d solution_) :
        convergent(convergent_),
        stability(stability_),
        solution(solution_) {}

};
//typedef std::vector<fixPoint> fixPoints;

//Transpose a matrix of made with std::vector's
template <class T>
void transpose(const T u,
               std::vector < stateType > & state )
{
    //Copy the data as transpose
    //for some stupid reason in windows we need i < u.size in ubuntu ins i <= u.size
    for(int i = 0; i < u.size() ; i++)
        for (int j = 0; j < STATE_SIZE ; j++ )
            state[j][i] = u[i][j];

}





//Reshape a vector Size N^2 into a Matrix NxN,
//vector is assumed to come in row order

Matrix3d reshapeVectorToMatrix(const stateType x)
{

    int N = sqrt(x.size());
    Matrix3d A(N,N);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A(i, j) = x[i*N + j];
        }
    }

    return A;
}


//Calculate all eigenvalues of matrix A, if all of them are less tahn 1
// the matrix is stable()true.
bool isStable(MatrixXd A)
{
    VectorXcd eivals = A.eigenvalues();
    for (int i = 0; i < eivals.size(); ++i)
        if (std::abs(eivals[i]) >= 1)
            return false;

    return true;
};


//Cheking for equality,
bool equal(double A, double B, double epsilon = 0.005f)
{
    return (fabs(A - B) < epsilon);
}
//get a std vector from eigen Vector3d
std::vector<double> toStdVectorD(const Vector3d v)
{
    std::vector<double> v2;
    v2.resize(v.size());
    VectorXd::Map(&v2[0], v.size()) = v;
    return v2;
}
//Copy a Vector3d into a std vector
void toStdVectorD(const Vector3d v, stateType &w)
{
	w.resize(v.size());
	VectorXd::Map(&w[0], v.size()) = v;
}
std::vector<double> toStdVectorD(vectorBoost v)
{
	std::vector<double> w(v.size());
	std::copy(v.begin(), v.end(), w.begin());
	return w;
}


vectorBoost toBoostVectorD(const stateType v)
{
	vectorBoost x(v.size());
	std::copy(v.begin(), v.end(), x.begin());
	return x;
}

Vector3d toEigenVector(const stateType v) {
	Vector3d v2(v.data());
	return v2;
}

//Integrate a ODE system and returns the last value of the integration in x
//In the stiff version the system should be a class that have the Jacobian as a functor
template <class T>
void integrateStiffSystem(T &system, stateType initialCondition,stateType &x,double t0,double tf){

	vectorBoost v = toBoostVectorD(initialCondition);;

	//std::copy(initialCondition.begin(), initialCondition.end(), v.begin());

	size_t num_of_steps = integrate_const(make_dense_output< rosenbrock4< double > >(1.0e-6, 1.0e-6),
		std::make_pair(system, system),
		v, t0, tf, 0.01);

		//std::vector<double> res(STATE_SIZE);
		std::copy(v.begin(), v.end(), x.begin());

}

//Integrate a ODE system and returns the last value of the integration in x
//In the stiff version the system should be a class that have the Jacobian as a functor
//save the data into u and t
template <class T>
void integrateStiffSystem(T &system, stateType initialCondition,
	 double t0, double tf, std::vector<stateType> &u, stateType &t) {

	vectorBoost v = toBoostVectorD(initialCondition);

	//std::copy(initialCondition.begin(), initialCondition.end(), v.begin());

	size_t num_of_steps = integrate_const(make_dense_output< rosenbrock4< double > >(1.0e-6, 1.0e-6),
		std::make_pair(system, system),
		v, t0, tf, 0.01, push_back_state_and_time(u, t));

}


//[ integrate_observer
struct push_back_state_and_time
{
	std::vector< stateType>& m_states;
	std::vector<double> & m_times;

	push_back_state_and_time(std::vector< stateType > &states, std::vector<double> &times)
		: m_states(states), m_times(times) { }

	void operator()(const stateType &x, double t)
	{
		m_states.push_back(x);
		m_times.push_back(t);
	}

	void operator()(const vectorBoost &x, double t)
	{
		m_states.push_back(toStdVectorD(x));
		m_times.push_back(t);
	}
};

#endif
