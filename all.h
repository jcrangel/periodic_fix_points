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

const std::vector<double> PARAMETERS = { 1,2.5,0.5,1.5,4.5,1,0.2,0.5 };

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



//[ integrate_observer

template <class T>
struct push_back_state_and_time
{
    std::vector< T >& m_states;
    std::vector<double> & m_times;

    push_back_state_and_time( std::vector< T > &states, std::vector<double> &times )
        : m_states( states ), m_times( times ) { }

    void operator()( const T &x, double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};
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

std::vector<double> toStdVectorD(Vector3d v)
{
    std::vector<double> v2;
    v2.resize(v.size());
    VectorXd::Map(&v2[0], v.size()) = v;
    return v2;
}


#endif
