/*
Util.h

TODO: All common headers, google styles says this is not good

Common data type definition and functions for vector manipulation.

*/


#ifndef ALL_H
#define ALL_H



#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
//Boosst lib
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>

#define STATE_SIZE 3
#define CONTROL_POS 2
#define DEBUG0 true
#define DEBUG1 false

//const std::vector<double> PARAMETERS = { 1,2.5,0.5,1.5,4.5,1,0.2,0.5 };
//This should not be global
//std::vector<double> PARAMETERS;

using namespace Eigen;
using namespace boost::numeric::odeint;

//An array [x1(t),x2(t),x3(t)]
//typedef std::array< double, STATE_SIZE > stateType;

typedef std::vector< double > stateType;
typedef boost::numeric::ublas::vector< double > vectorBoost;
typedef boost::numeric::ublas::matrix< double > matrixBoost;



//The class representing a fixpoint [x,y,z] and if info about 
//convergence and stability
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

template <class T>
Matrix3d reshapeVectorToMatrix(const T x)
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
bool equal(double A, double B, double epsilon = 0.000005f)
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



// Check if the fixpoint is in the Set S 
bool pointIsInSet(fixPoint p, std::vector<fixPoint> S)
{
	for (fixPoint i: S) {
		//Check
		int j;
		for (j = 0; j < STATE_SIZE; j++) {

			if ( !equal(i.solution[j],p.solution[j]) )
				break;	//Stop comparing, vectors ain't equal
		}
		if (j == STATE_SIZE) // never break, then all were equal
			return true;

	}
	return false;
}

bool pointHaveNegatives(fixPoint p) {
	for (int i = 0; i < p.solution.size(); i++)
		if (p.solution[i] < 0)
			return true;

	return false;
}

/*
Warper for writing both to the console and to a file.

work only for one << at time
TODO: this dont work with std::endl
*/
class LogAndStdout {
private:
	std::ofstream file;
public:

	LogAndStdout(const std::string fileName) {
		file.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);
	}

	template<class T>
	LogAndStdout& operator<<(const T data) {
		std::cout << data;
		file << data;

		return *this; // TODO : This clone the object? if this is true, then this is bad 
	}
};


//%Evaluate F(x,y,z)  at the final time.  
template <class T>
Vector3d evalFunInLast(T &functionName, stateType initialCondition, double tau, double d)
{
	initialCondition[CONTROL_POS] = initialCondition[CONTROL_POS] + d;

	std::vector<double> res(STATE_SIZE);
	integrateSystem(functionName, initialCondition, res, 0, tau);

	Vector3d v(res.data());
	return v;
	//initialCondition has the last
}


#endif
