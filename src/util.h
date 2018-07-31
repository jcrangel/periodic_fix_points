/*
Util.h

TODO: All common headers, google styles says this is not good

Common data type definition and functions for vector manipulation.

*/


#ifndef ALL_H
///< .
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

#define STATE_SIZE 3 //TODO: remove this
#define CONTROL_POS 2
#define DEBUG0 true
#define DEBUG1 false

//const std::vector<double> PARAMETERS = { 1,2.5,0.5,1.5,4.5,1,0.2,0.5 };
//This should not be global
//std::vector<double> PARAMETERS;

/**********************************************************************************************//**
 * @namespace	Eigen
 *
 * @brief	.
 **************************************************************************************************/

using namespace Eigen;

/**********************************************************************************************//**
 * @namespace	boost::numeric::odeint
 *
 * @brief	.
 **************************************************************************************************/

using namespace boost::numeric::odeint;

//An array [x1(t),x2(t),x3(t)]
//typedef std::array< double, STATE_SIZE > stateType;

/**********************************************************************************************//**
 * @typedef	std::vector< double > stateType
 *
 * @brief	Defines an alias representing type of the state
 **************************************************************************************************/

typedef std::vector< double > stateType;

/**********************************************************************************************//**
 * @typedef	boost::numeric::ublas::vector< double > vectorBoost
 *
 * @brief	Defines an alias representing the vector boost
 **************************************************************************************************/

typedef boost::numeric::ublas::vector< double > vectorBoost;

/**********************************************************************************************//**
 * @typedef	boost::numeric::ublas::matrix< double > matrixBoost
 *
 * @brief	Defines an alias representing the matrix boost
 **************************************************************************************************/

typedef boost::numeric::ublas::matrix< double > matrixBoost;

/**********************************************************************************************//**
 * @typedef	const std::vector<double> vecDoub_I
 *
 * @brief	Defines an alias representing the vector doub i
 **************************************************************************************************/

typedef const std::vector<double> vecDoub_I; // Utility vector for input

/**********************************************************************************************//**
 * @typedef	std::vector<double> vecDoub_IO
 *
 * @brief	Defines an alias representing the vector doub i/o
 **************************************************************************************************/

typedef std::vector<double>  vecDoub_IO;// Utility vector for input & ouput

//The class representing a fixpoint [x,y,z] and if info about 
//convergence and stability

/**********************************************************************************************//**
 * @class	fixPoint
 *
 * @brief	A fix point.
 *
 * @author	Iron
 * @date	7/25/2018
 **************************************************************************************************/

class fixPoint
{
public:
    /** @brief	True to convergent */
    bool convergent;
    /** @brief	True to stability */
    bool stability;
    /** @brief	The solution */
    VectorXd solution;

    fixPoint(bool convergent_,bool stability_,VectorXd solution_) :
        convergent(convergent_),
        stability(stability_),

        /**********************************************************************************************//**
         * @fn	fixPoint::solution(solution_)
         *
         * @brief	Constructor
         *
         * @author	Iron
         * @date	7/25/2018
         *
         * @param	parameter1	The first parameter.
         **************************************************************************************************/

        solution(solution_) {}

};

//typedef std::vector<fixPoint> fixPoints;

//Transpose a matrix of made with std::vector's
template <class T>
void transpose(const T u,
               /** @brief	The ) */
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

/**********************************************************************************************//**
 * @fn	Matrix3d reshapeVectorToMatrix(const T x)
 *
 * @brief	Reshape vector to matrix
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	x	A T to process.
 *
 * @return	A Matrix3d.
 **************************************************************************************************/

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

/**********************************************************************************************//**
 * @fn	bool isStable(MatrixXd A)
 *
 * @brief	Query if 'A' is stable
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	A	A MatrixXd to process.
 *
 * @return	True if stable, false if not.
 **************************************************************************************************/

bool isStable(MatrixXd A)
{
    VectorXcd eivals = A.eigenvalues();
    for (int i = 0; i < eivals.size(); ++i)
        if (std::abs(eivals[i]) >= 1)
            return false;

    return true;
};


//Cheking for equality,

/**********************************************************************************************//**
 * @fn	bool equal(double A, double B, double epsilon = 0.000005f)
 *
 * @brief	Equals
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	A	   	A double to process.
 * @param	B	   	A double to process.
 * @param	epsilon	(Optional) The epsilon.
 *
 * @return	True if it succeeds, false if it fails.
 **************************************************************************************************/

bool equal(double A, double B, double epsilon = 0.000005f)
{
    return (fabs(A - B) < epsilon);
}
//get a std vector from eigen Vector3d

/**********************************************************************************************//**
 * @fn	std::vector<double> toStdVectorD(const Vector3d v)
 *
 * @brief	Converts a v to a standard vector d
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	v	A Vector3d to process.
 *
 * @return	V as a std::vector&lt;double&gt;
 **************************************************************************************************/

stateType toStateType(const VectorXd v)
{
    std::vector<double> v2;
    v2.resize(v.size());
    VectorXd::Map(&v2[0], v.size()) = v;
    return v2;
}
//Copy a Vector3d into a std vector

/**********************************************************************************************//**
 * @fn	void toStdVectorD(const Vector3d v, stateType &w)
 *
 * @brief	Converts this object to a standard vector d
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param 		  	v	A Vector3d to process.
 * @param [in,out]	w	A stateType to process.
 **************************************************************************************************/

void toStateType(const VectorXd v, stateType &w)
{
	w.resize(v.size());
	VectorXd::Map(&w[0], v.size()) = v;
}

/**********************************************************************************************//**
 * @fn	std::vector<double> toStdVectorD(vectorBoost v)
 *
 * @brief	Converts a v to a standard vector d
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	v	A vectorBoost to process.
 *
 * @return	V as a std::vector&lt;double&gt;
 **************************************************************************************************/

stateType toStateType(vectorBoost v)
{
	std::vector<double> w(v.size());
	std::copy(v.begin(), v.end(), w.begin());
	return w;
}

/**********************************************************************************************//**
 * @fn	vectorBoost toBoostVectorD(const stateType v)
 *
 * @brief	Converts a v to a boost vector d
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	v	A stateType to process.
 *
 * @return	V as a vectorBoost.
 **************************************************************************************************/

vectorBoost toBoostVectorD(const stateType v)
{
	vectorBoost x(v.size());
	std::copy(v.begin(), v.end(), x.begin());
	return x;
}

/**********************************************************************************************//**
 * @fn	Vector3d toEigenVector(const stateType v)
 *
 * @brief	Converts a v to an eigen vector
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	v	A stateType to process.
 *
 * @return	V as a Vector3d.
 **************************************************************************************************/

VectorXd toEigenVector(stateType v) {
	//double* ptr = &v[0];
	//Eigen::Map< Eigen::VectorXd> v2(ptr, v.size());
	////Vector3d v2(v.data());
	//return v2;
	return VectorXd::Map(v.data(), v.size());
}



// Check if the fixpoint is in the Set S 

/**********************************************************************************************//**
 * @fn	bool pointIsInSet(fixPoint p, std::vector<fixPoint> S)
 *
 * @brief	Point is in set
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	p	A fixPoint to process.
 * @param	S	A std::vector&lt;fixPoint&gt; to process.
 *
 * @return	True if it succeeds, false if it fails.
 **************************************************************************************************/

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

/**********************************************************************************************//**
 * @fn	bool pointHaveNegatives(fixPoint p)
 *
 * @brief	Point have negatives
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	p	A fixPoint to process.
 *
 * @return	True if it succeeds, false if it fails.
 **************************************************************************************************/

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

/**********************************************************************************************//**
 * @class	LogAndStdout
 *
 * @brief	A log and stdout.
 *
 * @author	Iron
 * @date	7/25/2018
 **************************************************************************************************/

class LogAndStdout {
private:
	/** @brief	The file */
	std::ofstream file;
public:

	/**********************************************************************************************//**
	 * @fn	TODO::LogAndStdout(const std::string fileName)
	 *
	 * @brief	Constructor
	 *
	 * @author	Iron
	 * @date	7/25/2018
	 *
	 * @param	fileName	Filename of the file.
	 **************************************************************************************************/

	LogAndStdout(const std::string fileName) {
		file.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);
	}

	template<class T>
///< .
	LogAndStdout& operator<<(const T data) {
///< .
		std::cout << data;
///< .
		file << data;

		/** @brief	TODO : This clone the object? if this is true, then this is bad */
		return *this;
	}
};


//%Evaluate F(x,y,z)  at the final time.  


/**********************************************************************************************//**
 * @fn	Vector3d evalFunInLast(T &functionName, stateType initialCondition, double tau, double d)
 *
 * @brief	Eval fun in last
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param [in,out]	functionName		Name of the function.
 * @param 		  	initialCondition	The initial condition.
 * @param 		  	tau					The tau.
 * @param 		  	d					A double to process.
 *
 * @return	A Vector3d.
 **************************************************************************************************/
template <class T>
VectorXd evalFunInLast(T &functionName, stateType initialCondition, double tau, double d)
{
	initialCondition[CONTROL_POS] = initialCondition[CONTROL_POS] + d;
	int N = functionName.getSystemSize();
	std::vector<double> res(N);
	integrateSystem(functionName, initialCondition, res, 0, tau);

	//To create the VectorXd from the std::vector
	double* ptr = &res[0];
	Eigen::Map<Eigen::VectorXd> v(ptr, N);

	return v;
	//initialCondition has the last
}


/**
Create the cartesian set A x B
*/
template<class T>

/**********************************************************************************************//**
 * @fn	std::vector< std::vector<T> > cartesianProduct(const std::vector<T> A, const std::vector<T> B)
 *
 * @brief	Cartesian product
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	A	A std::vector&lt;T&gt; to process.
 * @param	B	A std::vector&lt;T&gt; to process.
 *
 * @return	A std::vector&lt;std::vector&lt;T&gt; &gt;
 **************************************************************************************************/

std::vector< std::vector<T> > cartesianProduct(const std::vector<T> A, const std::vector<T> B)
{
	std::vector<std::vector<T> > prod;
	for (int i = 0; i < A.size(); i++)
		for (int j = 0; j < B.size(); j++)
			prod.push_back(std::vector<T> {A[i], B[j]});

	return prod;
}
/**

*/

/**********************************************************************************************//**
 * @fn	template<class T> std::vector< std::vector<T> > cartesianProduct(const std::vector<std::vector<T>> A, const std::vector<T> B)
 *
 * @brief	Create the cartesian set A x B = (a x b) x B, where A is already a cartesian set.
 * 			Then to create the cartesian set of {x,y,z} x {x,y,z} x {x,y,z}. We have to : s = {x,
 * 			y,z};
 * 			std::vector&lt;std::vector&lt;double&gt; &gt; prod = cartesianProduct(s,s);
 * 			prod = cartesianProduct(prod,s);
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @tparam	T	Generic type parameter.
 * @param	A	A std::vector&lt;std::vector&lt;T&gt;&gt; to process.
 * @param	B	A std::vector&lt;T&gt; to process.
 *
 * @return	A std::vector&lt;std::vector&lt;T&gt; &gt;
 **************************************************************************************************/

template<class T>
std::vector< std::vector<T> > cartesianProduct(const std::vector<std::vector<T>> A, const std::vector<T> B)
{
	std::vector<std::vector<T> > prod;
	for (int i = 0; i < A.size(); i++)
		for (int j = 0; j < B.size(); j++) {
			std::vector<T> AB;
			AB.reserve(A.size() + B.size());
			AB.insert(AB.end(), A[i].begin(), A[i].end());
			AB.insert(AB.end(), B[j]);

			prod.push_back(AB);

		}
	return prod;
}


template<class T>
void printVector(std::vector<T> M)
{
		for (T j : M)
			std::cout << j << " ";
		std::cout << "\n";
}

/**********************************************************************************************//**
 * @fn	template<class T> void printVectorVector(std::vector< std::vector<T> > M)
 *
 * @brief	Print vector of vectors
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @tparam	T	Generic type parameter.
 * @param	M	A std::vector&lt;std::vector&lt;T&gt;&gt; to process.
 **************************************************************************************************/

template<class T>
void printVectorVector(std::vector< std::vector<T> > M)
{

	for (std::vector<T> i : M) {
		for (T j : i)
			std::cout << j << " ";
		std::cout << "\n";
	}

}





#endif
