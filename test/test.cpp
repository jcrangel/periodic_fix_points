#include <cmath>
#include <iostream>
// #include <complex>

// #include <Eigen/Core>
// #include <Eigen/Eigenvalues> 

// using namespace Eigen;
// double Exp(double x) // the functor we want to apply
// {
// 	return std::exp(x);
// }

// bool isStable(MatrixXd A){
// 	VectorXcd eivals = A.eigenvalues();
// 	for (int i = 0; i < eivals.size(); ++i)
// 		if (std::abs(eivals[i]) >= 1)
// 			return false;
	
// 	return true;
// }

int main()
{
	// Eigen::Vector4d A (1,2,3,4);
	
	// std::cout << A << std::endl << "becomes: ";
	// std::cout << std::endl << A.unaryExpr(&Exp) << std::endl;

	// MatrixXd ones(2,2);
	// ones << 0.5, 0, 0, 0.2;
	// VectorXcd eivals = ones.eigenvalues();
	// //std::cout << "The eigenvalues of the 3x3 matrix of ones are:" 
	// //		  << std::endl << std::abs(eivals[0]) << std::endl;
	// //std::cout << "The absolute value are" 
	// //	<< std::endl << eivals.unaryExpr([](std::complex<double> z) {return std::abs(z); })
	// //	<< std::endl;

	// std::cout << "Da matrix is " << (isStable(ones) ? "Stable" : "Unstable");

	// std::cin.get();
	std::cout << 10e-4;
}