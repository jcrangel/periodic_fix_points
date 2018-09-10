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
#include <map>
#include <iomanip>
#include <sstream>

String keystr(Doub a, Doub b, int precision = 6){
	strings<<tream stream;
	stream << fixed << setprecision(precision) <<
	 String as = std::to_string(a);
	String bs =  std::to_string(b);

	return
}
// function key = keystr(a, b)
// 	key = rounddecimals(hash(a, b), 4);
// key = num2str(key);

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
	// std::cout << 10e-4;

	std::map<double,int> a;

	a.insert(std::pair<double,int>(3.1416,3) );
	a.insert(std::pair<double,int>(3.14,32) );
	a.insert(std::pair<double,int>(4.5,76) );

	std::cout<<a[4.5]<<std::endl;
	std::cout << a[3.14] << std::endl;
	std::cout << a[3.1416] << std::endl;

	std::map<std::pair<double,double>, int> b;

	b.insert(std::pair< std::pair<double,double >, int>( , 3));

}