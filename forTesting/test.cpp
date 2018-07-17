#include <cmath>
#include <iostream>

#include <Eigen/Core>

double Exp(double x) // the functor we want to apply
{
	return std::exp(x);
}

int main()
{
	Eigen::Vector3f m(1,2,3);

	std::cout << m << std::endl << "becomes: ";
	std::cout << std::endl << m.unaryExpr(&Exp) << std::endl;
}