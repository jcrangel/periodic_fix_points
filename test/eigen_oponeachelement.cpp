#include <cmath>
#include <iostream>

#include <Eigen/Core>

double Exp(double x) // the functor we want to apply
{
    return std::exp(x);
}

int main()
{
    Eigen::Vector3f m= {1,2,3,4,5};
    
    std::cout << m << std::endl << "becomes: ";
    std::cout << std::endl << m.unaryExpr(&Exp) << std::endl;
 }