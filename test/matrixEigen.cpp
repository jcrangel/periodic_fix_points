#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
int main()
{
  Matrix3d m = Matrix3d::Identity(3,3);
  m = MatrixXd::Constant(3,3,1.2);
  cout << "m =" << endl << m << endl;
  cout << "I =" << endl << Matrix3d::Identity(3,3)<< endl;

  VectorXd v(3);
  v << 1, 2, 3;
  cout << "m * v =" << endl << m * v << endl;
}