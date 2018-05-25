/*
To do:
1. Define StateTye as Vector3d so is congruent with all the program.
*/

#include "all.h"


int main(){

	stateType IC = {0,1,0.06};
	std::vector<double> param={1,2.5,0.5,1.5,4.5,1,0.2,0.5};
	itikBanks fun(param);
	Matrix3d m = DFode(fun,IC,10,0.5);
	std::cout << m << std::endl;

	fixPoint a = newtonPoincare(fun, IC, 10, 0.5);

	std::cout <<"Convergent?:\n"<< a.convergent <<"\n Sol: \n"<< a.solution <<"\n Stable?:\n" << a.stability<<std::endl;
	std::cin.get();//VS windows only
	return 0;
}