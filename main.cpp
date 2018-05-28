/*
To do:
1. Define StateTye as Vector3d so is congruent with all the program.
    -- Create a function that integrate inside use std::vector outside return Vector3d
    --
2. Theres lots code chucks were a ode is solve it, but only the last value is needed. And save all
the state values is unnecesary

Symbols:
        #_# means that part of the should be rewritten to support systems of several size
        !@_@ incongruence with Vector and std::vector
        *_* unnecesary state saving

*/

#include "all.h"
#include "itikBanks.h"
#include "newtonPoincare.h"
#include "fixPointsInSpace"


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
