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
Note:
 The algorithm finds the same points several times. How this can be avoided?
 -Why to disthinguis beetween std::vector and boost::vector how about just use boost vector
  for everything

*/

#include "all.h"
#include "itikBanks.h"
#include "newtonPoincare.h"
#include "fixPointsInSpace.h"


//typedef std::vector<fixPoint> fixPoints;
int main(){

	std::vector<fixPoint> S;
		itikBanks fun(PARAMETERS);
		S = al21(0.1, 0.1, 0.1, 2, 2, 2, fun, 10, 0.5);

	for (fixPoint i : S) {
		std::cout << i.solution[0] << " " << i.solution[1] << " " << i.solution[2]
				  << " " << i.stability<<std::endl;
	}
	std::cout<<"program finished!"<< std::endl;
	std::cin.get();//VS windows only
	return 0;
}
