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

  Why is not working?
  -- DFode and DF, DF calculates the jacobian only using the last data. I can't find any reference were
		this should be done like that
	--%algorithm 2.2 line 88 step1




*/
#include "al21Manager.h"

//#include "all.h"
//typedef std::vector<fixPoint> fixPoints;
//extern std::vector<double> PARAMETERS;// TODO , delete this, PARAMETERS should be global
int main(){


	runAl21();

	return 0;
}
