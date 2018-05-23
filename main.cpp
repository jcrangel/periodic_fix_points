#include "all.h"


int main(){

	stateType IC = {0,1,0.06};
	std::vector<double> param={1,2.5,0.5,1.5,4.5,1,0.2,0.5};
	itikBanks fun(param);
    matrix<double> m = DFode(fun,IC,10,0.5);
    std::cout<< m <<std::endl;
	std::cin.get();//VS windows only
}