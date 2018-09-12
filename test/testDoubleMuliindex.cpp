/**
 * Multindex double index by string contatenation 
 */

#include <cmath>
#include <iostream>

// 	return true;
// }
#include <map>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

std::string strtokey(double a, double b, int precision = 5){
	std::stringstream stream;
	stream << std::fixed << std::setprecision(precision) << a; 
	std::string as = stream.str();
    //https://stackoverflow.com/questions/20731/how-do-you-clear-a-stringstream-variable
    stream.str(std::string());
    stream << std::fixed << std::setprecision(precision) << b;
    std::string bs = stream.str();

	return as+bs;
}
// function key = strtokey(a, b)
// 	key = rounddecimals(hash(a, b), 4);
// key = num2str(key);
class FixPoint
{
public:
	bool convergent;
	bool stability;
	std::vector<double> solution;

	FixPoint(bool convergent_, bool stability_, std::vector<double> solution_) :
		convergent(convergent_),
		stability(stability_),
		solution(solution_) {}

	FixPoint(const FixPoint &point){
		convergent = point.convergent;
		stability = point.stability;
		solution = point.solution;
	}	
};

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
// With presition 5, they are the same
    std::cout<<strtokey(234.234213123,3.1412312)<<std::endl;  
    std::cout<<strtokey(234.2342131,3.1412312)<<std::endl;
//Not the same
    std::cout << strtokey(234.234213123, 3.1412312,8) << std::endl;
    std::cout << strtokey(234.2342131, 3.1412312,8) << std::endl;

    std::map<std::string, int> g;
    double x = 3.1416123254; 
    double y = x/10; 

    g.insert(std::pair<std::string,int> (strtokey(x,y),1 ) );
    g.insert(std::pair<std::string,int> (strtokey(x+2,y),2 ) );
    g.insert(std::pair<std::string,int> (strtokey(x+3,y),3 ) );
    g.insert(std::pair<std::string,int> (strtokey(x+4,y),4 ) );

    std::cout<<g[strtokey(x,y)]<<std::endl;
    std::cout<<g[strtokey(x+2,y)]<<std::endl;
    std::cout<<g[strtokey(x+3,y)]<<std::endl;
    std::cout<<g[strtokey(x+5,y)]<<std::endl;
    std::cout<<g[strtokey(x+4,10)]<<std::endl;
    std::cout<<g[strtokey(x+4,3.1416123254 / 10)]<<std::endl;
    //You can only store points to objects in a map:
    //https://stackoverflow.com/questions/2281420/c-inserting-a-class-into-a-map-container
    std::map<std::string,FixPoint*> S;
    FixPoint p(true,true, std::vector<double> {1.1,2.2,3.3});

    S.insert(std::make_pair<std::string,FixPoint*> (strtokey(x,y), &p));

    std::cout<<S[strtokey(x,y)]->solution[0]<<std::endl;

    if(1==1==1
    	==1==1)
    	std::cout<<"aa";
}  
