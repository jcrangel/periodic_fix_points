#include <iostream>
#include <cmath>
using namespace std;

inline double round( double val )
{
    if( val < 0 ) return ceil(val - 0.5);
    return floor(val + 0.5);
}

int main(){
	double a = 0.0000123123123129731892;

	cout<<round(a);
}

