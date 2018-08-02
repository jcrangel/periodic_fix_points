#include <iostream>
#include <vector>
#include <cmath>


bool equal(double A, double B, double epsilon = 0.000005f)
{
    return (std::fabs(A - B) < epsilon);
}

void nextPointSubdomain(std::vector<double> &b, double min, double max, double step)
{
    int i = 0;

    while (equal(b[i], max - step) ) // do it with equals
    {
        b[i] = min;
        if(++i > b.size()) return;
    }
    b[i] += step;
}

template<class T>
void printVectorsReverse(const std::vector<T> xi)
{
    for (int i = xi.size() - 1; i >= 0; i--)
    {
        std::cout << xi[i] << " ";
    }
    std::cout << "\n";
}

template <class T>
void printVector(std::vector<T> M)
{
    for (T j : M)
        std::cout << j << " ";
    std::cout << "\n";
}

template <class T>
void printVectorInterval(std::vector<T> M,T step)
{
    for (T j : M)
        std::cout << "["<<j << ","<< j + step <<"]   ";
    std::cout << "\n";
}
int main()
{
    double N = 3;
    int divisions = 10;
    double min = 0, max = 1;
    double step = (max - min) / divisions;


    std::vector<double> point(N,min);


    int i;
    for( i = 0;i < std::pow(divisions, N) ; i++ ){
        printVectorInterval(point,step);
        nextPointSubdomain(point, min, max, step);
    }
    std::cout<< "Num subdomains : " <<divisions << "^" <<N << "= " << i <<std::endl;
}
