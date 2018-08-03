#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

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
        if(++i >= b.size()) return; // prevents to go further of array size
    }
    b[i] += step;
}

template <class T>
void nextPointSubdomain(T &b, T xmin, T xmax)
{
    int i = 0;
    double step;
    while (equal(b[i], xmax[i])) // do it with equals
    {
        b[i] = xmin[i];
        if (++i >= b.size())
            return; // prevents to go further of array size
    }
    step = (xmax[i] - xmin[i]) / 2;
    b[i] += step;
}

template<class T>
void printVectorsReverse(const std::vector<T> xi)
{
    for (int i = xi.size() - 1; i >= 0; i--)
    {
        std::cout << xi[i] << "\t";
    }
    std::cout << "\n";
}

template <class T>
void printVector(std::vector<T> M)
{
    for (T j : M)
        std::cout << j << "\t";
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
    std::vector<double> xi={0,0.5,0};
    std::vector<double> xf={1,1.5,1};

    std::vector<double> point2(xi);
                                          //3 there's always three options
    for (i = 0; i < std::pow(point2.size(), 3); i++)
    {
        std::vector<double> reversePoint(point2.size());
        std::reverse_copy(point2.begin(),point2.end(), reversePoint.begin());
        printVector(reversePoint);

        nextPointSubdomain(point2, xi, xf);
    }
    std::cout << "Num subdomains : " << i << std::endl;
}
