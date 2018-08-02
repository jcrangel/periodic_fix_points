/**
g++ permutations.cpp  -std=c++11
*/
#include <algorithm>
#include <vector>
#include <iostream>
#include <utility>
 
template <class T>
void printArr(T a[],int n)
{
    for (int i=0; i<n; i++)
        std::cout << a[i] << " ";
    std::cout << "\n" ;
}
/**
Create the cartesian set A x B
*/
template<class T>
std::vector< std::vector<T> > cartesianProduct(const std::vector<T> A,const std::vector<T> B)
{
	std::vector<std::vector<T> > prod;
	for(int i = 0; i < A.size(); i ++)
		for(int j = 0; j < B.size(); j++)
				prod.push_back(std::vector<T> {A[i],B[j]});

	return prod;		
}
/**
Create the cartesian set A x B = (a x b) x B, where A is already a cartesian set.
Then to create the cartesian set of {x,y,z} x {x,y,z} x {x,y,z}. We have to :
 s = {x,y,z};	
 std::vector<std::vector<double> > prod = cartesianProduct(s,s);
 prod = cartesianProduct(prod,s);
*/
template<class T>
std::vector< std::vector<T> > cartesianProduct(const std::vector<std::vector<T>> A,const std::vector<T> B)
{
	std::vector<std::vector<T> > prod;
	for(int i = 0; i < A.size(); i ++)
		for(int j = 0; j < B.size(); j++){
				std::vector<T> AB;
				AB.reserve(A.size()+B.size());
				AB.insert(AB.end(), A[i].begin(),A[i].end());
				AB.insert(AB.end(), B[j]);

				prod.push_back(AB);

			}
	return prod;		
}

template<class T>
void printVectorVector(std::vector< std::vector<T> > M)
{

	for(std::vector<T> i: M){
		for(T j: i)
			std::cout<< j <<" ";
	std::cout<<"\n";
	}

}
// template<class T>
// void swap(T &a, T &b)
// {
// T temp;
// temp = a;
// a=b;
// b=temp;
// } 

// Generating permutation using Heap Algorithm
template <class T>
void heapPermutation(T a[], int size, int n)
{
    // if size becomes 1 then prints the obtained
    // permutation
    if (size == 1)
    {
        printArr(a, n);
        return;
    }
 
    for (int i=0; i<size; i++)
    {
        heapPermutation(a,size-1,n);
 
        // if size is odd, swap first and last
        // element
        if (size%2==1)
            std::swap(a[0], a[size-1]);
 
        // If size is even, swap ith and last
        // element
        else
            std::swap(a[i], a[size-1]);
    }
}



int main()
{
    std::vector<double> s;
    s.push_back(1);
    s.push_back(2);
    s.push_back(3);
    std::sort(s.begin(), s.end());
    do {
        std::cout << s[0] <<","<< s[1] <<"," << s[2] << '\n';
    } while(std::next_permutation(s.begin(), s.end()));

    
    double a[] = {1, 0.5, 1};
    // int n = sizeof a/sizeof a[0];
    std::cout << "\n";
    heapPermutation(a, 2, 3);
    std::cout << "\n";


    std::vector<std::vector<double> > prod = cartesianProduct(s,s);
    for(int i = 1 ; i <= 4 ; i++ ){
    	prod = cartesianProduct(prod,s);
    }
    printVectorVector(prod);
    std::cout << prod.size() <<std::endl;
}