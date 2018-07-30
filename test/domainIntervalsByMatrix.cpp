#include <iostream>
#include <vector>
#include <bitset>
#include <cmath>
#define N 4
static int cprints = 0;

void shiftRowleft(std::vector<double> &xi, std::vector<double> &xm, std::vector<double> &xf, int row) {
	double temp1;
	temp1 = xi[row];
	xi[row] = xm[row];
	xm[row] = xf[row];
	xf[row] = temp1;

}

void shiftRowRight(std::vector<double> &xi, std::vector<double> &xm, std::vector<double> &xf, int row) {
	double temp1;
	temp1 = xf[row];
	xf[row] = xm[row];
	xm[row] = xi[row];
	xi[row] = temp1;

}

void printVectors(std::vector<double> &xi, std::vector<double> &xm, std::vector<double> &xf) {
	for (int i = 0; i < N; i++)
	{
		std::cout << xi[i] << " " << xm[i] << " " << xf[i] << "\n";
	}
	std::cout <<++cprints<< "--------------\n";
}

void printVectorsReverse(const std::vector<int> xi) {
	for (int i = xi.size() - 1; i >= 0; i--)
	{
		std::cout << xi[i] <<  " ";
	}
}
int nextBitSequence(std::vector<int> &b) {
	int n = b.size();
	int i = 0;
	while (b[i] == 1) {
		b[i] = 0;
		i++;
	}
	b[i] = 1;
	return i;
}

int main(){


	//std::vector<double> xiTemp;
	//std::vector<double> xmTemp;
	//std::vector<double> xfTemp;
	//printVectors(xi, xm, xf);
	//shiftRowleft(xi, xm, xf, 2);
	//std::cout << "\n";
	//printVectors(xi, xm, xf);

	//shiftRowleft(xi, xm, xf, 6);
	//std::cout << "\n";
	//printVectors(xi, xm, xf);

	//shiftRowleft(xi, xm, xf, 9);
	//shiftRowleft(xi, xm, xf, 9);
	//std::cout << "\n";
	//printVectors(xi, xm, xf);

	//xiTemp = xi;
	//xmTemp = xm;
	//xfTemp = xf;
	//printVectors(xiTemp, xmTemp, xfTemp);
	//for (int i = N - 1; i >= 0; i--) {
	//	int k = i;
	//		shiftRowleft(xiTemp, xmTemp, xfTemp, k);

	//		do {
	//			printVectors(xiTemp, xmTemp, xfTemp);
	//			if (k >= N - 1) break;
	//			shiftRowleft(xiTemp, xmTemp, xfTemp, k + 1);
	//			k++;
	//		} while (k <= N);


	//		xiTemp = xi;
	//		xmTemp = xm;
	//		xfTemp = xf;
	//}

	//std::string s = std::bitset< 64 >(12345).to_string(); // string conversion
	//std::cout << s[0];
	//std::cout << std::bitset< 8 >(255) << ' '; // direct output
	
	//std::bitset< 64 > input;
	//std::cin >> input;
	//unsigned long ul = input.to_ulong();


	//for (int i = 0; i < pow(2,N) -1 ; i++)
	//{
	//	int j;
	//	printVectors(xi, xm, xf);
	//	j = nextBitSequence(bitArray);
	//	shiftRowleft(xi, xm, xf, j);

	//}
	std::vector<double> xi(N, 1);
	std::vector<double> xm(N, 2);
	std::vector<double> xf(N, 3);
	std::vector<int>  bitArray(N, 0);

	for (int i = 0; i < pow(2, N) - 1 ; i++) {
		printVectors(xi, xm, xf);
		int j = 0;
		while (bitArray[j] == 1) {
			bitArray[j] = 0;
			shiftRowRight(xi, xm, xf, j);
			j++;
		}
		bitArray[j] = 1;
		shiftRowleft(xi, xm, xf, j);
	}
	printVectors(xi, xm, xf);

	std::cin.get();
}