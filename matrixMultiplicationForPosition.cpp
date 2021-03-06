/*
Matrix multiplication algorithm. Naive version.

time complexity = O(n^3)
space complexity = O(1)

*** A bunch faster, because the inner most loop iterates through columns, instead of rows. (less cache misses)
*/

#include <string>
#include <cstring>
#include "matrix_utility.cpp"


template <typename T>
T** matrixMultiplicationForPosition(T** a, int r_a, int c_a, T** b, int r_b, int c_b) {
	if (c_a!=r_b) return NULL;
	int r_c = r_a;
	int c_c = c_b;

	T** c = allocateM<T>(r_c,c_c);
	memset(c[0], 0, r_c*c_c*sizeof(T));

	for (int i=0; i<r_a; i++)
		for (int k=0; k<r_b; k++)
			for (int j=0; j<c_b; j++)
				c[i][j] += a[i][k]*b[k][j];

	return c;
}


#ifndef MATRIX_TEST_CPP
// Testing functions
#include <iostream>
#include <ctime>

using namespace std;

int main() {
	// srand(time(0));
	int squared = 512;
	int r_a = squared;
	int c_b = squared;
	int n = squared;
	
	int c_a = n;
	int r_b = n;

	float **a = allocateM<float>(r_a, c_a);

	for (int i=0; i<r_a*c_a; i++)
		a[i/c_a][i%c_a] = i;

	// printM(a, r_a, c_a, "a");

	float **b = allocateM<float>(r_b, c_b);
	randomizeM(b, r_b, c_b, 4.0, 1.0);
	// printM(b, r_b, c_b, "b");

	float **c;
	double t = measure_time_of(matrixMultiplicationForPosition, a, r_a, c_a, b, r_b, c_b, &c);
	// cout << "matrix multiplication in " << t << " sec " << endl;

	// printM(c, r_a, c_b, "c");
}

#endif

