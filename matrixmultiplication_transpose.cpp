// #include <cstdlib>
#include <string>
#include <cstring>
#include "matrix_utility.cpp"

template <typename T>
T** createTranspose(T** m, int r, int c) {
	T** m_T = allocateM<T>(c,r);

	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m_T[j][i] = m[i][j];

	return m_T;
}

template <typename T>
T** matrixMultiplication_transpose(T** a, int r_a, int c_a, T** b, int r_b, int c_b) {
	if (c_a!=r_b) return NULL;
	int r_c = r_a;
	int c_c = c_b;

	T** c = allocateM<T>(r_c, c_c);
	memset(c[0], 0, r_c*c_c*sizeof(T));
	
	T** b_T = createTranspose(b, r_b, c_b);

		for (int j=0; j<c_b; j++)
	for (int i=0; i<r_a; i++)
			for (int k=0; k<r_b; k++)
				c[i][j] += a[i][k]*b_T[j][k];

	free(b_T[0]);
	free(b_T);

	return c;
}


#ifndef MATRIX_TEST_CPP
// Testing functions
#include <iostream>
#include <ctime>

using namespace std;

int main() {
	// srand(time(0));
	int squared = 1024;
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
	double t = measure_time_of(matrixMultiplication_transpose, a, r_a, c_a, b, r_b, c_b, &c);
	cout << "matrix multiplication transpose in " << t << " sec " << endl;
	
	// if (c!= NULL) printM(c, r_a, c_b, "c");
}

#endif

 
 