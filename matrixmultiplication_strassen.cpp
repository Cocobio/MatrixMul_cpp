/*
Strassen's matrix multiplication algorithm. It uses divide and conquer paradigm.
It divides the data into 7 matrix product of half size, with a cuadratic work for add and rest.
T(n) = 7*T(n/2) + O(n^2)

a = [ A,  B,
	  C,  D ]

b = [ E,  F,
	  G,  H ]

P1 = A(F-H)
P2 = (A+B)H
P3 = (C+D)E
P4 = D(G-E)
P5 = (A+D)(E+H)
P6 = (B-D)(G+H)
P7 = (A-C)(E+F)

c = [ P5+P4-P2+P6,	P1+P2,
	  P3+P4,		P1+P5-P3-P7 ]

time complexity = O(n^2.8074)
space complexity = O(n^2)       ### temporal matrices for the 7 products, and the sums of the different divisions

*** This is a purely recursive implementation. By creating a^(log(n)/log(2)) instances, the complexity reduction gain is mitigated by the recursive constant.
*** It squares the input to the closes power of 2 size. ------> This makes the algorithm slower.
*/

#include <string>
#include <cstring>
#include "matrix_utility.cpp"
#include <cmath>

// Strassen purely recursive implementation
template <typename T>
void _matrixMultiplicationStrassen(T** a, int r_a, int c_a, T** b, int r_b, int c_b, T** c, T* pool, T** pool_row_ref) {
	// Base case return normal multiplication
	if (r_a<=1) { // || c_a<=1 || r_b<=1 || c_b<=1) {
		for (int i=0; i<r_a; i++)
			memset(c[i], 0, c_b*sizeof(T));

		for (int i=0; i<r_a; i++)
			for (int k=0; k<r_b; k++)
				for (int j=0; j<c_b; j++)
					c[i][j] += a[i][k]*b[k][j];
		return;
	}

	// For future non squared matrices
	int half_ra = r_a/2;
	int half_ca = c_a/2;
	int half_rb = r_b/2;
	int half_cb = c_b/2;

	int r_c = r_a;
	int c_c = c_b;

	int half_rc = r_c/2;
	int half_cc = c_c/2;

	// It creates references into the sub-matrices of a and b
	T** A = pool_row_ref;
	getViewOfMatrixWithPool(a, 0, half_ra, 0, pool_row_ref);
	pool_row_ref += half_ra;
	T** B = pool_row_ref;
	getViewOfMatrixWithPool(a, 0, half_ra, half_ca, pool_row_ref);
	pool_row_ref += half_ra;
	T** C = pool_row_ref;
	getViewOfMatrixWithPool(a, half_ra,r_a,	0, pool_row_ref);
	pool_row_ref += half_ra;
	T** D = pool_row_ref;
	getViewOfMatrixWithPool(a, half_ra,r_a,	half_ca, pool_row_ref);
	pool_row_ref += half_ra;

	T** E = pool_row_ref;
	getViewOfMatrixWithPool(b,0,half_rb,0,pool_row_ref);
	pool_row_ref += half_rb;
	T** F = pool_row_ref;
	getViewOfMatrixWithPool(b,0,half_rb,half_cb,pool_row_ref);
	pool_row_ref += half_rb;
	T** G = pool_row_ref;
	getViewOfMatrixWithPool(b,half_rb,r_b,0,pool_row_ref);
	pool_row_ref += half_rb;
	T** H = pool_row_ref;
	getViewOfMatrixWithPool(b,half_rb,r_b,half_cb,pool_row_ref);
	pool_row_ref += half_rb;


	// It creates references into the sub-matrices for result
	T** I   = pool_row_ref;
	getViewOfMatrixWithPool(c, 0, half_rc, 0, pool_row_ref);
	pool_row_ref += half_rc;
	T** II  = pool_row_ref;
	getViewOfMatrixWithPool(c, 0, half_rc, half_cc, pool_row_ref);
	pool_row_ref += half_rc;
	T** III = pool_row_ref;
	getViewOfMatrixWithPool(c, half_rc, r_c, 0, pool_row_ref);
	pool_row_ref += half_rc;
	T** IV  = pool_row_ref;
	getViewOfMatrixWithPool(c, half_rc, r_c, half_cc, pool_row_ref);
	pool_row_ref += half_rc;

	// Get view on pool data
	T** tmp_1 = pool_row_ref;
	for (int i=0; i<half_rc; i++)
		tmp_1[i] = &pool[half_cc*i];
	pool_row_ref += half_rc;

	T** tmp_2 = pool_row_ref;
	for (int i=0; i<half_rc; i++)
		tmp_2[i] = &pool[half_cc*half_rc + half_cc*i];
	pool_row_ref += half_rc;

	pool += r_c*c_c/2;


	// We calculate the different P1-P7, using the result matrix as temporal sub-matrices. So we use less memory.
	// P6 to I
	matrixRestTo(B, D, II,  half_ra, half_ca);
	matrixAddTo (G, H, III, half_rb, half_cb);
	_matrixMultiplicationStrassen(II, half_ra, half_ca, III, half_rb, half_cb, I, pool, pool_row_ref);		// P6

	// P5 to IV
	matrixAddTo(A, D, II,  half_ra, half_ca);
	matrixAddTo(E, H, III, half_rb, half_cb);
	_matrixMultiplicationStrassen(II, half_ra, half_ca, III, half_rb, half_cb, IV, pool, pool_row_ref);		// P5

	// I = P6+P5
	addMatrix(I, IV, half_rc, half_cc);

	// P7 to tmp_1
	matrixRestTo(A, C, II,  half_ra,half_ca);
	matrixAddTo (E, F, III, half_rb,half_cb);
	_matrixMultiplicationStrassen(II, half_ra, half_ca, III, half_rb, half_cb, tmp_1, pool, pool_row_ref);	// P7

	// IV = P5-P7
	restMatrix(IV, tmp_1, half_rc, half_cc);

	// P1 to II
	matrixRestTo(F, H, III, half_rb, half_cb);
	_matrixMultiplicationStrassen(A, half_ra, half_ca, III, half_rb, half_cb, II, pool, pool_row_ref);		// P1

	// IV = P5-P7+P1
	addMatrix(IV, II, half_rc, half_cc);

	// P3 to III
	matrixAddTo(C, D, tmp_1, half_ra, half_ca);
	_matrixMultiplicationStrassen(tmp_1, half_ra, half_ca, E, half_rb, half_cb, III, pool, pool_row_ref);		// P3

	// IV = P5-P7+P1-P3
	restMatrix(IV, III, half_rc, half_cc);

	// P2 to tmp_2
	matrixAddTo(A, B, tmp_1, half_ra, half_ca);
	_matrixMultiplicationStrassen(tmp_1, half_ra, half_ca, H, half_rb, half_cb, tmp_2, pool, pool_row_ref);	//P_2

	// I = P6+P5-P2
	restMatrix(I, tmp_2, half_rc, half_cc);

	// II = P1+P2
	addMatrix(II, tmp_2, half_rc, half_cc);

	// P4 to tmp_2
	matrixRestTo(G, E, tmp_1, half_rb, half_cb);
	_matrixMultiplicationStrassen(D, half_ra, half_ca, tmp_1, half_rb, half_cb, tmp_2, pool, pool_row_ref);	//P_4

	// III = P3+P4
	addMatrix(III, tmp_2, half_rc, half_cc);

	// I = P6+P5-P2+P4
	addMatrix(I, tmp_2, half_rc, half_cc);
}

// Drive function, for allocation memory of temporal matrices.
// It also squares the input matrices to the closes power of 2 number.
template <typename T>
T** matrixMultiplicationStrassen(T** a, int r_a, int c_a, T** b, int r_b, int c_b) {
	// Returns NULL if the matrices can't be multiplied
	if (c_a!=r_b) { 
		std::cout << "error, matrices don't align" << std::endl;
		return NULL;
	}

	// It checks wether the input matrices are squared or not
	// creates new variables, that will be pass as parameters for the algorithm.
	bool squared=true;
	T** a_tmp, **b_tmp, **c_tmp;
	int r_a_tmp, c_a_tmp, r_b_tmp, c_b_tmp;
	r_a_tmp = r_a;
	
	if (r_a!=c_a || r_b!=c_b) {
		squared = false;

		if (r_a_tmp<c_a) r_a_tmp = c_a;
		if (r_a_tmp<c_b) r_a_tmp = c_b;
	} 

	c_a_tmp = r_a_tmp;
	r_b_tmp = r_a_tmp;
	c_b_tmp = r_a_tmp;

	// Checks wheter the size of the matrix is a power of 2
	bool is_a_power_of_2 = ceil(log(r_a_tmp)/log(2)) == log(r_a_tmp)/log(2);


	// Whether is not squared or not a power of 2
	// make a squared with a size that is power of 2
	if (!squared || !is_a_power_of_2) {
		int exp = ceil(log(r_a_tmp)/log(2));
		r_a_tmp = std::pow(2, exp);
		c_a_tmp = r_a_tmp;
		r_b_tmp = r_a_tmp;
		c_b_tmp = r_a_tmp;

		a_tmp = allocateM<T>(r_a_tmp, c_a_tmp);
		b_tmp = allocateM<T>(r_b_tmp, c_b_tmp);

		memset(a_tmp[0], 0, r_a_tmp*c_a_tmp*sizeof(T));
		memset(b_tmp[0], 0, r_b_tmp*c_b_tmp*sizeof(T));

		for (int i=0; i<r_a; i++)
			for (int j=0; j<c_a; j++)
				a_tmp[i][j] = a[i][j];

		for (int i=0; i<r_b; i++)
			for (int j=0; j<c_b; j++)
				b_tmp[i][j] = b[i][j];
	}
	else {
		a_tmp = a;
		b_tmp = b;
	}

	// allocate result memory
	c_tmp = allocateM<T>(r_a_tmp, c_b_tmp);

	// allocate pool memory for temporal files inside Strassen
	T* pool = (T*)malloc(r_a_tmp*c_b_tmp*sizeof(T)); 
	T** pool_row_ref = (T**) malloc(14*r_a_tmp*sizeof(T*));

	// STRASSEN
	_matrixMultiplicationStrassen(a_tmp, r_a_tmp, c_a_tmp, b_tmp, r_b_tmp, c_b_tmp, c_tmp, pool, pool_row_ref);

	// Free the memory asked for tmp files
	free(pool);
	free(pool_row_ref);

	// Return the result, but if we created a squared matrix, we need to release the memory.
	if (!squared || !is_a_power_of_2) {
		T** c = allocateM<T>(r_a, c_b);

		for (int i=0; i<r_a; i++)
			for (int j=0; j<c_b; j++)
				c[i][j] = c_tmp[i][j];

		free(a_tmp[0]);
		free(b_tmp[0]);
		free(c_tmp[0]);
		free(a_tmp);
		free(b_tmp);
		free(c_tmp);

		return c;
	}
	else return c_tmp;
}


#ifndef MATRIX_TEST_CPP
// Testing functions
#include <iostream>

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
	// for (int i=0; i<r_b*c_b; i++)
	// 	b[i/c_b][i%c_b] = i;
	// printM(b, r_b, c_b, "b");

	float **c; // = matrixMultiplicationStrassen(a, r_a, c_a, b, r_b, c_b);
	double t = measure_time_of(matrixMultiplicationStrassen, a, r_a, c_a, b, r_b, c_b, &c);
	cout << "matrix multiplication strassen in " << t << " sec " << endl;
	
	// if (c!= NULL) printM(c, r_a, c_b, "c");
}

#endif

