// #include <cstdlib>
#include <string>
#include <cstring>
#include "matrix_utility.cpp"

template <typename T>
void getViewOfMatrixWithPool(T** a, int r_begin, int r_end, int c_begin, T** pool_row_ref) {
	for (int i=r_begin; i<r_end; i++)
		pool_row_ref[i-r_begin] = &a[i][c_begin];
}

template <typename T>
void _matrixMultiplication_strassen(T** a, int r_a, int c_a, T** b, int r_b, int c_b, T** c, T* pool, T** pool_row_ref) {
	if (c_a!=r_b) { 
		std::cout << "error, matrix dont align" << std::endl;
		return;
	}

	// Base case return normal multiplication
	if (r_a<=1 || c_a<=1 || r_b<=1 || c_b<=1) {
		for (int i=0; i<r_a; i++)
			memset(c[i], 0, c_b*sizeof(T));

		for (int i=0; i<r_a; i++)
			for (int k=0; k<r_b; k++)
				for (int j=0; j<c_b; j++)
					c[i][j] += a[i][k]*b[k][j];
		return;
	}

	int half_ra = r_a/2;
	int half_ca = c_a/2;
	int half_rb = r_b/2;
	int half_cb = c_b/2;

	int r_c = r_a;
	int c_c = c_b;

	int half_rc = r_c/2;
	int half_cc = c_c/2;

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


	// P6 to I
	matrixRestTo(B, D, II,  half_ra, half_ca);
	matrixAddTo (G, H, III, half_rb, half_cb);
	_matrixMultiplication_strassen(II, half_ra, half_ca, III, half_rb, half_cb, I, pool, pool_row_ref);		// P6

	// P5 to IV
	matrixAddTo(A, D, II,  half_ra, half_ca);
	matrixAddTo(E, H, III, half_rb, half_cb);
	_matrixMultiplication_strassen(II, half_ra, half_ca, III, half_rb, half_cb, IV, pool, pool_row_ref);		// P5

	// I = P6+P5
	addMatrix(I, IV, half_rc, half_cc);

	// P7 to tmp_1
	matrixRestTo(A, C, II,  half_ra,half_ca);
	matrixAddTo (E, F, III, half_rb,half_cb);
	_matrixMultiplication_strassen(II, half_ra, half_ca, III, half_rb, half_cb, tmp_1, pool, pool_row_ref);	// P7

	// IV = P5-P7
	restMatrix(IV, tmp_1, half_rc, half_cc);

	// P1 to II
	matrixRestTo(F, H, III, half_rb, half_cb);
	_matrixMultiplication_strassen(A, half_ra, half_ca, III, half_rb, half_cb, II, pool, pool_row_ref);		// P1

	// IV = P5-P7+P1
	addMatrix(IV, II, half_rc, half_cc);

	// P3 to III
	matrixAddTo(C, D, tmp_1, half_ra, half_ca);
	_matrixMultiplication_strassen(tmp_1, half_ra, half_ca, E, half_rb, half_cb, III, pool, pool_row_ref);		// P3

	// IV = P5-P7+P1-P3
	restMatrix(IV, III, half_rc, half_cc);

	// P2 to tmp_2
	matrixAddTo(A, B, tmp_1, half_ra, half_ca);
	_matrixMultiplication_strassen(tmp_1, half_ra, half_ca, H, half_rb, half_cb, tmp_2, pool, pool_row_ref);	//P_2

	// I = P6+P5-P2
	restMatrix(I, tmp_2, half_rc, half_cc);

	// II = P1+P2
	addMatrix(II, tmp_2, half_rc, half_cc);

	// P4 to tmp_2
	matrixRestTo(G, E, tmp_1, half_rb, half_cb);
	_matrixMultiplication_strassen(D, half_ra, half_ca, tmp_1, half_rb, half_cb, tmp_2, pool, pool_row_ref);	//P_4

	// III = P3+P4
	addMatrix(III, tmp_2, half_rc, half_cc);

	// I = P6+P5-P2+P4
	addMatrix(I, tmp_2, half_rc, half_cc);
}

template <typename T>
T** matrixMultiplication_strassen(T** a, int r_a, int c_a, T** b, int r_b, int c_b) {
	T** c = allocateM<T>(r_a, c_b);

	T* pool = (T*)malloc(r_a*c_b*sizeof(T)); 
	T** pool_row_ref = (T**) malloc(14*r_a*sizeof(T*));

	_matrixMultiplication_strassen(a, r_a, c_a, b, r_b, c_b, c, pool, pool_row_ref);

	free(pool);
	free(pool_row_ref);

	return c;
}


#ifndef MATRIX_TEST_CPP
// Testing functions
#include <iostream>

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
	// for (int i=0; i<r_b*c_b; i++)
	// 	b[i/c_b][i%c_b] = i;
	// printM(b, r_b, c_b, "b");

	float **c; // = matrixMultiplication_strassen(a, r_a, c_a, b, r_b, c_b);
	double t = measure_time_of(matrixMultiplication_strassen, a, r_a, c_a, b, r_b, c_b, &c);
	cout << "matrix multiplication strassen in " << t << " sec " << endl;
	
	// if (c!= NULL) printM(c, r_a, c_b, "c");
}

#endif

/*
strassen:
512 = 19.069 sec	// 8.831 sec		// 2.686 sec
1024 = 76.821 sec 	// 62.172 sec		// 18.714 sec 
2048 = 132.02 sec

Transpose:
512 = 0.59 sec
1024 = 4.762 sec
2048 = 36.252 sec

Naive:
512 = 1.01 sec 
1024 = 18.557 sec
2048 = 206.439 sec

*/