#ifndef MATRIX_UTILITY
#define MATRIX_UTILITY

#include <ctime>
#include <cstdlib>
#include <string>
#include <iostream>

using namespace std;

template <typename T>
void addMatrix(T** a, T** b, int r, int c) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			a[i][j] += b[i][j];
}

template <typename T>
void restMatrix(T** a, T** b, int r, int c) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			a[i][j] -= b[i][j];
}

template <typename T>
T** getViewOfMatrix(T** m, int r_begin, int r_end, int c_begin, int c_end) {
	T **view = (T**) malloc((r_end-r_begin)*sizeof(T*));

	for (int i=r_begin; i<r_end; i++)
		view[i-r_begin] = &m[i][c_begin];

	return view;
}

template <typename T>
T** allocateM(int r, int c) {
	T* m_pool = (T*) malloc(r*c*sizeof(T));
	T** m = (T**) malloc(r*sizeof(T*));

	for (int i=0; i<r; i++)
		m[i] = &m_pool[c*i];

	return m;
}

template <typename T>
T** copyMatrix(T** a, int r_begin, int r_end, int c_begin, int c_end) {
	int r = r_end-r_begin;
	int c = c_end-c_begin;

	T** m = allocateM<T>(r, c);
	
	for (int i=0; i<r; i++)
		// memcpy
		for (int j=0; j<c; j++)
			m[i][j] = a[r_begin+i][c_begin+j];

	return m;
}

template <typename T>
T** matrixSum(T** a, T** b, int r, int c) {
	T** m = allocateM<T>(r, c);

	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m[i][j] = a[i][j]+b[i][j];

	return m;
}

template <typename T>
T** matrixRest(T** a, T** b, int r, int c) {
	T** m = allocateM<T>(r, c);

	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m[i][j] = a[i][j]-b[i][j];

	return m;
}

template <typename T>
void matrixSumTo(T** a, T** b, int r, int c, T** m, int r_begin, int c_begin) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m[r_begin+i][c_begin+j] = a[i][j]+b[i][j];
}

template <typename T>
void matrixAddTo(T** x, T** y, T** z, int r, int c) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			z[i][j] = x[i][j]+y[i][j];
}

template <typename T>
void matrixRestTo(T** x, T** y, T** z, int r, int c) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			z[i][j] = x[i][j]-y[i][j];
}

template <typename T>
void matrixSumFromTo(T** a, int r_a_b, int c_a_b, T** b, int r_b_b, int c_b_b, int r, int c, T** m, int r_m_b, int c_m_b) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m[r_m_b+i][c_m_b+j] = a[r_a_b+i][c_a_b+j]+b[r_b_b+i][c_b_b+j];
}

template <typename T>
void matrixRestTo(T** a, T** b, int r, int c, T** m, int r_begin, int c_begin) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m[r_begin+i][c_begin+j] = a[i][j]-b[i][j];
}

template <typename T>
void matrixRestFromTo(T** a, int r_a_b, int c_a_b, T** b, int r_b_b, int c_b_b, int r, int c, T** m, int r_m_b, int c_m_b) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m[r_m_b+i][c_m_b+j] = a[r_a_b+i][c_a_b+j]-b[r_b_b+i][c_b_b+j];
}

template <typename T>
double measure_time_of(T** f (T**,int,int,T**,int,int), T **a, int r_a, int c_a, T **b, int r_b, int c_b, T*** c) {
	clock_t t; 
	t = clock(); 
	*c = f(a,r_a,c_a,b,r_b,c_b);
	t = clock() - t; 
	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
	return time_taken;
}

template <typename T>
double measure_time_of(void f (T**,int,int,T**,int,int,T**), T **a, int r_a, int c_a, T **b, int r_b, int c_b, T** c) {
	clock_t t; 
	t = clock(); 
	f(a,r_a,c_a,b,r_b,c_b,c);
	t = clock() - t; 
	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
	return time_taken;
}

// template <typename T>
// double measure_time_of(void f (T**,int,int,T**,int,int,T**,int,int), T **a, int r_a, int c_a, T **b, int r_b, int c_b, T** c, int r_c_b, int c_c_b) {
// 	clock_t t; 
// 	t = clock(); 
// 	f(a,r_a,c_a,b,r_b,c_b,c,r_c_b,c_c_b);
// 	t = clock() - t; 
// 	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
// 	return time_taken;
// }

template <typename T>
void printM(T** m, int r_m, int c_m, string name) {
	cout << name << ":";
	for (int i=0;i<r_m; i++) {
		cout << endl;
		for (int j=0; j<c_m; j++) {
			cout << " " << m[i][j];
		}
	}
	cout << endl;
}

template <typename T>
void randomizeM(T** m, int r, int c, float range, float offset) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m[i][j] = ((double) rand())/RAND_MAX*range+offset;
}

#endif
