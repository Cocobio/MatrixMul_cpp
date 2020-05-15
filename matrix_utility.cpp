#ifndef MATRIX_UTILITY
#define MATRIX_UTILITY

#include <ctime>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

// Gets a view for Matrices into a vector for rows
template <typename T>
void getViewOfMatrixWithPool(T** a, int r_begin, int r_end, int c_begin, T** pool_row_ref) {
	for (int i=r_begin; i<r_end; i++)
		pool_row_ref[i-r_begin] = &a[i][c_begin];
}

// Adds b to a, using r as rows and c as columns a = a+b
template <typename T>
void addMatrix(T** a, T** b, int r, int c) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			a[i][j] += b[i][j];
}

// Rests b on a, using r as rows and c as columns. a = a-b
template <typename T>
void restMatrix(T** a, T** b, int r, int c) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			a[i][j] -= b[i][j];
}

// Gets a view for matrices, and allocates de memory needed for the rows
template <typename T>
T** getViewOfMatrix(T** m, int r_begin, int r_end, int c_begin, int c_end) {
	T **view = (T**) malloc((r_end-r_begin)*sizeof(T*));

	for (int i=r_begin; i<r_end; i++)
		view[i-r_begin] = &m[i][c_begin];

	return view;
}

// Allocates the memory for a matrix
template <typename T>
T** allocateM(int r, int c) {
	T* m_pool = (T*) malloc(r*c*sizeof(T));
	T** m = (T**) malloc(r*sizeof(T*));

	for (int i=0; i<r; i++)
		m[i] = &m_pool[c*i];

	return m;
}

// Copies an entire matrix, using boundaries
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

// Adds two matrices with same size and returns a new matrix with the result
template <typename T>
T** matrixSum(T** a, T** b, int r, int c) {
	T** m = allocateM<T>(r, c);

	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m[i][j] = a[i][j]+b[i][j];

	return m;
}

// Rests two matrices with the same size and returns a new matrix with the result
template <typename T>
T** matrixRest(T** a, T** b, int r, int c) {
	T** m = allocateM<T>(r, c);

	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m[i][j] = a[i][j]-b[i][j];

	return m;
}

// Adds two matrices with same size and saves the result in c
template <typename T>
void matrixAddTo(T** x, T** y, T** z, int r, int c) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			z[i][j] = x[i][j]+y[i][j];
}

// Rests two matrices with same size and saves the result in C
template <typename T>
void matrixRestTo(T** x, T** y, T** z, int r, int c) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			z[i][j] = x[i][j]-y[i][j];
}

// Creates random values from offset to offset+range and populates a matrix with them.
template <typename T>
void randomizeM(T** m, int r, int c, float range, float offset) {
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m[i][j] = ((double) rand())/RAND_MAX*range+offset;
}

// Creates a transpose of a matrix.
template <typename T>
T** createTranspose(T** m, int r, int c) {
	T** m_T = allocateM<T>(c,r);

	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m_T[j][i] = m[i][j];

	return m_T;
}


// Binary file reader for matrix
template <typename T>
T** readMatrixFromFile(int *r, int *c, string file) {
	//we input the data from a binary file
	ifstream input;
	T **m;

	// It opens the file, in binary and reads the row and columns
	input.open(file, ios::binary | ios::in);
	input.read((char*)r, sizeof(int));
	input.read((char*)c, sizeof(int));
	
	// Allocates the needed memory
	m = allocateM<T>(*r, *c);

	// Reads the memory to the matrix
	input.read((char*)m[0], (*r)*(*c)*sizeof(T));
	input.close();

	return m;
}

// Binary file writer for matrix
template <typename T>
void writeFileWithMatrix(T** m, int r, int c, string file) {
	//we output the data into a binary file
	ofstream output;

	// It opens a binary file (if it exists it overwrites it)
	// It writes 2 ints with rows and columns for the matrix
	// Finally it writes the data
	output.open(file, ios::binary | ios::out | ios::trunc);
	output.write((char*)&r, sizeof(int));
	output.write((char*)&c, sizeof(int));
	output.write((char*)m[0], r*c*sizeof(T));
		
	// It closes the file
	output.close();
}

// Functions for measure performance of algoritms
template <typename T>
double measure_time_of(T** f (T**,int,int,T**,int,int), T **a, int r_a, int c_a, T **b, int r_b, int c_b, T*** c) {
	clock_t t; 
	t = clock(); 
	*c = f(a,r_a,c_a,b,r_b,c_b);
	t = clock() - t; 
	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
	return time_taken;
}

// Function for printing a matrix, for testing porpuses
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

#endif
