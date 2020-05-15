#define MATRIX_TEST_CPP

#include <iostream>
#include <string>
#include <vector>
#include "matrix_utility.cpp"

#include "matrixmultiplication_naive.cpp"
#include "matrixMultiplicationForPosition.cpp"
#include "matrixmultiplication_transpose.cpp"
#include "matrixmultiplication_strassen.cpp"
#include "matrixmultiplication_strassen_optimization.cpp"

using namespace std;

typedef float data_type;

int main() {
	//// List of algorithms
	vector<data_type** (*)(data_type**,int,int,data_type**,int,int)> algorithms;
	vector<string> f_names;

	// uncomment for sorting integer vectors created by create_data_set
	algorithms.push_back(matrixMultiplication);
	f_names.push_back("matrixMultiplication");
	
	algorithms.push_back(matrixMultiplicationForPosition);
	f_names.push_back("matrixMultiplicationForPosition");

	algorithms.push_back(matrixMultiplicationStrassen);
	f_names.push_back("matrixMultiplicationStrassen");
	
	algorithms.push_back(matrixMultiplicationStrassen_Opti);
	f_names.push_back("matrixMultiplicationStrassen_Opti");

	algorithms.push_back(matrixMultiplicationTranspose);
	f_names.push_back("matrixMultiplicationTranspose");

	// input and output data files
	string input_a, input_b, output;

	cout << "This programs will assume float matrices." << endl;

	cout << "Input file Matrix A: ";
	cin >> input_a;
	cout << "Input file Matrix B: ";
	cin >> input_b;

	cout << "Output file: ";
	cin >> output;

	int option;

	cout << "Algorithm list:" << endl;
	for (int i=0; i<f_names.size(); i++) {
		cout << i << ":\t" << f_names[i] << endl;
	}

	cout << "Select Algorithm: ";
	cin >> option;

	if (option<0 || option>=f_names.size()) {
		cout << "Option out of bound." << endl;
		return -1;
	}

	data_type **a, **b, **c;
	int r_a, c_a, r_b, c_b;

	cout << "Reading file: " << input_a << endl;	
	a = readMatrixFromFile<data_type>(&r_a, &c_a, input_a);
	cout << "Reading file: " << input_b << endl;
	b = readMatrixFromFile<data_type>(&r_b, &c_b, input_b);

	
	double t = measure_time_of(algorithms[option], a, r_a, c_a, b, r_b, c_b, &c);
	cout << f_names[option] << " in " << t << "sec for " << r_a << "x" << c_a << " elements." << endl;
	writeFileWithMatrix(c, r_a, c_b, output);

	return 0;
}

		