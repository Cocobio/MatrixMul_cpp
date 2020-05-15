/*
File for data creation, and time measurement on the algorithms
*/

#define MATRIX_TEST_CPP

#include <iostream>
#include <string>
#include <vector>

#include <sstream>
#include <iomanip>

#include "matrix_utility.cpp"

#include "matrixmultiplication_naive.cpp"
#include "matrixMultiplicationForPosition.cpp"
#include "matrixmultiplication_transpose.cpp"
#include "matrixmultiplication_strassen.cpp"
#include "matrixmultiplication_strassen_optimization.cpp"

#include <Windows.h>

using namespace std;

// Inspecting folder for a list of all files
vector<string> getFileListWithinFolder(std::string folder) {
    vector<string> names;
    string search_path = folder + "/*.*";
    WIN32_FIND_DATA fd; 
    HANDLE hFind = ::FindFirstFile(search_path.c_str(), &fd); 
    if(hFind != INVALID_HANDLE_VALUE) { 
        do { 
            // read all (real) files in current folder
            // , delete '!' read other 2 default folder . and ..
            if(! (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) ) {
                names.push_back(folder+"/"+fd.cFileName);
            }
        }while(::FindNextFile(hFind, &fd)); 
        ::FindClose(hFind); 
    }
    return names;
}

// Data creation
void createSquaredData() {
	int n = 3000;
	int step = 100;
	string name;
	string folder;
	float **m;

	folder = "A/";
	name = "mat_";

	for (int row=step; row<n; row+=step) {
		stringstream ss;
		ss << std::setw(4) << std::setfill('0') << row;

		m = allocateM<float>(row, row);
		randomizeM(m, row, row, 4.0, 1.0);

		writeFileWithMatrix(m, row, row, folder+name+ss.str()+".bin");

		free(m[0]);
		free(m);
	}

	folder = "B/";

	for (int row=step; row<n; row+=step) {
		stringstream ss;
		ss << std::setw(4) << std::setfill('0') << row;

		m = allocateM<float>(row, row);
		randomizeM(m, row, row, 4.0, 1.0);

		writeFileWithMatrix(m, row, row, folder+name+ss.str()+".bin");

		free(m[0]);
		free(m);
	}
}

// Data creation
void createNonSquaredData() {
	int n = 1000;
	int fixed = 250;
	int step = 100;
	string name;
	string folder;
	float **m;

	folder = "A/nonSquared/";
	name = "mat_";

	for (int row=step; row<n; row+=step) {
		stringstream ss;
		ss << std::setw(4) << std::setfill('0') << row;

		m = allocateM<float>(row, fixed);
		randomizeM(m, row, fixed, 4.0, 1.0);

		writeFileWithMatrix(m, row, fixed, folder+name+ss.str()+".bin");

		free(m[0]);
		free(m);
	}

	folder = "B/nonSquared/";

	for (int column=step; column<n; column+=step) {
		stringstream ss;
		ss << std::setw(4) << std::setfill('0') << column;

		m = allocateM<float>(fixed, column);
		randomizeM(m, fixed, column, 4.0, 1.0);

		writeFileWithMatrix(m, fixed, column, folder+name+ss.str()+".bin");

		free(m[0]);
		free(m);
	}
}

// Testing an algorithms
template <typename T>
void runTMatrixTestsWriteDataTo(T** (*s)(T**,int,int,T**,int,int), string folder_a, string folder_b,string output_file) {
	// We define the matrix variables
	T **a, **b, **c;
	int r_a,r_b,c_a,c_b;
	double t;
	int n = 10;
	float t_avr;	

	// we search the file list for the datasets and iterate over them
	vector<string> files_a = getFileListWithinFolder(folder_a);
	vector<string> files_b = getFileListWithinFolder(folder_b);
	cout.precision(6);

	ofstream output;

	output.open(output_file, ios::out | ios::trunc);

	for (int i=0; i<files_a.size(); i++) {
		t_avr = 0;
		
		// We read the files as matrices
		cout << "Reading file: " << files_a[i] << endl;
		a = readMatrixFromFile<float>(&r_a, &c_a, files_a[i]);
		cout << "Reading file: " << files_b[i] << endl;
		b = readMatrixFromFile<float>(&r_b, &c_b, files_b[i]);

		for (int j=0; j<n; j++) {
			// we measure the time
			t = measure_time_of(s, a, r_a, c_a, b, r_b, c_b, &c);
			t_avr += t/n;

			free(c[0]);
			free(c);
		}
		output << t_avr << ";";

		// we free de memory that the vector and temporal vector were using
		free(a[0]);
		free(b[0]);
		free(a);
		free(b);
	}

	// we close the file
	output.close();
}

// Test all algorithms
void runTest(string folder_a, string folder_b, string folder_c) {
	vector<float** (*)(float**,int,int,float**,int,int)> algorithms;
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

	for (int i=0; i<algorithms.size(); i++) {
		cout << "Testing with " << f_names[i] << endl;
		runTMatrixTestsWriteDataTo<float>(algorithms[i], folder_a, folder_b, folder_c+f_names[i]+".txt");
	}
}

// Test only measures time. It doesn't write the output matrix to a file.
int main() {
	// // Data creation functions
	createSquaredData();
	createNonSquaredData();

	// Squared matrices
	runTest("A/","B/","C/");

	// Non squared matrices
	runTest("A/nonSquared/","B/nonSquared/","C/nonSquared/");
	
	return 0;
}
