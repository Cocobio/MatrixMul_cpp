# Matrix multiplication algoritms implemented in C++

### About
Work done for asignment at UdeC - Master in computer science.

### Implementations:

The algoritms included in this repository are the following: 
  - naive multiplication function
  - naive with different for order
  - transpose
  - Strassen purely recursive
  - Strassen hybrid
 
### Input data:

The algorithms receive a pointer to an array, containing pointers to the rows.
The vector_utility.cpp contains a template allocate function that returns an allocated matrix with rows 'r' and columns 'c'.

It also can read and write to binary files using a template function.

Function for creating random matrices are included in the utility file.

#### Binary file structure
The files contain 4 bytes representing the row size, stored as an integer.
4 bytes represeting the column size, also stored as an integer.
And the data in a row-major order.

_Read function requires specify data type._

### How to test:

For testing simply compile and run the "run.cpp" file. You will need 2 matrices stored as binary files to use as input.
