// matrix.h: Header file and its contents.

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>

class Matrix{
	private:
		int dimension; // size of matrix
		std::vector< std::vector<double> > matrix; // the coefficient matrix
		void helperFunc(std::string &strnum); // helper function for correct input.
	public:
		Matrix(int size = 1); // constructor for an augmented matrix
		Matrix Multiply(Matrix matrix2);
		double getValue(int i, int j) const;
		void setValue(int i, int j, double value);
};