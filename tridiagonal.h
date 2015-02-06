// tridiagonal.h: Header file and its contents. Works for tridiagonal matrices, usually when
// the matrix has diagonally dominant entries.

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>

class TriDiagonal{
	private:
		int numberOfx; // number of unknowns
		std::vector< std::vector<double> > matrix; // the coefficient matrix
		std::vector<double> bVector;
		std::vector<double> diagMatrix;
		std::vector< std::vector<double> > remMatrix;
		void helperFunc(std::string &strnum); // helper function for correct input.
	public:
		TriDiagonal(int numberOfx = 2); // constructor for an augmented matrix
		void SolveJacobi();
};