// matrix.h: Header file and its contents.

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>

class AugMatrix{
	private:
		int numberOfx; // number of unknowns
		std::vector< std::vector<double> > matrix; // the augmented coefficient matrix
		std::vector<double> xArray; // the solution to the augmented matrix
		void helperFunc(std::string &strnum); // helper function for correct input.
	public:
		AugMatrix(int numx = 1); // constructor for an augmented matrix
		void findX();
};