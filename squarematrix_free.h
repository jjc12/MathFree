// matrix.h: SquareMatrix class for non-augmented matrices.
//this class still needs to implement all six "freebies" 
//i.e. move constructor, move assignment operator, etc.

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <cmath>

class SquareMatrix{
	
	friend std::ostream& operator<<(std::ostream& s, SquareMatrix m);

	private:

		int dimension; // size of matrix
		std::vector< std::vector<double> > matrix; // the coefficient matrix
		void helperFunc(std::string &strnum); // helper function for correct input.

	public:

		// constructor for a matrix
		//initializes everything to zero.
		explicit SquareMatrix(int dimension = 1);
		explicit SquareMatrix(double&& dimension) = delete;
		explicit SquareMatrix(float&& dimension) = delete;
		explicit SquareMatrix(char&& dimension) = delete;

		//copy constructor
		SquareMatrix(const SquareMatrix& m);

		//assignment operator overload
		SquareMatrix& operator=(const SquareMatrix& m);

		//accessor operator overload
		//gets value at position - NOT index
		//returns value 0 if out of bounds.
		const double& operator()(int i, int j) const;

		//mutator operator overload
		//sets value at position - NOT index
		//returns false if out of bounds.
		double& operator()(int i, int j);

		//other operator overloads.
		SquareMatrix operator+(const SquareMatrix& matrix2) const;
		SquareMatrix& operator+=(const SquareMatrix& matrix2);
		SquareMatrix operator-(const SquareMatrix& matrix2) const;
		SquareMatrix& operator-=(const SquareMatrix& matrix2);
		SquareMatrix operator*(double scalar);
		SquareMatrix operator*(const SquareMatrix& matrix2) const;
		std::vector<double> operator*(std::vector<double> vVec) const;
		SquareMatrix& operator*=(const SquareMatrix& matrix2);

		//allows the user to input values into the matrix.
		void consoleInput();
		void fileInput(const std::string& filename);

		//displays the matrix on the console
		void show() const;

		//gets dimension of matrix
		int getDimension() const;

		//returns transpose
		SquareMatrix transpose() const;

		//row operations; multiplier multiplies rownum2 and performs
		//a specific operation and modifies the SAME matrix.
		SquareMatrix& rowAdd(int rownum1, int rownum2, double multiplier = 1);
		SquareMatrix& rowNeg(int rownum1, int rownum2, double multiplier = 1);
		SquareMatrix& rowMul(int rownum1, double multiplier);
		SquareMatrix& rowSwap(int rownum1, int rownum2);

		//obtain submatrix by eliminating row, column
		SquareMatrix submatrix(int row, int column) const;

		//transformations for 2d vectors represented in homogenous co-ordinates
		static std::vector<double>
			translation(double x, double y, const std::vector<double>& myVector);
		static std::vector<double>
			horizontalShear(double multiplier, const std::vector<double>& myVector);
		static std::vector<double> 
			verticalShear(double multiplier, const std::vector<double>& myVector);
		static std::vector<double> 
			horizontalFlip(const std::vector<double>& myVector);
		static std::vector<double> 
			verticalFlip(const std::vector<double>& myVector);
		static std::vector<double> 
			squeezeMap(double multiplier, const std::vector<double>& myVector);
		static std::vector<double> 
			rotation(double angle, const std::vector<double>& myVector);

		//row reduce the matrix
		SquareMatrix rowEchelon() const;
		SquareMatrix rowReducedEchelon() const;

		//determinant for nxn square matrices
		//calculates it recursively
		double getDeterminant();

		//returns true if matrix contains L.I. columns
		//false if a column is a multiple of another.
		bool isLinearlyIndependent();

		//calculates the inverse of matrix
		//returns original if no inverse exists
		SquareMatrix inverse();

		//given an nx1 vector b
		//calculates x-transpose in Ax = b
		std::vector<double> xVec(std::vector<double> bVec);

		//calculates the null space (kernel) of the matrix
		std::vector<double> nullSpace();
		

		//LU-factorization
		//std::vector<SquareMatrix> factorLU();

};
