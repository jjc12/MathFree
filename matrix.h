// matrix.h: Matrix class for non-augmented matrices.
//this class still needs to implement all six "freebies" 
//i.e. move constructor, move assignment operator, etc.

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

class Matrix{

	private:

		int dimension; // size of matrix
		std::vector< std::vector<double> > matrix; // the coefficient matrix
		void helperFunc(std::string &strnum); // helper function for correct input.

	public:

		// constructor for a matrix
		//initializes everything to zero.
		explicit Matrix(int size = 1);

		//copy constructor
		Matrix(const Matrix& m);

		//assignment operator overload
		Matrix operator=(const Matrix& m);

		//accessor operator overload
		//gets value at position - NOT index
		//returns value 0 if out of bounds.
		const double& operator()(int i, int j) const;

		//mutator operator overload
		//sets value at position - NOT index
		//returns false if out of bounds.
		double& operator()(int i, int j);

		//other operator overloads.
		Matrix operator+(const Matrix& matrix2) const;
		Matrix& operator+=(const Matrix& matrix2);
		Matrix operator-(const Matrix& matrix2) const;
		Matrix& operator-=(const Matrix& matrix2);
		Matrix operator*(double scalar);
		Matrix operator*(const Matrix& matrix2) const;
		std::vector<double> operator*(std::vector<double> vVec) const;
		Matrix& operator*=(const Matrix& matrix2);

		//allows the user to input values into the matrix.
		void consoleInput();
		void fileInput(std::string filename);

		//displays the matrix on the console
		void Show() const;

		//gets dimension of matrix
		int getDimension() const;

		//returns transpose
		Matrix Transpose() const;

		//row operations; multiplier multiplies rownum2 and performs
		//a specific operation and modifies the SAME matrix.
		Matrix& rowAdd(int rownum1, int rownum2, double multiplier = 1);
		Matrix& rowNeg(int rownum1, int rownum2, double multiplier = 1);
		Matrix& rowMul(int rownum1, double multiplier);
		Matrix& rowSwap(int rownum1, int rownum2);

		//obtain submatrix by eliminating row, column
		Matrix Submatrix(int row, int column) const;

		//transformations for 2d matrices
		Matrix horizontalShear(double multiplier) const;
		Matrix verticalShear(double multiplier) const;
		Matrix horizontalFlip(double multiplier) const;
		Matrix verticalFlip(double multiplier) const;
		Matrix squeezeMap(double multiplier) const;
		Matrix rotation(double angle) const;

		//row reduce the matrix
		Matrix rowEchelon() const;
		Matrix reducedRowEchelon() const;

		//determinant for nxn square matrices
		//calculates it recursively
		double getDeterminant();

		//returns true if matrix contains L.I. columns
		//false if a column is a multiple of another.
		bool isLinearlyIndependent();

		//calculates the inverse of matrix
		//returns original if no inverse exists
		Matrix inverse();

		//given an nx1 vector b
		//calculates x-transpose in Ax = b
		std::vector<double> xVec(std::vector<double> bVec);
		

		//LU-factorization
		//std::vector<Matrix> factorLU();

};