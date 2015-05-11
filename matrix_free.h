//matrix_free.h: Matrix class for non-augmented matrices.
//this class still needs to implement all six "freebies" 
//i.e. move constructor, move assignment operator, etc.

#ifndef _MATRIX_FREE__H
#define _MATRIX_FREE__H

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <cmath>

class Matrix{

	//output overloads for ostream
	friend std::ostream& operator<<(std::ostream& s, const Matrix& m);
	friend std::ostream& operator<<(std::ostream& s, const std::vector<double>& v);

	//left-scalar multiplication e.g. 4*matrixObject;
	friend Matrix operator*(double scalar, const Matrix& matrix2);

protected:

	// dimensions of matrix
	int rows;
	int columns;
	
	// the coefficient matrix
	std::vector< std::vector<double> > matrix;

	// helper function for correct input.
	void isValidInput(std::string &strnum);

	// helper function to change numbers of magnitude <1e-12 to zero.
	void changeToZerosIfNeeded(Matrix&);
	// override
	void changeToZerosIfNeeded(std::vector<double>& myV);

	// tolerance for said change
	const double tolerance = 1e-12;

public:

	//constructor for a matrix
	//initializes everything to zero.
	Matrix();
	Matrix(int rows, int columns);
	//disallow non-integer dimensions
	explicit Matrix(double&& rows, double&& columns) = delete;
	explicit Matrix(float&& rows, float&& columns) = delete;
	explicit Matrix(char&& rows, char&& columns) = delete;

	//destructor
	~Matrix();

	//copy constructor
	Matrix(const Matrix& m);

	//conversion constructor (int vector to matrix)
	//NOTE: Doesn't have to be explicit.
	Matrix(std::vector<int> myintVec);

	//conversion constructor (double vector to matrix)
	//NOTE: Doesn't have to be explicit.
	Matrix(std::vector<double> mydoubVec);

	//conversion constructor (std::initializer_list to matrix)
	Matrix::Matrix(std::initializer_list<double> myList);

	//assignment operator overload
	Matrix& operator=(const Matrix& m);
	//comparison operator overload
	bool operator==(const Matrix& m) const;


	//accessor operator overload
	//gets value at position - NOT index
	//returns value at row 1, column 1 if out of bounds
	const double& operator()(int i, int j) const;

	//mutator operator overload
	//sets value at position - NOT index
	//modifies value at row 1, column 1 if out of bounds
	double& operator()(int i, int j);

	//other operator overloads.
	Matrix operator+(const Matrix& matrix2) const;
	Matrix& operator+=(const Matrix& matrix2);
	Matrix operator-(const Matrix& matrix2) const;
	Matrix& operator-=(const Matrix& matrix2);
	Matrix operator*(double scalar);
	Matrix& operator*=(double scalar);
	//left-scalar multiplication written as friend
	Matrix operator*(const Matrix& matrix2) const;
	Matrix& operator*=(const Matrix& matrix2);
	Matrix operator*(std::vector<double> vVec) const;

	//allows the user to input values into the matrix.
	void consoleInput();
	void fileInput(const std::string& filename);

	//displays the matrix on the console
	void show() const;

	//gets the number of rows of matrix
	int getRows() const;

	//gets the number of columns of matrix
	int getColumns() const;

	//returns transpose
	Matrix transpose() const;

	//row operations; multiplier multiplies rownum2 and performs
	//a specific operation and modifies the SAME matrix.
	Matrix& rowAdd(int rownum1, int rownum2, double multiplier = 1);
	Matrix& rowNeg(int rownum1, int rownum2, double multiplier = 1);
	Matrix& rowMul(int rownum1, double multiplier);
	Matrix& rowSwap(int rownum1, int rownum2);

	//obtain submatrix by eliminating row, column
	Matrix submatrix(int row, int column) const;

	//transformations for 2d vectors represented in homogenous co-ordinates
	static Matrix
		translation(double x, double y, const std::vector<double>& myVector);
	static Matrix
		horizontalShear(double multiplier, const std::vector<double>& myVector);
	static Matrix
		verticalShear(double multiplier, const std::vector<double>& myVector);
	static Matrix
		horizontalFlip(const std::vector<double>& myVector);
	static Matrix
		verticalFlip(const std::vector<double>& myVector);
	static Matrix
		squeezeMap(double multiplier, const std::vector<double>& myVector);
	static Matrix
		rotation(double angle, const std::vector<double>& myVector);

	//row reduce the matrix
	Matrix rowEchelon();
	Matrix rowReducedEchelon();

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
	Matrix xVec(const Matrix& bVec);

	//calculates the null space (kernel) of the matrix
	Matrix nullSpace();


	//LU-factorization
	//std::vector<Matrix> factorLU();

};

#endif
