// matrix.cpp: Source file and its contents.

#include "matrix_free.h"

std::ostream& operator<<(std::ostream& s, const Matrix& m){

	std::cout << "[";

	//output of matrix
	for (int i = 0; i < m.rows; ++i){
		for (int j = 0; j < m.columns; ++j){
			std::cout << m.matrix[i][j];
			if (j != m.columns - 1)
				std::cout << ", ";
		}
		if (i != m.rows - 1)
			std::cout << "; ";
	}
	std::cout << "]";

	return s;

}

std::ostream& operator<<(std::ostream& s, const std::vector<double>& v){

	std::cout << "[";

	//output of matrix
	for (unsigned i = 0; i < v.size(); ++i){
		std::cout << v[i];
		if (i != v.size() - 1)
			std::cout << ", ";
	}
	std::cout << "]";

	return s;

}

void Matrix::isValidInput(std::string &strnum){
	std::string tryagain;
	int bad = 0;
	int watchout = 0;
	int strsize = strnum.size();

	//if invalid number
	for (int j = 0; j < strsize; ++j){
		//if first digit is not a '-' or a '.' or if any digit after that
		//is not a period otherwise invalid. if two periods detected, invalid.
		if ((j == 0 && !isdigit(strnum[j]) && !(strnum[j] == '-')
			&& !(strnum[j] == '.')) || (j != 0
			&& !isdigit(strnum[j]) && !(strnum[j] == '.'))){
			bad = 1;
		}
		if (strnum[j] == '.'){
			++watchout;
		}
	}

	//if bad = 1 or watchout <= 2,
	//redirects you to the beginning of the function, recursively.
	if (bad == 1 || 2 <= watchout){
		std::cout << "Invalid input; please enter numbers. Try again: ";
		std::cin >> tryagain;
		isValidInput(tryagain);
		strnum = tryagain;
	}

	return;
}

void Matrix::changeToZerosIfNeeded(Matrix& sm){

	//changes values if they are below a certain tolerance (due to
	//floating-point precision errors)
	for (int i = 0; i < sm.rows; ++i){
		for (int j = 0; j < sm.columns; ++j){
			if (fabs(sm.matrix[i][j]) < tolerance){
				sm.matrix[i][j] = 0;
			}
		}
	}

}

void Matrix::changeToZerosIfNeeded(std::vector<double>& myV){

	//changes values if they are below a certain tolerance (due to
	//floating-point precision errors)
	for (unsigned i = 0; i < myV.size(); ++i){
		if (fabs(myV[i]) < tolerance){
			myV[i] = 0;
		}
	}

}

Matrix::Matrix(){

	rows = columns = 1;

	matrix = std::vector< std::vector<double> >(rows, std::vector<double>(columns, 0));

}

//constructs a rows-by-columns box of numbers.
Matrix::Matrix(int rows, int columns){

	if (rows <= 0 || columns <= 0){
		this->rows = 1;
		this->columns = 1;
	}
	else{
		this->rows = rows;
		this->columns = columns;
	}

	matrix = std::vector< std::vector<double> >(rows, std::vector<double>(columns, 0));

}

Matrix::~Matrix(){

}

//copy constructs a size-by-size matrix.
Matrix::Matrix(const Matrix& m){

	//in case copy constructor copies to itself.
	if (&m == this){
		rows = 1;
		columns = 1;
	}
	else{
		rows = m.rows;
		columns = m.columns;
	}

	matrix = m.matrix;

}

//conversion constructor (int vector to matrix)
Matrix::Matrix(std::vector<int> myintVec){

	rows = 1;
	columns = myintVec.size();

	std::vector<double> row;
	for (int j = 0; j < columns; ++j){
		row.push_back(myintVec[j]);
	}
	matrix.push_back(row);

}

//conversion constructor (double vector to matrix)
Matrix::Matrix(std::vector<double> mydoubVec){

	rows = 1;
	columns = mydoubVec.size();

	std::vector<double> row = mydoubVec;
	matrix.push_back(row);

}

//conversion constructor (std::initializer_list to matrix)
Matrix::Matrix(std::initializer_list<double> myList){

	rows = 1;
	columns = myList.size();

	matrix.push_back(std::vector<double>(columns));

	std::copy(myList.begin(), myList.end(), matrix[0].begin());

}

//sets two matrices equal (deep copy)
Matrix& Matrix::operator=(const Matrix& m){

	if (&m == this){
		return *this;
	}

	rows = m.rows;
	columns = m.columns;

	//clear the original matrix
	matrix = m.matrix;

	return *this;

}

bool Matrix::operator==(const Matrix& m) const{

	//comparing equality of 2-d vectors
	if (matrix != m.matrix){
		return false;
	}
	else{
		return true;
	}

}

const double& Matrix::operator()(int i, int j) const{

	if (i < 1 || rows < i){
		std::cout << "Position " << i << " is out of bounds; bound is "
			<< rows << std::endl;
		std::cout << "Returning element in row 1, column 1..." << std::endl;
		return matrix[0][0];
	}
	else if (j < 1 || columns < j){
		std::cout << "Position " << j << " is out of bounds; bound is "
			<< columns << std::endl;
		std::cout << "Returning element in row 1, column 1..." << std::endl;
		return matrix[0][0];
	}

	return matrix[i - 1][j - 1];

}

double& Matrix::operator()(int i, int j){

	if (i < 1 || rows < i){
		std::cout << "Position " << i << " is out of bounds; bound is "
			<< rows << std::endl;
		std::cout << "Returning element in row 1, column 1..." << std::endl;
		return matrix[0][0];
	}
	else if (j < 1 || columns < j){
		std::cout << "Position " << j << " is out of bounds; bound is "
			<< columns << std::endl;
		std::cout << "Returning element in row 1, column 1..." << std::endl;
		return matrix[0][0];
	}

	return matrix[i - 1][j - 1];

}

Matrix Matrix::operator+(const Matrix& matrix2) const{

	if (rows != matrix2.rows || columns != matrix2.columns){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return *this;
	}
	//the end summand
	Matrix mProd(rows, columns);

	//performing the matrix multiplication.
	for (int i = 0; i < rows; ++i){
		for (int j = 0; j < columns; ++j){
			double tempvar = matrix[i][j] + matrix2.matrix[i][j];
			mProd.matrix[i][j] = tempvar;
		}
	}
	return mProd;
}

Matrix& Matrix::operator+=(const Matrix& matrix2){
	*this = *this + matrix2;
	return *this;
}

Matrix Matrix::operator-(const Matrix& matrix2) const{
	if (rows != matrix2.rows || columns != matrix2.columns){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return *this;
	}
	//the end subtractand
	Matrix mProd(rows, columns);

	//performing the matrix multiplication.
	for (int i = 0; i < rows; ++i){
		for (int j = 0; j < columns; ++j){
			mProd.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
		}
	}
	return mProd;
}

Matrix& Matrix::operator-=(const Matrix& matrix2){
	*this = *this - matrix2;
	return *this;
}

Matrix Matrix::operator*(double scalar){

	Matrix mProd = *this;

	//performing the scalar multiplication
	for (int i = 0; i < rows; ++i){
		for (int j = 0; j < columns; ++j){
			mProd.matrix[i][j] *= scalar;
		}
	}

	return mProd;

}

Matrix& Matrix::operator*=(double scalar){
	*this = *this * scalar;
	return *this;
}

Matrix operator*(double scalar, const Matrix& matrix2){
	Matrix temp = matrix2;
	return temp * scalar;
}

Matrix Matrix::operator*(const Matrix& matrix2) const{

	if (matrix2.rows == 1 && matrix2.columns == 1){
		return matrix2.matrix[0][0] * *this;
	}
	else if (columns != matrix2.rows){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return *this;
	}

	//the end product
	Matrix mProd(rows, matrix2.columns);

	//performing the matrix multiplication.
	for (int i = 0; i < mProd.rows; ++i){
		for (int j = 0; j < mProd.columns; ++j){
			double tempvar = 0;
			//k < columns; or k < matrix2.rows;
			for (int k = 0; k < columns; ++k){
				tempvar += matrix[i][k] * matrix2.matrix[k][j];
			}
			mProd.matrix[i][j] = tempvar;
		}
	}

	return mProd;

}

Matrix Matrix::operator*(std::vector<double> vVec) const{

	Matrix myVec = vVec;
	return (*this)*myVec;

}

Matrix& Matrix::operator*=(const Matrix& matrix2){
	*this = *this * matrix2;
	return *this;
}


void Matrix::consoleInput(){ //inputs for the nxn matrix.

	std::string acoeff;
	int rowCounter = 0;

	do {
		std::cout << "Enter coefficients of linear expression "
			<< rowCounter + 1 << ":" << std::endl;
		for (int j = 0; j < columns; ++j){
			std::cin >> acoeff;
			isValidInput(acoeff);
			matrix[rowCounter][j] = stod(acoeff);
		}
		++rowCounter;
	} while (rowCounter < rows);

}

//inputs for the nxn matrix.
void Matrix::fileInput(const std::string& filename){

	std::string acoeff;
	int rowCounter = 0;
	std::ifstream inFile(filename);

	do {
		for (int j = 0; j < columns; ++j){
			if (inFile >> acoeff){
				matrix[rowCounter][j] = stod(acoeff);
			}
			else
				break;
		}
		++rowCounter;
	} while (rowCounter < rows);

}

void Matrix::show() const{
	std::cout << "The matrix is " << std::endl;
	for (int i = 0; i < rows; ++i){ //output of matrix
		for (int j = 0; j < columns; ++j){
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	return;
}

int Matrix::getRows() const{
	return rows;
}

int Matrix::getColumns() const{
	return columns;
}


Matrix Matrix::transpose() const{

	Matrix trans(columns, rows);
	for (int i = 0; i < rows; ++i){
		for (int j = 0; j < columns; ++j){
			trans.matrix[j][i] = matrix[i][j];
		}
	}

	return trans;

}

Matrix&
Matrix::rowAdd(int rownum1, int rownum2, double multiplier){

	if (rownum1 < 1 || rows < rownum1){
		std::cout << "Error: Parameter 1 out of bounds" << std::endl;
		return *this;
	}
	else if (rownum2 < 1 || columns < rownum2){
		std::cout << "Error: Parameter 2 out of bounds" << std::endl;
		return *this;
	}

	for (int i = 0; i < columns; ++i){
		matrix[rownum1 - 1][i] += matrix[rownum2 - 1][i] * multiplier;
	}

	return *this;

}

Matrix& Matrix::rowMul(int rownum1, double multiplier){

	if (rownum1 < 1 || rows < rownum1){
		std::cout << "Error: Parameter out of bounds" << std::endl;
		return *this;
	}

	for (int i = 0; i < columns; ++i){
		matrix[rownum1 - 1][i] *= multiplier;
	}

	return *this;

}

Matrix& Matrix::rowSwap(int rownum1, int rownum2){

	if (rownum1 < 1 || rows < rownum1){
		std::cout << "Error: Parameter 1 out of bounds" << std::endl;
		return *this;
	}
	else if (rownum2 < 1 || columns < rownum2){
		std::cout << "Error: Parameter 2 out of bounds" << std::endl;
		return *this;
	}
	else if (rownum1 == rownum2){
		std::cout << "No swap performed: parameters are equal" << std::endl;
		return *this;
	}

	std::vector<double> temp(columns, 0);

	for (int i = 0; i < columns; ++i){
		temp[i] = matrix[rownum1 - 1][i];
		matrix[rownum1 - 1][i] = matrix[rownum2 - 1][i];
		matrix[rownum2 - 1][i] = temp[i];
	}

	return *this;

}

Matrix Matrix::submatrix(int row, int column) const{

	if (row < 1 || rows < row){
		std::cout << "Index " << row << " is out of bounds; bound is "
			<< rows << std::endl;
		return *this;
	}
	else if (column < 1 || columns < column){
		std::cout << "Index " << column << " is out of bounds; bound is "
			<< columns << std::endl;
		return *this;
	}

	if (rows == 1 || columns == 1){
		Matrix temp;
		return temp;
	}

	Matrix subM(rows - 1, columns - 1);

	//skips the deleted row/column
	for (int i = 0, count_one = 0; i < rows; ++i){
		if (i == row - 1){
			++count_one;
			continue;
		}
		for (int j = 0, count_two = 0; j < columns; ++j){
			if (j == column - 1){
				++count_two;
				continue;
			}
			else{
				subM.matrix[i - count_one][j - count_two] = matrix[i][j];
			}
		}
	}

	return subM;

}

Matrix Matrix::translation(double x, double y,
	const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	Matrix mNew(3, 3);

	mNew.matrix[0][0] = 1;
	mNew.matrix[0][1] = 0;
	mNew.matrix[0][2] = x;
	mNew.matrix[1][0] = 0;
	mNew.matrix[1][1] = 1;
	mNew.matrix[1][2] = y;
	mNew.matrix[2][0] = 0;
	mNew.matrix[2][1] = 0;
	mNew.matrix[2][2] = 1;
	return mNew*myVector;
}

Matrix Matrix::horizontalShear(double multiplier,
	const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	Matrix mNew(3, 3);

	mNew.matrix[0][0] = 1;
	mNew.matrix[0][1] = multiplier;
	mNew.matrix[0][2] = 0;
	mNew.matrix[1][0] = 0;
	mNew.matrix[1][1] = 1;
	mNew.matrix[1][2] = 0;
	mNew.matrix[2][0] = 0;
	mNew.matrix[2][1] = 0;
	mNew.matrix[2][2] = 1;

	return mNew*myVector;

}

Matrix Matrix::verticalShear(double multiplier,
	const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	Matrix mNew(3, 3);

	mNew.matrix[0][0] = 1;
	mNew.matrix[0][1] = 0;
	mNew.matrix[0][2] = 0;
	mNew.matrix[1][0] = multiplier;
	mNew.matrix[1][1] = 1;
	mNew.matrix[1][2] = 0;
	mNew.matrix[2][0] = 0;
	mNew.matrix[2][1] = 0;
	mNew.matrix[2][2] = 1;

	return mNew*myVector;

}

Matrix Matrix::horizontalFlip(const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	Matrix mNew(3, 3);

	mNew.matrix[0][0] = -1;
	mNew.matrix[0][1] = 0;
	mNew.matrix[0][2] = 0;
	mNew.matrix[1][0] = 0;
	mNew.matrix[1][1] = 1;
	mNew.matrix[1][2] = 0;
	mNew.matrix[2][0] = 0;
	mNew.matrix[2][1] = 0;
	mNew.matrix[2][2] = 1;

	return mNew*myVector;

}

Matrix Matrix::verticalFlip(const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	Matrix mNew(3, 3);

	mNew.matrix[0][0] = 1;
	mNew.matrix[0][1] = 0;
	mNew.matrix[0][2] = 0;
	mNew.matrix[1][0] = 0;
	mNew.matrix[1][1] = -1;
	mNew.matrix[1][2] = 0;
	mNew.matrix[2][0] = 0;
	mNew.matrix[2][1] = 0;
	mNew.matrix[2][2] = 1;

	return mNew*myVector;

}

Matrix Matrix::squeezeMap(double multiplier,
	const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}
	else if (multiplier == 0){
		std::cout << "Error: Multiplier must be non-zero" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	Matrix mNew(3, 3);

	mNew.matrix[0][0] = multiplier;
	mNew.matrix[0][1] = 0;
	mNew.matrix[0][2] = 0;
	mNew.matrix[1][0] = 0;
	mNew.matrix[1][1] = 1 / multiplier;
	mNew.matrix[1][2] = 0;
	mNew.matrix[2][0] = 0;
	mNew.matrix[2][1] = 0;
	mNew.matrix[2][2] = 1;

	return mNew*myVector;

}

Matrix Matrix::rotation(double angle,
	const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	angle = fmod(angle, 360);

	Matrix mNew(3, 3);

	mNew.matrix[0][0] = cos(angle);
	mNew.matrix[0][1] = sin(angle);
	mNew.matrix[0][2] = 0;
	mNew.matrix[1][0] = -sin(angle);
	mNew.matrix[1][1] = cos(angle); 
	mNew.matrix[1][2] = 0;
	mNew.matrix[2][0] = 0;
	mNew.matrix[2][1] = 0;
	mNew.matrix[2][2] = 1;

	return mNew*myVector;

}

Matrix Matrix::rowEchelon(){

	//temp and *this have the same dimensions
	Matrix temp = *this;

	//inefficient bubble-sort sorting algorithm to place all zero rows
	//at the bottom

	for (int i = 0; i < rows; ++i){

		for (int j = i; j < rows; ++j){

			//breaks if row i is not a zero vector.
			if (temp.matrix[i] != std::vector<double>(temp.matrix[i].size(), 0)
				){
				break;
			}

			//may swap zero vectors with themselves once all zero vectors are
			//placed at the bottom; not a problem, still.
			if (temp.matrix[i] != temp.matrix[j]){
				temp.rowSwap(i + 1, j + 1);
			}

		}

	}

	//the matrix temp transforms to row-echelon form.
	for (int i = 0; i < temp.rows - 1; ++i){

		//swaps the rows if pivot point is zero
		if (temp.matrix[i][i] == 0){
			int switch_ = 0;
			for (int k = i + 1; k < rows; ++k){
				if (temp.matrix[k][i] != 0){
					temp.rowSwap(k + 1, i + 1);
					switch_ = 1;
					break;
				}
			}
			if (switch_ == 0){

				//this means all numbers below and including "pivot" are zero.

				//preliminary check to see if columns to the right of pivot
				//(and below) are equal to zero 
				//(we must check, or else infinite loop of zero matrix occurs)
				bool switch_ = false;
				for (int row = 0; i + row < rows; ++row){
					if (temp.matrix[i + row]
						!= std::vector<double>(temp.matrix[i].size(), 0)){
						switch_ = true;
					}
				}
				if (switch_ == false){
					return temp;
				}

				//move to the next column, treat that element as pivot.
				Matrix temp_(temp.rows - i, temp.columns - i);

				//shift the column by 1 to the right, copy it to temp_.
				for (int row = 0; row < temp_.rows; ++row){
					for (int column = 0; column < temp_.columns - 1;
						++column){
						temp_.matrix[row][column]
							= temp.matrix[i + row][i + 1 + column];
					}
				}

				//perform row echelon on this modified temp matrix
				temp_ = temp_.rowEchelon();

				//set temp's elements equal to temp_'s by shifting
				//the column by 1 to the left.
				for (int row = 0; row < temp_.rows; ++row){
					for (int column = 0; column < temp_.columns - 1;
						++column){
						temp.matrix[i + row][i + 1 + column]
							= temp_.matrix[row][column];
					}
				}

			}

		}
		
		changeToZerosIfNeeded(temp);

		//multiplies entire row by reciprocal of pivot.
		for (int j = i; j < temp.rows; ++j){
			if (temp.matrix[j][i] != 0){
				temp.rowMul(j + 1, 1 / temp.matrix[j][i]);
			}
		}

		//subtracts rows from each other
		//this loop removes all numbers below pivot at row j.
		for (int j = i + 1; j < temp.rows; ++j){
			if (temp.matrix[j][i] != 0){
				temp.rowAdd(j + 1, i + 1, -1);
			}
		}

	}

	changeToZerosIfNeeded(temp);

	//ensures the bottom pivot is equal to one.
	//this check is invalid. change so that it equals the bottom pivot (not the bottom-corner).
	for (int j = 0; j < temp.columns; ++j){

		changeToZerosIfNeeded(temp);

		if (temp.matrix[temp.rows - 1][j] != 0){
			temp.rowMul(temp.rows, 1 / temp.matrix[temp.rows - 1][j]);
			break;
		}
	}

	changeToZerosIfNeeded(temp);

	return temp;

}

Matrix Matrix::rowReducedEchelon(){

	//Get the matrix in row-echelon form
	Matrix temp = rowEchelon();

	//preliminary check to see if top rows are zero rows
	//(we must check, or else infinite loop of zero matrix occurs)
	int zeroRows = 0;
	for (int row = rows - 1; 0 <= row; --row){
		if (temp.matrix[row]
			== std::vector<double>(temp.matrix[row].size(), 0)){
			++zeroRows;
		}
	}
	//it's a zero matrix
	if (zeroRows == temp.rows){
		return temp;
	}

	//create zeros above each pivot
	for (int k = temp.rows - 1 - zeroRows; 1 <= k; --k){
		int zerosBeforePivot = 0;
		//checks where the pivot is
		for (int j = 0; j < temp.columns - 1; ++j){
			if (temp.matrix[k][j] == 0){
				++zerosBeforePivot;
			}
			else{
				break;
			}
		}
		//eliminates numbers above pivot
		for (int l = k - 1, j = zerosBeforePivot; 0 <= l; --l){
			if (temp.matrix[l][j] != 0){
				temp.rowAdd(l + 1, k + 1, -1 * temp.matrix[l][j]);
			}
		}
	}

	return temp;

}

double Matrix::getDeterminant(){

	if (rows != columns){
		std::cout << "Determinant for non-square matrices is undefined.\n";
		std::cout << "Returning zero by default..." << std::endl;
		return 0;
	}

	//assumes rows == columns
	int dimension = rows;

	if (dimension == 1){
		return matrix[0][0];
	}
	else if (dimension == 2){
		return matrix[0][0] * matrix[1][1]
			- matrix[0][1] * matrix[1][0];
	}

	double determinant = 0;

	for (int i = 0; i < dimension; ++i){
		if (matrix[0][i] == 0){
			continue;
		}
		else{
			determinant += matrix[0][i] * submatrix(1, i + 1).getDeterminant()
				* pow(-1, i);
		}
	}

	return determinant;

}

bool Matrix::isLinearlyIndependent(){

	return getDeterminant() == 0 ? false : true;

}

Matrix Matrix::inverse(){

	if (getDeterminant() == 0){
		std::cout << "Determinant is zero ---> Inverse does not exist"
			<< std::endl;
		return *this;
	}
	//Assumes matrix is square
	//thus rows == columns.
	else if (rows == 1){
		Matrix temp;
		temp(1, 1) = 1 / matrix[0][0];
		return temp;
	}

	int dimension = rows;

	//calculate minors
	//also calculate cofactors
	Matrix minors(dimension, dimension);

	for (int i = 0; i < dimension; ++i){
		for (int j = 0; j < dimension; ++j){
			minors.matrix[i][j] = submatrix(i + 1, j + 1).getDeterminant();
			if (i % 2 == 0){
				minors.matrix[i][j] *= pow(-1, j);
			}
			else{
				minors.matrix[i][j] *= pow(-1, j + 1);
			}
		}
	}

	//calculate adjoint (adjugate)
	//multiply by 1/determinant
	//return this inverse
	return minors.transpose()*(1 / getDeterminant());

}

Matrix Matrix::xVec(const Matrix& bVec){

	//here we assume bVec is a 1xn matrix, the transpose of b.

	std::vector<double> xVector(columns, 0);

	//A^-1*b
	if (getDeterminant() != 0){
		return inverse()*bVec.transpose();
	}
	//incompatible dimensions
	else if (rows != bVec.columns){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return xVector;
	}
	//cx = 0
	else if (bVec.rows == 1 && bVec.columns == 1 && bVec.matrix[0][0] == 0){
		xVector.resize(1);
		//0x = 0
		if (matrix[0][0] == 0){
			std::cout << "Any solution works." << std::endl;
		}
		//c*0 = 0
		xVector[0] = 0;
		return xVector;
	}
	//either system is overdetermined or underdetermined.
	//modifies the matrix so rows = columns and calls xVec().
	if (rows != columns){

		//case where system is overdetermined (rows > columns)
		if (rows > columns){
			//find a solution for the first two rows of this matrix.
			Matrix solnMatrix(columns, columns);
			for (int i = 0; i < columns; ++i){
				for (int j = 0; j < columns; ++j){
					solnMatrix.matrix[i][j] = matrix[i][j];
				}
			}

			Matrix modVec(1, columns);
			//copying the vector b.
			for (int i = 0; i < modVec.columns; ++i){
				modVec.matrix[0][i] = bVec.matrix[0][i];
			}
			Matrix soln = solnMatrix.xVec(modVec);

			//substitute this solution into original matrix
			//need to ensure that floating-point error does not occur.
			Matrix testMatrix(rows, 1);
			for (int i = 0; i < rows; ++i){

				double tempvar = 0;
				for (int j = 0; j < columns; ++j){
					tempvar += matrix[i][j] * soln.matrix[j][0];
				}
				//store the residual
				testMatrix.matrix[i][0] = tempvar - bVec.matrix[0][i];

			}
			changeToZerosIfNeeded(testMatrix);

			//if residual of sol'n is not zero, then zero-vector is returned.
			//this implies an inconsistent system.
			for (int i = 0; i < rows; ++i){
				if (testMatrix.matrix[i][0] != 0){
					//elements are zero by default
					std::cout << "Inconsistent system detected;";
					std::cout << " matrix has no solution. " << std::endl;
					std::cout << "Returning default zero vector." << std::endl;
					return Matrix(soln.rows, 1);
				}
			}

			//else it returns the original soln, since it satisfies the system.
			return soln;

		}


		//case where system is underdetermined (rows < columns)
		if (rows < columns){

			//copying smaller matrix to larger one.
			Matrix tempbeforetemp(columns, columns);
			for (int i = 0; i < rows; ++i){
				for (int j = 0; j < columns; ++j){
					tempbeforetemp.matrix[i][j] = matrix[i][j];
				}
			}

			//appending zero-rows to matrix.
			for (int i = rows; i < columns; ++i){
				for (int j = 0; j < columns; ++j){
					tempbeforetemp.matrix[i][j] = 0;
				}
			}

			Matrix modVec(1, columns);

			for (int i = 0; i < bVec.columns; ++i){
				modVec.matrix[0][i] = bVec.matrix[0][i];
			}
			for (int i = bVec.columns; i < columns; ++i){
				modVec.matrix[0][i] = 0;
			}

			Matrix soln = tempbeforetemp.xVec(modVec);

			return soln;

		}

	}

	//**This also assumes the matrix is square**
	//let dimension == rows since rows == columns.
	//Assume the matrix is a (n+1)x(n+1) with a zero row at the bottom.
	//the nx1th column is the b-vector.
	//copying the nxn portion
	int dimension = rows;

	Matrix temp(dimension + 1, dimension + 1);
	for (int i = 0; i < dimension; ++i){
		for (int j = 0; j < dimension; ++j){
			temp.matrix[i][j] = matrix[i][j];
		}
	}
	//copying the nx1th column
	for (int i = 0; i < temp.columns - 1; ++i){
		temp.matrix[i][temp.columns - 1] = bVec.matrix[0][i];
	}

	temp = temp.rowEchelon();

	//Augmented matrix A is in row-echelon
	//now we are ready for back substitution

	//the last index
	int start = temp.rows - 2;

	//detect if A is consistent by making sure no pivot exists in b-column.
	//initialize a zero-vector for testing
	std::vector<double> zero_vector(temp.columns - 1, 0);
	for (int i = 0; i < temp.columns - 1; ++i){
		//this vector gets the elements in the range 0 to n-1 of the matrix
		//sees if any of them equal zero for row i.
		std::vector<double> test_vector(temp.columns - 1, 0);
		for (int j = 0; j < temp.columns - 1; ++j){
			test_vector[j] = temp.matrix[i][j];
		}
		if (test_vector == zero_vector
			&& temp.matrix[i][temp.columns - 1] != 0){
			std::cout << "Inconsistent system detected;";
			std::cout << " matrix has no solution. " << std::endl;
			std::cout << "Returning default zero vector." << std::endl;
			std::vector<double> newVec(temp.columns - 1, 0);
			return newVec;
		}
	}

	//detect which variables of A are free
	//by checking which columns aren't pivot columns
	std::vector<bool> pivotColumns(temp.columns - 1, false);
	for (int i = 0; i < temp.rows - 1; ++i){
		for (int j = 0; j < temp.columns - 1; ++j){
			if (temp.matrix[i][j] != 0){
				pivotColumns[j] = true;
				break;
			}
		}
	}

	//runs through pivotColumns; false means the variable is free
	for (int i = 0; i < temp.columns - 1; ++i){
		if (pivotColumns[i] == false){
			std::cout << "x_" << i + 1 << " is free;";
			std::cout << " setting it equal to one..." << std::endl;
			xVector[i] = 1;
		}
	}

	//b - Ax to find the next x-value
	//finds first non-zero row from the bottom.
	int subtractedRow = 0;
	for (int i = start; 0 <= i; --i){
		if (temp.matrix[i] == std::vector<double>(temp.columns, 0)){
			++subtractedRow;
		}
	}

	//row i
	//pivot is assumed to equal 1
	bool breaker = false;
	for (int i = start - subtractedRow, k = xVector.size() - 1; 0 <= i; --i, --k){
		//if x_k is free, skip over that x.
		//make sure k is not negative.
		while (xVector[k] == 1){
			--k;
			if (k < 0){
				breaker = true;
				break;
			}
		}
		if (breaker == true){
			break;
		}
		double _Ax = 0;

		//calculate _Ax in row i
		for (int j = start; k < j; --j){
			_Ax += temp.matrix[i][j] * xVector[j];

		}

		xVector[k] = temp.matrix[i][temp.columns - 1] - _Ax;
	
	}

	changeToZerosIfNeeded(xVector);

	return xVector;

}

Matrix Matrix::nullSpace(){
	std::vector<double> nullVec(columns, 0);
	return xVec(nullVec);
}

/*std::vector<Matrix> Matrix::factorLU(){

Matrix L(dimension);
Matrix U(dimension);

for (int i = 0; i < dimension; ++i){
L.matrix[i][i] = 1;
}

std::vector<Matrix> mArray;

mArray.push_back(L);
mArray.push_back(U);

return mArray;

}*/
