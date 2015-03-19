// matrix.cpp: Source file and its contents.

#include "matrix.h"

void Matrix::helperFunc(std::string &strnum){
	std::string tryagain;
	int bad = 0;
	int watchout = 0;
	int strsize = strnum.size();

	for (int i = 0; i < strsize; ++i){
		if (!((isdigit(strnum[i])) || (strnum[i] == '.'))){
			bad = 1;
			if ((strnum[0] == '-') && (1 < strsize)){
				for (int j = 1; j < strsize; ++j){
					if (strnum[j] == '-'){
						goto label;
					}
				}
				bad = 0;
				continue;
			}
		label:
			break;
		}
		if (strnum[i] == '.'){
			++watchout;
			if (watchout == 2){
				break;
			}
		}
	}

	if ((bad == 1) || (2 == watchout)){
		std::cout << "Invalid input; please enter numbers. Try again: ";
		std::cin >> tryagain;
		helperFunc(tryagain); 
		strnum = tryagain;
	}
	return;
}


Matrix::Matrix(int size){ //constructs a size-by-size matrix.

	if (size <= 0){
		dimension = 1;
	}
	else{
		dimension = size;
	}

	int rowCounter = 0;

	do {
		std::vector<double> row;
		for (int j = 0; j < dimension; ++j){
			row.push_back(0);
		}
		matrix.push_back(row);
		++rowCounter;
	} while (rowCounter < dimension);

}

Matrix::Matrix(const Matrix& m){ //copy constructs a size-by-size matrix.

	dimension = m.dimension;
	int rowCounter = 0;

	do {
		std::vector<double> row;
		for (int j = 0; j < dimension; ++j){
			row.push_back(0);
		}
		matrix.push_back(row);
		++rowCounter;
	} while (rowCounter < dimension);

	for (int i = 0; i < dimension; i++){
		for (int j = 0; j < dimension; j++){
			matrix[i][j] = m.matrix[i][j];
		}
	}

}

Matrix Matrix::operator=(const Matrix& m){ //sets two matrices equal (deep copy)

	dimension = m.dimension;
	Matrix rM(dimension);

	for (int i = 0; i < dimension; i++){
		for (int j = 0; j < dimension; j++){
			rM.matrix[i][j] = m.matrix[i][j];
		}
	}

	return rM;

}

const double& Matrix::operator()(int i, int j) const{

	if (i < 1 || dimension < i){
		std::cout << "Index " << i << " is out of bounds; bound is " << dimension << std::endl;
		std::cout << "Returning element in row 1, column 1..." << std::endl;
		return matrix[0][0];
	}
	else if (j < 1 || dimension < j){
		std::cout << "Index " << j << " is out of bounds; bound is " << dimension << std::endl;
		std::cout << "Returning element in row 1, column 1..." << std::endl;
		return matrix[0][0];
	}

	return matrix[i - 1][j - 1];

}

double& Matrix::operator()(int i, int j){

	if (i < 1 || dimension < i){
		std::cout << "Index " << i << " is out of bounds; bound is " << dimension << std::endl;
		std::cout << "Returning element in row 1, column 1..." << std::endl;
		return matrix[0][0];
	}
	else if (j < 1 || dimension < j){
		std::cout << "Index " << j << " is out of bounds; bound is " << dimension << std::endl;
		std::cout << "Returning element in row 1, column 1..." << std::endl;
		return matrix[0][0];
	}

	return matrix[i - 1][j - 1];

}

Matrix Matrix::operator+(const Matrix& matrix2) const{
	
	if (dimension != matrix2.dimension){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return *this;
	}
	//the end summand
	Matrix mProd(dimension);

	//performing the matrix multiplication.
	for (int i = 0; i < dimension; ++i){
		for (int j = 0; j < dimension; ++j){
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
	if (dimension != matrix2.dimension){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return *this;
	}
	//the end subtractand
	Matrix mProd(dimension);

	for (int i = 0; i < dimension; ++i){ //performing the matrix multiplication.
		for (int j = 0; j < dimension; ++j){
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
	for (int i = 0; i < dimension; ++i){
		for (int j = 0; j < dimension; ++j){
			mProd.matrix[i][j] *= scalar;
		}
	}

	return mProd;

}

Matrix Matrix::operator*(const Matrix& matrix2) const{

	if (matrix2.dimension == 1){
		double scalar = matrix2.matrix[0][0];
		Matrix temp = *this;
		return temp*scalar;
	}
	else if (dimension != matrix2.dimension){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return *this;
	}

	//the end product
	Matrix mProd(dimension);

	//performing the matrix multiplication.
	for (int i = 0; i < dimension; ++i){
		for (int j = 0; j < dimension; ++j){
			double tempvar = 0;
			for (int k = 0; k < dimension; ++k){
				tempvar += matrix[i][k] * matrix2.matrix[k][j];
			}
			mProd.matrix[i][j] = tempvar;
		}
	}
	return mProd;

}

std::vector<double> Matrix::operator*(std::vector<double> vVec) const{

	std::vector<double> mVec(dimension, 0);

	if (dimension != vVec.size()){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return mVec;
	}

	//the end product
	Matrix mProd(dimension);

	//performing the matrix multiplication.
	for (int i = 0; i < dimension; ++i){
		double tempvar = 0;
		for (int j = 0; j < dimension; ++j){
			tempvar += matrix[i][j] * vVec[j];
		}
		mVec[i] = tempvar;
	}

	return mVec;

}

Matrix& Matrix::operator*=(const Matrix& matrix2){
	*this = *this * matrix2;
	return *this;
}


void Matrix::consoleInput(){ //inputs for the nxn matrix.

	std::string acoeff;
	int rowCounter = 0;

	do {
		std::cout << "Enter coefficients of linear equation " << rowCounter + 1 << ":" << std::endl;
		for (int j = 0; j < dimension; ++j){
			std::cin >> acoeff;
			helperFunc(acoeff);
			matrix[rowCounter][j] = stod(acoeff);
		}
		++rowCounter;
	} while (rowCounter < dimension);

}

void Matrix::fileInput(std::string filename){ //inputs for the nxn matrix.

	std::string acoeff;
	int rowCounter = 0;
	std::ifstream inFile(filename);

	do {
		for (int j = 0; j < dimension; ++j){
			if (inFile >> acoeff){
				matrix[rowCounter][j] = stod(acoeff);
			}
			else
				break;
		}
		++rowCounter;
	} while (rowCounter < dimension);

}

void Matrix::Show() const{
	std::cout << "The matrix is " << std::endl;
	for (int i = 0; i < dimension; ++i){ //output of matrix
		for (int j = 0; j < dimension; ++j){
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	return;
}

int Matrix::getDimension() const{
	return dimension;
}

Matrix Matrix::Transpose() const{

	Matrix trans(dimension);
	for (int i = 0; i < dimension; ++i){
		for (int j = 0; j < dimension; ++j){
			trans.matrix[j][i] = matrix[i][j];
		}
	}

	return trans;

}

Matrix& Matrix::rowAdd(int rownum1, int rownum2, double multiplier){

	if (rownum1 < 1 || dimension < rownum1){
		std::cout << "Error: Parameter 1 out of bounds" << std::endl;
		return *this;
	}
	else if (rownum2 < 1 || dimension < rownum2){
		std::cout << "Error: Parameter 2 out of bounds" << std::endl;
		return *this;
	}

	for (int i = 0; i < dimension; ++i){
		matrix[rownum1 - 1][i] += matrix[rownum2 - 1][i] * multiplier;
	}

	return *this;

}

Matrix& Matrix::rowNeg(int rownum1, int rownum2, double multiplier){

	if (rownum1 < 1 || dimension < rownum1){
		std::cout << "Error: Parameter 1 out of bounds" << std::endl;
		return *this;
	}
	else if (rownum2 < 1 || dimension < rownum2){
		std::cout << "Error: Parameter 2 out of bounds" << std::endl;
		return *this;
	}

	for (int i = 0; i < dimension; ++i){
		matrix[rownum1 - 1][i] -= matrix[rownum2 - 1][i] * multiplier;
	}

	return *this;

}

Matrix& Matrix::rowMul(int rownum1, double multiplier){

	if (rownum1 < 1 || dimension < rownum1){
		std::cout << "Error: Parameter out of bounds" << std::endl;
		return *this;
	}

	for (int i = 0; i < dimension; ++i){
		matrix[rownum1 - 1][i] *= multiplier;
	}

	return *this;

}

Matrix& Matrix::rowSwap(int rownum1, int rownum2){

	if (rownum1 < 1 || dimension < rownum1){
		std::cout << "Error: Parameter 1 out of bounds" << std::endl;
		return *this;
	}
	else if (rownum2 < 1 || dimension < rownum2){
		std::cout << "Error: Parameter 2 out of bounds" << std::endl;
		return *this;
	}
	else if (rownum1 == rownum2){
		std::cout << "No swap performed: parameters are equal" << std::endl;
	}

	std::vector<double> temp(dimension, 0);

	for (int i = 0; i < dimension; ++i){
		temp[i] = matrix[rownum1 - 1][i];
		matrix[rownum1 - 1][i] = matrix[rownum2 - 1][i];
		matrix[rownum2 - 1][i] = temp[i];
	}

	return *this;

}

Matrix Matrix::Submatrix(int row, int column) const{

	if (row < 1 || dimension < row){
		std::cout << "Index " << row << " is out of bounds; bound is " << dimension << std::endl;
		return *this;
	}
	else if (column < 1 || dimension < column){
		std::cout << "Index " << column << " is out of bounds; bound is " << dimension << std::endl;
		return *this;
	}

	Matrix subM(dimension - 1);

	for (int i = 0, count_one = 0; i < dimension; ++i, ++count_one){
		if (i == row - 1){
			--count_one;
			continue;
		}
		for (int j = 0, count_two = 0; j < dimension; ++j, ++count_two){
			if (j == column - 1){
				--count_two;
				continue;
			}
			else{
				subM.matrix[count_one][count_two] = matrix[i][j];
			}
		}
	}

	return subM;

}

Matrix Matrix::horizontalShear(double multiplier) const{

	if (dimension != 2){
		std::cout << "Error: Must be a 2-dimensional matrix" << std::endl;
		return *this;
	}

	Matrix mNew(dimension);

	mNew.matrix[0][0] = 1;
	mNew.matrix[0][1] = multiplier;
	mNew.matrix[1][0] = 0;
	mNew.matrix[1][1] = 1;

	return mNew**this;

}

Matrix Matrix::verticalShear(double multiplier) const{

	if (dimension != 2){
		std::cout << "Error: Must be a 2-dimensional matrix" << std::endl;
		return *this;
	}

	Matrix mNew(dimension);

	mNew.matrix[0][0] = 1;
	mNew.matrix[0][1] = 0;
	mNew.matrix[1][0] = multiplier;
	mNew.matrix[1][1] = 1;

	return mNew**this;

}

Matrix Matrix::horizontalFlip(double multiplier) const{

	if (dimension != 2){
		std::cout << "Error: Must be a 2-dimensional matrix" << std::endl;
		return *this;
	}

	Matrix mNew(dimension);

	mNew.matrix[0][0] = -1;
	mNew.matrix[0][1] = 0;
	mNew.matrix[1][0] = 0;
	mNew.matrix[1][1] = 1;

	return mNew**this;

}

Matrix Matrix::verticalFlip(double multiplier) const{

	if (dimension != 2){
		std::cout << "Error: Must be a 2-dimensional matrix" << std::endl;
		return *this;
	}

	Matrix mNew(dimension);

	mNew.matrix[0][0] = 1;
	mNew.matrix[0][1] = 0;
	mNew.matrix[1][0] = 0;
	mNew.matrix[1][1] = -1;

	return mNew**this;

}

Matrix Matrix::squeezeMap(double multiplier) const{

	if (dimension != 2){
		std::cout << "Error: Must be a 2-dimensional matrix" << std::endl;
		return *this;
	}
	else if (multiplier == 0){
		std::cout << "Error: Must be a non-zero parameter" << std::endl;
		return *this;
	}

	Matrix mNew(dimension);

	mNew.matrix[0][0] = multiplier;
	mNew.matrix[0][1] = 0;
	mNew.matrix[1][0] = 0;
	mNew.matrix[1][1] = 1 / multiplier;

	return mNew**this;

}

Matrix Matrix::rotation(double angle) const{

	if (dimension != 2){
		std::cout << "Error: Must be a 2-dimensional matrix" << std::endl;
		return *this;
	}
	else if (fmod(angle, 360) == 0){
		std::cout << "No rotation took place; returning original..." << std::endl;
		return *this;
	}

	angle *= atan(1.0) * 4 / 180;

	Matrix mNew(dimension);

	mNew.matrix[0][0] = cos(angle);
	mNew.matrix[0][1] = -sin(angle);
	mNew.matrix[1][0] = sin(angle);
	mNew.matrix[1][1] = cos(angle);

	return mNew**this;

}

Matrix Matrix::rowEchelon() const{

	Matrix temp = *this;

	//row reduces the matrix temp.
	for (int i = 0; i < dimension - 1; ++i){

		//swaps the rows if pivot point is zero
		if (temp.matrix[i][i] == 0){
			for (int k = i + 1; k < dimension; ++k){
				if (temp.matrix[k][i] != 0){
					temp.rowSwap(k + 1, i + 1);
					break;
				}
			}
		}

		//multiplies entire row by reciprocal of pivot.
		for (int j = i; j < dimension; ++j){
			if (temp.matrix[j][i] != 0){
				temp.rowMul(j + 1, 1 / temp.matrix[j][i]);
			}
		}

		//subtracts rows from each other
		for (int j = i + 1; j < dimension; ++j){
			if (temp.matrix[j][i] != 0){
				temp.rowNeg(j + 1, i + 1);
			}
		}

	}

	if (temp.matrix[dimension - 1][dimension - 1] != 0)
		temp.matrix[dimension - 1][dimension - 1] /= temp.matrix[dimension - 1][dimension - 1];

	return temp;

}

Matrix Matrix::reducedRowEchelon() const{

	Matrix temp = rowEchelon();
	Matrix tempT = temp.Transpose();
	return tempT.rowEchelon();
	
}

double Matrix::getDeterminant(){

	if (dimension == 1){
		return matrix[0][0];
	}
	else if (dimension == 2){
		return matrix[0][0] * matrix[1][1]
			- matrix[0][1] * matrix[1][0];
	}

	double determinant = 0;

	for (int i = 0; i < dimension; ++i){
		determinant += matrix[0][i] * Submatrix(1, i + 1).getDeterminant() * pow(-1, i);
	}

	return determinant;

}

bool Matrix::isLinearlyIndependent(){

	Matrix temp = reducedRowEchelon();
	for (int i = 0; i < temp.dimension; ++i){
		if (temp.matrix[i][i] != 1){
			return false;
		}
	}

	return true;

}

Matrix Matrix::inverse(){

	if (getDeterminant() == 0){
		std::cout << "Determinant is zero ---> Inverse does not exist"
			<< std::endl;
		return *this;
	}

	//calculate minors
	//also calculate cofactors
	Matrix minors(dimension);

	for (int i = 0; i < dimension; ++i){
		for (int j = 0; j < dimension; ++j){
			minors.matrix[i][j] = Submatrix(i + 1, j + 1).getDeterminant();
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
	return minors.Transpose()*(1/getDeterminant());

}

std::vector<double> Matrix::xVec(std::vector<double> bVec){

	Matrix temp = *this;
	std::vector<double> xVector(dimension, 0);

	//A^-1*b
	if (getDeterminant() != 0){
		Matrix temp = *this;
		return temp.inverse()*bVec;
	}
	else if (dimension != bVec.size()){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return xVector;
	}
	else if (bVec.size() == 1 && bVec[0] == 0){
		xVector.resize(1);
		xVector[0] = 0;
		return xVector;
	}

	//row reduces the matrix temp.
	for (int i = 0; i < dimension - 1; ++i){

		//swaps the rows if pivot point is zero
		if (temp.matrix[i][i] == 0){
			for (int k = i + 1; k < dimension; ++k){
				if (temp.matrix[k][i] != 0){
					temp.rowSwap(k + 1, i + 1);
					std::swap(bVec[k], bVec[i]);
					break;
				}
			}
		}

		//multiplies entire row by reciprocal of pivot.
		for (int j = i; j < dimension; ++j){
			if (temp.matrix[j][i] != 0){
				//bVec[j] reassignment MUST go before
				//temp reassignment.
				bVec[j] /= temp.matrix[j][i];
				temp.rowMul(j + 1, 1 / temp.matrix[j][i]);
			}
		}

		//subtracts rows from each other
		for (int j = i + 1; j < dimension; ++j){
			if (temp.matrix[j][i] != 0){
				temp.rowNeg(j + 1, i + 1);
				bVec[j] -= bVec[i];
			}
		}

	}

	if (temp.matrix[dimension - 1][dimension - 1] != 0){
		bVec[dimension - 1] /= temp.matrix[dimension - 1][dimension - 1];
		temp.matrix[dimension - 1][dimension - 1] /= temp.matrix[dimension - 1][dimension - 1];
	}

	//matrix A is row-reduced and b is consistent with these operations
	//now we are ready for back substitution

	//the last index
	int start = dimension - 1;

	for (int k = 0; k < dimension; ++k){
		double restOfRowValue = 0;
		int counter = 0;
		//detects if equation has a free variable.
		if (bVec[start - k] == 0){
			for (int i = 0; i < dimension; ++i){
				if (temp.matrix[start - k][i] == 0){
					++counter;
				}
				else
					break;
			}
			if (counter == dimension){
				std::cout << "x_" << start - k + 1 << " is free;";
				std::cout << " setting it equal to zero..." << std::endl;
				xVector[start - k] = 0;
				continue;
			}
		}

		//detects if equation is inconsistent
		else if (temp.matrix[start - k][start - k] == 0){
			std::cout << "Inconsistent system detected;";
			std::cout << " matrix has only trivial solution." << std::endl;
			std::vector<double> newVec(dimension, 0);
			return newVec;
		}

		//b - Ax to find the next x-value
		for (int j = start; start - k < j; --j){
			restOfRowValue += temp.matrix[start - k][j]*xVector[j];
		}

		xVector[start - k] = bVec[start - k] - restOfRowValue;

	}

	return xVector;

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