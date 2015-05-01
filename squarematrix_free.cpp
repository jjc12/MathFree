// matrix.cpp: Source file and its contents.

#include "squarematrix_free.h"

std::ostream& operator<<(std::ostream& s, SquareMatrix m){

	std::cout << "[";

	for (int i = 0; i < m.dimension; ++i){ //output of matrix
		for (int j = 0; j < m.dimension; ++j){
			std::cout << m.matrix[i][j];
			if (j != m.dimension - 1)
				 std::cout << ", ";
		}
		if (i != m.dimension - 1)
			std::cout << "; ";
	}
	std::cout << "]";

	return s;

}

std::ostream& operator<<(std::ostream& s, std::vector<double> v){

	std::cout << "[";

	for (unsigned i = 0; i < v.size(); ++i){ //output of matrix
		std::cout << v[i];
		if (i != v.size() - 1)
			std::cout << ", ";
	}
	std::cout << "]";

	return s;

}

void SquareMatrix::helperFunc(std::string &strnum){
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

	//if bad = 1 or watchout <= 2, redirects you to the beginning of the function.
	if (bad == 1 || 2 <= watchout){
		std::cout << "Invalid input; please enter numbers. Try again: ";
		std::cin >> tryagain;
		helperFunc(tryagain); 
		strnum = tryagain;
	}

	return;
}

void SquareMatrix::changeToZerosIfNeeded(SquareMatrix& sm){

	//changes values if they are below a certain tolerance (due to
	//floating-point precision errors)
	for (int i = 0; i < sm.dimension; ++i){
		for (int j = 0; j < sm.dimension; ++j){
			if (fabs(sm.matrix[i][j]) < tolerance){
				sm.matrix[i][j] = 0;
			}
		}
	}

}


SquareMatrix::SquareMatrix(int dimension){ //constructs a size-by-size matrix.

	if (dimension <= 0){
		this->dimension = 1;
	}
	else{
		this->dimension = dimension;
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

SquareMatrix::~SquareMatrix(){
	
}

//copy constructs a size-by-size matrix.
SquareMatrix::SquareMatrix(const SquareMatrix& m){

	//in case copy constructor copies to itself.
	if (&m == this){
		dimension = 1;
	}
	else{
		dimension = m.dimension;
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

	for (int i = 0; i < dimension; i++){
		for (int j = 0; j < dimension; j++){
			matrix[i][j] = m.matrix[i][j];
		}
	}

}

//sets two matrices equal (deep copy)
SquareMatrix& SquareMatrix::operator=(const SquareMatrix& m){

	dimension = m.dimension;
	
	//clear the original matrix
	matrix.clear();
	matrix.resize(dimension, std::vector<double>(dimension));

	for (int i = 0; i < dimension; i++){
		for (int j = 0; j < dimension; j++){
			matrix[i][j] = m.matrix[i][j];
		}
	}

	return *this;

}

bool SquareMatrix::operator==(const SquareMatrix& m) const{
	if (matrix != m.matrix){
		return false;
	}
	else{
		return true;
	}
}

const double& SquareMatrix::operator()(int i, int j) const{

	if (i < 1 || dimension < i){
		std::cout << "Index " << i << " is out of bounds; bound is "
			<< dimension << std::endl;
		std::cout << "Returning element in row 1, column 1..." << std::endl;
		return matrix[0][0];
	}
	else if (j < 1 || dimension < j){
		std::cout << "Index " << j << " is out of bounds; bound is "
			<< dimension << std::endl;
		std::cout << "Returning element in row 1, column 1..." << std::endl;
		return matrix[0][0];
	}

	return matrix[i - 1][j - 1];

}

double& SquareMatrix::operator()(int i, int j){

	if (i < 1 || dimension < i){
		std::cout << "Index " << i << " is out of bounds; bound is "
			<< dimension << std::endl;
		std::cout << "Returning element in row 1, column 1..." << std::endl;
		return matrix[0][0];
	}
	else if (j < 1 || dimension < j){
		std::cout << "Index " << j << " is out of bounds; bound is "
			<< dimension << std::endl;
		std::cout << "Returning element in row 1, column 1..." << std::endl;
		return matrix[0][0];
	}

	return matrix[i - 1][j - 1];

}

SquareMatrix SquareMatrix::operator+(const SquareMatrix& matrix2) const{
	
	if (dimension != matrix2.dimension){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return *this;
	}
	//the end summand
	SquareMatrix mProd(dimension);

	//performing the matrix multiplication.
	for (int i = 0; i < dimension; ++i){
		for (int j = 0; j < dimension; ++j){
			double tempvar = matrix[i][j] + matrix2.matrix[i][j];
			mProd.matrix[i][j] = tempvar;
		}
	}
	return mProd;
}

SquareMatrix& SquareMatrix::operator+=(const SquareMatrix& matrix2){
	*this = *this + matrix2;
	return *this;
}

SquareMatrix SquareMatrix::operator-(const SquareMatrix& matrix2) const{
	if (dimension != matrix2.dimension){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return *this;
	}
	//the end subtractand
	SquareMatrix mProd(dimension);

	//performing the matrix multiplication.
	for (int i = 0; i < dimension; ++i){
		for (int j = 0; j < dimension; ++j){
			mProd.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
		}
	}
	return mProd;
}

SquareMatrix& SquareMatrix::operator-=(const SquareMatrix& matrix2){
	*this = *this - matrix2;
	return *this;
}

SquareMatrix SquareMatrix::operator*(double scalar){

	SquareMatrix mProd = *this; 

	//performing the scalar multiplication
	for (int i = 0; i < dimension; ++i){
		for (int j = 0; j < dimension; ++j){
			mProd.matrix[i][j] *= scalar;
		}
	}

	return mProd;

}

SquareMatrix SquareMatrix::operator*(const SquareMatrix& matrix2) const{

	if (matrix2.dimension == 1){
		double scalar = matrix2.matrix[0][0];
		SquareMatrix temp = *this;
		return temp*scalar;
	}
	else if (dimension != matrix2.dimension){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return *this;
	}

	//the end product
	SquareMatrix mProd(dimension);

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

std::vector<double> SquareMatrix::operator*(std::vector<double> vVec) const{

	std::vector<double> mVec(dimension, 0);

	if (dimension != vVec.size()){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return mVec;
	}

	//the end product
	SquareMatrix mProd(dimension);

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

SquareMatrix& SquareMatrix::operator*=(const SquareMatrix& matrix2){
	*this = *this * matrix2;
	return *this;
}


void SquareMatrix::consoleInput(){ //inputs for the nxn matrix.

	std::string acoeff;
	int rowCounter = 0;

	do {
		std::cout << "Enter coefficients of linear expression "
			<< rowCounter + 1 << ":" << std::endl;
		for (int j = 0; j < dimension; ++j){
			std::cin >> acoeff;
			helperFunc(acoeff);
			matrix[rowCounter][j] = stod(acoeff);
		}
		++rowCounter;
	} while (rowCounter < dimension);

}

//inputs for the nxn matrix.
void SquareMatrix::fileInput(const std::string& filename){

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

void SquareMatrix::show() const{
	std::cout << "The matrix is " << std::endl;
	for (int i = 0; i < dimension; ++i){ //output of matrix
		for (int j = 0; j < dimension; ++j){
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	return;
}

int SquareMatrix::getDimension() const{
	return dimension;
}

SquareMatrix SquareMatrix::transpose() const{

	SquareMatrix trans(dimension);
	for (int i = 0; i < dimension; ++i){
		for (int j = 0; j < dimension; ++j){
			trans.matrix[j][i] = matrix[i][j];
		}
	}

	return trans;

}

SquareMatrix& 
	SquareMatrix::rowAdd(int rownum1, int rownum2, double multiplier){

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

SquareMatrix& 
	SquareMatrix::rowNeg(int rownum1, int rownum2, double multiplier){

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

SquareMatrix& SquareMatrix::rowMul(int rownum1, double multiplier){

	if (rownum1 < 1 || dimension < rownum1){
		std::cout << "Error: Parameter out of bounds" << std::endl;
		return *this;
	}

	for (int i = 0; i < dimension; ++i){
		matrix[rownum1 - 1][i] *= multiplier;
	}

	return *this;

}

SquareMatrix& SquareMatrix::rowSwap(int rownum1, int rownum2){

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
		return *this;
	}

	std::vector<double> temp(dimension, 0);

	for (int i = 0; i < dimension; ++i){
		temp[i] = matrix[rownum1 - 1][i];
		matrix[rownum1 - 1][i] = matrix[rownum2 - 1][i];
		matrix[rownum2 - 1][i] = temp[i];
	}

	return *this;

}

SquareMatrix SquareMatrix::submatrix(int row, int column) const{

	if (dimension == 1){
		SquareMatrix temp;
		return temp;
	}

	if (row < 1 || dimension < row){
		std::cout << "Index " << row << " is out of bounds; bound is "
			<< dimension << std::endl;
		return *this;
	}
	else if (column < 1 || dimension < column){
		std::cout << "Index " << column << " is out of bounds; bound is "
			<< dimension << std::endl;
		return *this;
	}

	SquareMatrix subM(dimension - 1);

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

std::vector<double> SquareMatrix::translation(double x, double y,
	const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	SquareMatrix mNew(3);

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

std::vector<double> SquareMatrix::horizontalShear(double multiplier,
	const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	SquareMatrix mNew(3);

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

std::vector<double> SquareMatrix::verticalShear(double multiplier,
	const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	SquareMatrix mNew(3);

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

std::vector<double> SquareMatrix::horizontalFlip(const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	SquareMatrix mNew(3);

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

std::vector<double> SquareMatrix::verticalFlip(const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	SquareMatrix mNew(3);

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

std::vector<double> SquareMatrix::squeezeMap(double multiplier,
	const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}
	else if (multiplier == 0){
		std::cout << "Error: Multiplier must be non-zero" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	SquareMatrix mNew(3);

	mNew.matrix[0][0] = multiplier;
	mNew.matrix[0][1] = 0;
	mNew.matrix[0][2] = 0;
	mNew.matrix[1][0] = 0;
	mNew.matrix[1][1] = 1/multiplier;
	mNew.matrix[1][2] = 0;
	mNew.matrix[2][0] = 0;
	mNew.matrix[2][1] = 0;
	mNew.matrix[2][2] = 1;

	return mNew*myVector;

}

std::vector<double> SquareMatrix::rotation(double angle,
	const std::vector<double>& myVector){

	if (myVector.size() != 3){
		std::cout << "Error: Must be a 3-dimensional vector" << std::endl;
		return std::vector < double > {0, 0, 0};
	}

	angle = fmod(angle, 360);

	SquareMatrix mNew(3);

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

SquareMatrix SquareMatrix::rowEchelon(){

	SquareMatrix temp = *this;

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

	changeToZerosIfNeeded(temp);

	if (temp.matrix[dimension - 1][dimension - 1] != 0)
		temp.matrix[dimension - 1][dimension - 1]
			/= temp.matrix[dimension - 1][dimension - 1];

	for (int i = 0; i < temp.dimension; ++i){
		for (int j = 0; j < dimension; j++){
			if (temp.matrix[i][j] != 1){
				continue;
			}
			else if (i != j){
				temp.rowSwap(i + 1, j + 1);
			}
		}
	}

	return temp;

}

SquareMatrix SquareMatrix::rowReducedEchelon(){

	SquareMatrix temp = rowEchelon();

	//back substitution

	for (int i = temp.dimension - 1; 1 <= i; --i){
		for (int j = i - 1; 0 <= j; --j){
			if (temp.matrix[i][i] != 0 && temp.matrix[j][i] != 0){
				temp.rowNeg(j + 1, i + 1, temp.matrix[j][i]);
			}
		}
	}
	

	return temp;
	
}

double SquareMatrix::getDeterminant(){

	if (dimension == 1){
		return matrix[0][0];
	}
	else if (dimension == 2){
		return matrix[0][0] * matrix[1][1]
			- matrix[0][1] * matrix[1][0];
	}

	double determinant = 0;

	for (int i = 0; i < dimension; ++i){
		determinant += matrix[0][i] * submatrix(1, i + 1).getDeterminant()
			* pow(-1, i);
	}

	return determinant;

}

bool SquareMatrix::isLinearlyIndependent(){

	SquareMatrix temp = rowReducedEchelon();
	for (int i = 0; i < temp.dimension; ++i){
		if (temp.matrix[i][i] != 1){
			return false;
		}
	}

	return true;

}

SquareMatrix SquareMatrix::inverse(){

	if (getDeterminant() == 0){
		std::cout << "Determinant is zero ---> Inverse does not exist"
			<< std::endl;
		return *this;
	}
	else if (dimension == 1){
		SquareMatrix temp;
		temp(1, 1) = 1 / matrix[0][0];
		return temp;
	}

	//calculate minors
	//also calculate cofactors
	SquareMatrix minors(dimension);

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
	return minors.transpose()*(1/getDeterminant());

}

std::vector<double> SquareMatrix::xVec(std::vector<double> bVec){

	SquareMatrix temp = *this;
	std::vector<double> xVector(temp.dimension, 0);

	//A^-1*b
	if (getDeterminant() != 0){
		SquareMatrix temp = *this;
		return temp.inverse()*bVec;
	}
	else if (temp.dimension != bVec.size()){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return xVector;
	}
	else if (bVec.size() == 1 && bVec[0] == 0){
		xVector.resize(1);
		if (temp(1, 1) == 0){
			std::cout << "Any solution works." << std::endl;
		}
		xVector[0] = 0;
		return xVector;
	}

	//row reduces the matrix temp.
	for (int i = 0; i < temp.dimension - 1; ++i){

		//swaps the rows if pivot point is zero
		if (temp.matrix[i][i] == 0){
			for (int k = i + 1; k < temp.dimension; ++k){
				if (temp.matrix[k][i] != 0){
					temp.rowSwap(k + 1, i + 1);
					std::swap(bVec[k], bVec[i]);
					break;
				}
			}
		}

		//multiplies entire row by reciprocal of pivot.
		for (int j = i; j < temp.dimension; ++j){
			if (temp.matrix[j][i] != 0){
				//bVec[j] reassignment MUST go before
				//temp reassignment.
				bVec[j] /= temp.matrix[j][i];
				temp.rowMul(j + 1, 1 / temp.matrix[j][i]);
			}
		}

		//subtracts rows from each other
		for (int j = i + 1; j < temp.dimension; ++j){
			if (temp.matrix[j][i] != 0){
				temp.rowNeg(j + 1, i + 1);
				bVec[j] -= bVec[i];
			}
		}

	}

	changeToZerosIfNeeded(temp);

	if (temp.matrix[dimension - 1][dimension - 1] != 0){

		bVec[dimension - 1] /= temp.matrix[dimension - 1][dimension - 1];
		temp.matrix[dimension - 1][dimension - 1] 
			/= temp.matrix[dimension - 1][dimension - 1];

	}

	//(inefficient) bubble-sort algorithm to 
	//place the pivots where they belong
	for (int i = 0; i < temp.dimension; ++i){
		for (int j = 0; j < dimension; j++){
			if (temp.matrix[i][j] != 1){
				continue;
			}
			else if (i != j){
				temp.rowSwap(i + 1, j + 1);
			}
		}
	}

	//matrix A is row-reduced and b is consistent with these operations
	//now we are ready for back substitution

	//the last index
	int start = temp.dimension - 1;

	for (int k = 0; k < temp.dimension; ++k){

		double restOfRowValue = 0;
		int counter = 0;
		//detects if equation has a free variable.
		if (bVec[start - k] == 0){
			for (int i = 0; i < temp.dimension; ++i){
				if (temp.matrix[start - k][i] == 0){
					++counter;
				}
				else
					break;
			}
			if (counter == temp.dimension){
				std::cout << "x_" << start - k + 1 << " is free;";
				std::cout << " setting it equal to one..." << std::endl;
				xVector[start - k] = 1;
				continue;
			}
		}

		//detects if equation is inconsistent
		else if (temp.matrix[start - k][start - k] == 0){
			std::cout << "Inconsistent system detected;";
			std::cout << " matrix has only trivial solution." << std::endl;
			std::vector<double> newVec(temp.dimension, 0);
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

std::vector<double> SquareMatrix::nullSpace(){
	std::vector<double> nullVec(dimension, 0);
	return xVec(nullVec);
}

/*std::vector<SquareMatrix> SquareMatrix::factorLU(){

	SquareMatrix L(dimension);
	SquareMatrix U(dimension);

	for (int i = 0; i < dimension; ++i){
		L.matrix[i][i] = 1;
	}

	std::vector<SquareMatrix> mArray;

	mArray.push_back(L);
	mArray.push_back(U);

	return mArray;

}*/
