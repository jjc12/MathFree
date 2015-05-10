// matrix.cpp: Source file and its contents.

#include "squarematrix_free.h"

std::ostream& operator<<(std::ostream& s, SquareMatrix m){

	std::cout << "[";

	//output of matrix
	for (int i = 0; i < m.dimension; ++i){
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

	//output of matrix
	for (unsigned i = 0; i < v.size(); ++i){
		std::cout << v[i];
		if (i != v.size() - 1)
			std::cout << ", ";
	}
	std::cout << "]";

	return s;

}

void SquareMatrix::isValidInput(std::string &strnum){
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
	//redirects you to the beginning of the function.
	if (bad == 1 || 2 <= watchout){
		std::cout << "Invalid input; please enter numbers. Try again: ";
		std::cin >> tryagain;
		isValidInput(tryagain); 
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

void SquareMatrix::changeToZerosIfNeeded(std::vector<double>& myV){

	//changes values if they are below a certain tolerance (due to
	//floating-point precision errors)
	for (unsigned i = 0; i < myV.size(); ++i){
		if (fabs(myV[i]) < tolerance){
			myV[i] = 0;
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
			isValidInput(acoeff);
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

std::vector<double>
SquareMatrix::horizontalFlip(const std::vector<double>& myVector){

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

std::vector<double>
SquareMatrix::verticalFlip(const std::vector<double>& myVector){

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

	//inefficient bubble-sort sorting algorithm to place all zero rows
	//at the bottom

	for (int i = 0; i < temp.dimension; ++i){
		for (int j = i; j < temp.dimension; ++j){

			//breaks if row j is not a zero vector.
			if (temp.matrix[i] != std::vector<double>(temp.matrix[i].size(), 0)
			){
				break;
			}

			//may swap zero vectors with themselves once all zero vectors are
			//placed at the bottom; not a problem, still.
			if (temp.matrix[i] == std::vector<double>(temp.matrix[i].size(), 0)
				&& temp.matrix[i] != temp.matrix[j]){
				temp.rowSwap(i + 1, j + 1);
			}

		}
	}

	//the matrix temp transforms to row-echelon form.
	for (int i = 0; i < temp.dimension - 1; ++i){

		//swaps the rows if pivot point is zero
		if (temp.matrix[i][i] == 0){
			int switch_ = 0;
			for (int k = i + 1; k < temp.dimension; ++k){
				if (temp.matrix[k][i] != 0){
					temp.rowSwap(k + 1, i + 1);
					switch_ = 1;
					break;
				}
			}
			if (switch_ == 0){

				//this means all numbers below and including pivot are zero.

				//preliminary check to see if columns to the right of pivot
				//(and below) are equal to zero 
				//(we must check, or else infinite loop of zero matrix occurs)
				bool switch_ = false;
				for (int row = 0; i + row < temp.dimension; ++row){
					if (temp.matrix[i + row]
						!= std::vector<double>(temp.matrix[i].size(), 0)){
						switch_ = true;
					}
				}
				if (switch_ == false){
					return temp;
				}

				//move to the next column, treat that element as pivot.
				SquareMatrix temp_(temp.dimension - i);

				//shift the column by 1 to the right, copy it to temp_.
				for (int row = 0; row < temp_.dimension; ++row){
					for (int column = 0; column < temp_.dimension - 1;
						++column){
						temp_.matrix[row][column]
							= temp.matrix[i + row][i + 1 + column];
					}
				}

				//perform row echelon on this modified temp matrix
				temp_ = temp_.rowEchelon();

				//set temp's elements equal to temp_'s by shifting
				//the column by 1 to the left.
				for (int row = 0; row < temp_.dimension; ++row){
					for (int column = 0; column < temp_.dimension - 1;
						++column){
						temp.matrix[i + row][i + 1 + column]
							= temp_.matrix[row][column];
					}
				}

			}
		}

		//multiplies entire row by reciprocal of pivot.
		for (int j = i; j < temp.dimension; ++j){
			if (temp.matrix[j][i] != 0){
				temp.rowMul(j + 1, 1 / temp.matrix[j][i]);
			}
		}
		//subtracts rows from each other
		//this loop removes all numbers below pivot at row j.
		for (int j = i + 1; j < temp.dimension; ++j){
			if (temp.matrix[j][i] != 0){
				temp.rowNeg(j + 1, i + 1);
			}
		}

	}

	changeToZerosIfNeeded(temp);

	if (temp.matrix[temp.dimension - 1][temp.dimension - 1] != 0)
		temp.matrix[temp.dimension - 1][temp.dimension - 1]
			/= temp.matrix[temp.dimension - 1][temp.dimension - 1];

	return temp;

}

SquareMatrix SquareMatrix::rowReducedEchelon(){

	//Get the matrix in row-echelon form
	SquareMatrix temp = rowEchelon();
	
	//preliminary check to see if top rows are zero rows
	//(we must check, or else infinite loop of zero matrix occurs)
	int zeroRows = 0;
	for (int row = dimension - 1; 0 <= row; --row){
		if (temp.matrix[row]
			== std::vector<double>(temp.matrix[row].size(), 0)){
			++zeroRows;
		}
	}
	//it's a zero matrix
	if (zeroRows == temp.dimension){
		return temp;
	}

	//create zeros above each pivot
	for (int k = temp.dimension - 1 - zeroRows; 1 <= k; --k){
		int zerosBeforePivot = 0;
		//checks where the pivot is
		for (int j = 0; j < temp.dimension - 1; ++j){
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
				temp.rowNeg(l + 1, k + 1, temp.matrix[l][j]);
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

	std::vector<double> xVector(dimension, 0);

	//A^-1*b
	if (getDeterminant() != 0){
		return (*this).inverse()*bVec;
	}
	//incompatible dimensions
	else if (dimension != bVec.size()){
		std::cout << "Error: incompatible dimensions" << std::endl;
		return xVector;
	}
	//cx = 0
	else if (bVec.size() == 1 && bVec[0] == 0){
		xVector.resize(1);
		//0x = 0
		if (matrix[0][0] == 0){
			std::cout << "Any solution works." << std::endl;
		}
		//c*0 = 0
		xVector[0] = 0;
		return xVector;
	}

	//Assume the matrix is a (n+1)x(n+1) with a zero row at the bottom.
	//copying the nxn portion
	SquareMatrix temp(this->dimension + 1);
	for (int i = 0; i < dimension; ++i){
		for (int j = 0; j < dimension; ++j){
			temp.matrix[i][j] = this->matrix[i][j];
		}
	}
	//copying the nx1th column
	for (int i = 0; i < temp.dimension - 1; ++i){
		temp.matrix[i][temp.dimension - 1] = bVec[i];
	}

	temp = temp.rowEchelon();

	//Augmented matrix A is in row-echelon
	//now we are ready for back substitution

	//the last index
	int start = temp.dimension - 2;

	//detect if A is consistent by making sure no pivot exists in b-column.
	std::vector<double> zero_vector(temp.dimension - 1, 0);
	for (int i = 0; i < temp.dimension - 1; ++i){
		std::vector<double> test_vector(temp.dimension - 1, 0);
		for (int j = 0; j < temp.dimension - 1; ++j){
			test_vector[j] = temp.matrix[i][j];
		}
		if (test_vector == zero_vector
			&& temp.matrix[i][temp.dimension - 1] != 0){
			std::cout << "Inconsistent system detected;";
			std::cout << " matrix has no solution. " << std::endl;
			std::cout << "Returning default zero vector." << std::endl;
			std::vector<double> newVec(temp.dimension - 1, 0);
			return newVec;
		}
	}

	//detect which variables of A are free
	//by checking which columns aren't pivot columns
	std::vector<bool> pivotColumns(temp.dimension - 1, false);
	for (int i = 0; i < temp.dimension - 1; ++i){
		for (int j = 0; j < temp.dimension - 1; ++j){
			if (temp.matrix[i][j] != 0){
				pivotColumns[j] = true;
				break;
			}
		}
	}

	//runs through pivotColumns; false means the variable is free
	for (int i = 0; i < temp.dimension - 1; ++i){
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
		if (temp.matrix[i] == std::vector<double>(temp.dimension, 0)){
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

		xVector[k] = temp.matrix[i][temp.dimension - 1] - _Ax;
	}

	changeToZerosIfNeeded(xVector);

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
