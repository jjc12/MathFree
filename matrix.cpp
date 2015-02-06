// matrix.cpp: Source file and its contents.

#include "matrix.h"

void Matrix::helperFunc(std::string &strnum){
	std::string tryagain;
	int bad = 0;
	int watchout = 0;
	int strsize = strnum.size();

	for (int i = 0; i < strsize; i++){
		if (!((isdigit(strnum[i])) || (strnum[i] == '.'))){
			bad = 1;
			if ((strnum[0] == '-') && (1 < strsize)){
				for (int j = 1; j < strsize; j++){
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
			watchout++;
			if (watchout == 2){
				break;
			}
		}
	}

	if ((bad == 1) || (2 == watchout)){
		std::cout << "Invalid input; please enter numbers. Try again: ";
		std::cin >> tryagain;
		helperFunc(tryagain); //indeed, what if I am reading from a file or webpage? Ask Dr. Myers about this.
		//"make it a boolean and if it's false, end the program right then and there" is another option
		strnum = tryagain;
	}
	return;
}

Matrix::Matrix(int size){ //constructs an nxn matrix.
	if (size <= 0){
		std::cout << "Number of unknowns must be at least 1; setting unknowns to 1..." << std::endl;
		size = 1;
	}
	dimension = size;
	std::string acoeff;
	int rowCounter = 0;
	do {
		std::vector<double> row;
		std::cout << "Enter coefficients of linear equation " << rowCounter + 1 << " (including the coefficient in the b vector)." << std::endl;
		for (int j = 0; j < dimension; j++){
			std::cin >> acoeff;
			helperFunc(acoeff);
			row.push_back(stod(acoeff));
		}
		matrix.push_back(row);
		rowCounter++;
	} while (rowCounter < dimension);

	std::cout << "The matrix is " << std::endl;
	for (int i = 0; i < dimension; i++){ //output of matrix "in question"
		for (int j = 0; j < dimension; j++){
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

Matrix Matrix::Multiply(Matrix matrix2){
	//the end product
	std::vector< std::vector<double> > mProd(dimension, std::vector<double>(dimension));

	for (int i = 0; i < dimension; i++){ //performing the matrix multiplication.
		for (int j = 0; j < dimension; j++){
			for (int k = 0; k < dimension; k++){
				mProd[i][j] += matrix[i][k] * matrix2.getValue(j, k);
			}
		}
	}

	return //mProd;
}

double Matrix::getValue(int i, int j) const{
	if (dimension <= i){
		std::cout << "Index " << i << " is out of bounds; bound is " << dimension << std::endl;
	}
	else if (dimension <= j){
		std::cout << "Index " << j << " is out of bounds; bound is " << dimension << std::endl;
	}
	return this[i][j];
}

void Matrix::setValue(int i, int j, double value){
	this[i][j] = value;
}