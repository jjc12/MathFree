// matrix.cpp: Source file and its contents.

#include "augmatrix.h"

void AugMatrix::helperFunc(std::string &strnum){
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

AugMatrix::AugMatrix(int numx){
	if (numx <= 0){
		std::cout << "Number of unknowns must be at least 1; setting unknowns to 1..." << std::endl;
		numx = 1;
	}
	numberOfx = numx;
	std::string acoeff;
	int rowCounter = 0;
	do {
		std::vector<double> row;
		std::cout << "Enter coefficients of linear equation " << rowCounter + 1 << " (including the coefficient in the b vector)." << std::endl;
		for (int j = 0; j <= numberOfx; j++){
			std::cin >> acoeff;
			helperFunc(acoeff);
			row.push_back(stod(acoeff));
		}
		matrix.push_back(row);
		rowCounter++;
	} while (rowCounter < numberOfx);

	for (int i = 0; i < numberOfx; i++){
		xArray.push_back(0);
	}
}

void AugMatrix::findX(){
	for (int start = 0; !(start == numberOfx - 1); start++){
		/*std::cout << "The matrix I am about to solve is " << std::endl;
		for (int i = start; i < numberOfx; i++){ //output of matrix "in question"
		for (int j = start; j < numberOfx + 1; j++){
		std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
		}
		std::cout << std::endl;*/
		double maxValue = abs(matrix[start][start]);
		int maxRow = start;
		for (int i = start; i < numberOfx; i++){ //finding the max (absolute) value in column "start" in preparation for 1st step of elimination.
			if (maxValue < abs(matrix[i][start])){
				maxValue = abs(matrix[i][start]);
				maxRow = i;
			}
		}
		/*std::cout << "Final (absolute) maxValue is " << maxValue << " at row " << maxRow << std::endl;
		std::cout << std::endl;
		std::cout << "Now swapping..." << std::endl;
		std::cout << std::endl;*/
		for (int i = start; i < numberOfx + 1; i++){ //swapping of the rows
			iter_swap(matrix[start].begin() + i, matrix[maxRow].begin() + i);
		}
		/*std::cout << "New matrix is" << std::endl;
		for (int i = start; i < numberOfx; i++){ //output of matrix with swapped rows, ready for operations.
		for (int j = start; j < numberOfx + 1; j++){
		std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
		}*/
		double bufferValue;
		for (int i = start; i < numberOfx; i++){ // the first step of the algorithm in order to begin Gaussian elimination. HERE is where we have the bug "division by zero exception"
			bufferValue = matrix[i][start];
			for (int j = start; j < numberOfx + 1; j++){
				if (bufferValue == 0)
					continue; // avoids division by zero
				matrix[i][j] /= (bufferValue / maxValue);
			}
		}
		/*std::cout << std::endl;
		std::cout << "New matrix is" << std::endl;
		for (int i = start; i < numberOfx; i++){ //output of matrix after preparation for Gaussian elimination.
		for (int j = start; j < numberOfx + 1; j++){
		std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
		}
		std::cout << std::endl;*/
		for (int i = start + 1; i < numberOfx; i++){ //now we subtract r_i-r_j;
			if (matrix[i][start] == 0)
				continue; // no need to subtract when the a[i][start] is a zero.
			for (int j = start; j < numberOfx + 1; j++){
				matrix[i][j] -= matrix[start][j];
			}
		}
		/*for (int i = 0; i < numberOfx; i++){ //output of matrix after row subtraction.
		for (int j = 0; j < numberOfx + 1; j++){
		std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
		}
		std::cout << std::endl;*/
	}
	int limit = numberOfx;
	double numerator;
	double denominator;
	for (int k = 0; k < limit; k++){
		numerator = matrix[numberOfx - 1 - k][numberOfx];
		int indexStart = numberOfx;
		int indexLimit = numberOfx - k;
		while (indexLimit < indexStart){
			numerator -= matrix[numberOfx - 1 - k][indexStart - 1] * xArray[indexStart - 1];
			indexStart--;
		}
		denominator = matrix[numberOfx - 1 - k][numberOfx - 1 - k];
		if ((numerator == 0) && (denominator == 0)){
			std::cout << "x_" << numberOfx - k << " is a free variable. Setting it equal to 0..." << std::endl;
			xArray[numberOfx - 1 - k] = 0;
			continue;
		}
		else if (denominator == 0){
			std::cout << "This system of linear equations has no solution other than the trivial solution." << std::endl;
			int xSize = xArray.size();
			for (int r = 0; r < xSize; r++){
				xArray[r] = 0;
			}
			return;
		}
		else
			xArray[numberOfx - 1 - k] = numerator / denominator;
	}
	std::vector<double> rArray(numberOfx);
	int xSize = xArray.size();

	std::cout << "The solution to the matrix is " << std::endl;
	for (int i = 0; i < xSize; i++){
		std::cout << "x_" << i + 1 << " = " << xArray[i] << std::endl;
	}

	for (int i = 0; i < numberOfx; i++){ //creation of residuals
		rArray[i] = matrix[i][numberOfx];
		for (int j = 0; j < numberOfx; j++){
			rArray[i] -= matrix[i][j] * xArray[j];
		}
	}
	std::cout << std::endl;
	std::cout << "The residual for this matrix is " << std::endl;
	for (int i = 0; i < xSize; i++){
		std::cout << "r_" << i + 1 << " = " << rArray[i] << std::endl;
	}
}