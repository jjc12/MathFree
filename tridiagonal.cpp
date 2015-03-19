// matrix.cpp: Source file and its contents.

#include "tridiagonal.h"

void TriDiagonal::helperFunc(std::string &strnum){
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

TriDiagonal::TriDiagonal(int numx){
	if (numx <= 1){
		std::cout << "Number of unknowns must be at least 2; setting unknowns to 2..." << std::endl;
		numx = 2;
	}
	numberOfx = numx;
	double acoeff;
	int rowCounter = 0;
	do {
		std::vector<double> row;
		for (int j = 0; j < numberOfx; j++){ //initialization of coefficient matrix to zero.
			acoeff = 0;
			row.push_back(acoeff);
		}
		matrix.push_back(row);
		rowCounter++;
	} while (rowCounter < 3);

	std::string tcoeff; //tester coefficient

	std::cout << "Enter coefficients of band 1 a[i][i-1], i.e. the band below the main diagonal." << std::endl;
	for (int i = 1; i < numberOfx; i++){ //band 1
		std::cin >> tcoeff;
		helperFunc(tcoeff);
		matrix[0][i] = stod(tcoeff);
	}

	std::cout << "Enter coefficients of band 2 a[i][i], i.e. the band along the main diagonal." << std::endl;
	for (int i = 0; i < numberOfx; i++){ //band 2
		std::cin >> tcoeff;
		helperFunc(tcoeff);
		matrix[1][i] = stod(tcoeff);
	}

	std::cout << "Enter coefficients of band 3 a[i][i+1], i.e. the band above the main diagonal." << std::endl;
	for (int i = 0; i < numberOfx - 1; i++){ //band 3
		std::cin >> tcoeff;
		helperFunc(tcoeff);
		matrix[2][i] = stod(tcoeff);
	}

	for (int i = 0; i < numberOfx; i++){ //b-vector initialization
		std::cout << "Enter b-vector coefficient " << i + 1 << ":" << std::endl;
		std::cin >> tcoeff;
		helperFunc(tcoeff);
		bVector.push_back(stod(tcoeff));
	}
}

void TriDiagonal::SolveJacobi(){  
	double acoeff;

	int rowCounter = 0;
	do {
		std::vector<double> row;
		for (int j = 0; j < numberOfx; j++){ //initialization of remainder matrix to zero. does not contain b-vector.
			acoeff = 0;
			row.push_back(acoeff);
		}
		remMatrix.push_back(row);
		rowCounter++;
	} while (rowCounter < 2); //each iter makes a new row for the matrix.

	for (int i = 0; i < numberOfx; i++){ //the inverse of the diagonal coefficients matrix (and initialization).
		diagMatrix.push_back(1 / matrix[1][i]);
	}

	for (int j = 1; j < numberOfx; j++){ //init of remMatrix to non-diag coeffs, first row
		remMatrix[0][j] = matrix[0][j];
	}

	for (int j = 0; j < numberOfx - 1; j++){ //init of remMatrix to non-diag coeffs, second row.
		remMatrix[1][j] = matrix[2][j];
	}

	/*std::cout << "The remainder matrix is" << std::endl;
	for (int i = 0; i <= 1; i++){
		for (int j = 0; j < numberOfx; j++){
			std::cout << remMatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "The inverse diagonal matrix is" << std::endl;
	for (int i = 0; i < numberOfx; i++){
		std::cout << diagMatrix[i] << std::endl;
	}*/

	double epsilon;
	std::string epsilonstr;
	std::cout << "What would you like the tolerance to be?" << std::endl;
	std::cin >> epsilonstr; // the distance between xIter elements between rows.
	helperFunc(epsilonstr);
	epsilon = fabs(stod(epsilonstr));

	std::vector< std::vector<double> > xIter(2, std::vector<double>(numberOfx, 0));

	for (int i = 0; i < numberOfx; i++){ //initialization of x^{(0)}.
		std::string estimates;
		std::cout << "Enter an estimate for x_" << i + 1 << ":" << std::endl;
		std::cin >> estimates;
		helperFunc(estimates);
		xIter[0][i] = stod(estimates);
		//the transpose of the solution, i.e. ([x_1, x_2, x_3, ...., x_n]^T)^(0)
	}

	int k;
	int switcher = 0;
	std::vector<double> AX(numberOfx, 0);

	for (k = 1; switcher != numberOfx; k++){
		std::vector<double> remTimesX(numberOfx, 0);
		for (int i = 0; i < numberOfx; i++){
			if (i != 0){
				remTimesX[i] += remMatrix[0][i] * xIter[0][i - 1];
			}
			if (i != numberOfx - 1){
				remTimesX[i] += remMatrix[1][i] * xIter[0][i + 1];
			}
		}

		for (int i = 0; i < numberOfx; i++){ //initialization of x^{(1)}.
			double nextIterElement = (diagMatrix[i])*(bVector[i] - remTimesX[i]);// INVERSE(D)*(b - Rx^(0));
			xIter[1][i] = nextIterElement; //the transpose of the solution, i.e. ([x_1, x_2, x_3, ...., x_n]^T)^(1)
		}

		for (int test = 0; test < numberOfx; test++){ // testing if the residual is less than epsilon.
			double thisCurrentElement = (1 / diagMatrix[test])*xIter[0][test] + remTimesX[test];
			AX[test] = thisCurrentElement;
			if (fabs(bVector[test] - AX[test]) <= epsilon){
				switcher++;
			}
			else{
				switcher = 0;
				break;
			}
		}

		for (int i = 0; i < numberOfx; i++){ //swap
			xIter[0][i] = xIter[1][i];
		}
	}

	std::cout << "The solution to this matrix after " << k - 1 << " iterations" << std::endl;
	std::cout << "with a tolerance of " << epsilon << " is" << std::endl;

	std::streamsize oldprecision = std::cout.precision();

	std::cout.precision(2);
	for (int r = 0; r < numberOfx; r++){
		std::cout << "x_" << r + 1 << " = " << xIter[0][r] << " and residual r_" << r + 1 << " = " << fabs(bVector[r] - AX[r]) << std::endl;
	}
	std::cout.precision(oldprecision);
	std::cout << "Done!" << std::endl;

}