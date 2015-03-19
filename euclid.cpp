#include "euclid.h"

int Euclid::gcd(int num1, int num2){
	if (num2 == 0){
		return num1;
	}
	return gcd(num2, num1%num2);
}

bool Euclid::isPrime(int num1){
	if ((num1 == 0) || num1 == 1){
		return false;
	}
	if ((num1 == 2) || num1 == 3){
		return true;
	}
	std::vector<int> primes(num1+1);
	for (int i = 0; i <= num1; i++){
		primes[i] = i;
	}
	for (int i = 2; i < sqrt(num1); i++){
		for (int multiplier = 2; i*multiplier <= num1; multiplier++){
			primes[i*multiplier] = 0;
		}
	}
	primes[1] = 0;
	int linecounter = 0;
	for (int k = 0; k <= num1; k++){
		if (num1 == primes[num1]){
			std::cout << num1 << " is prime.";
			return true;
		}
	}
	return false;
}

void Euclid::listOfPrimes(int num1){ //lists primes up to a given num1
	std::vector<int> primes(num1 + 1);
	for (int i = 0; i <= num1; i++){
		primes[i] = i;
	}
	for (int i = 2; i < sqrt(num1); i++){
		if (primes[i] != 0){
			for (int multiplier = 2; i*multiplier <= num1; multiplier++){
				primes[i*multiplier] = 0;
			}
		}
		else
			continue;
	}
	primes[1] = 0;
	int linecounter = 0;
	for (int i = 0; i <= num1; i++){
		if (primes[i] != 0){
			std::cout << i << " ";
			linecounter++;
			if (linecounter % 8 == 0)
				std::cout << std::endl;
		}
	}
}

int Euclid::lcm(int num1, int num2){
	int lcmnum = (num1*num2) / gcd(num1, num2);
	return lcmnum;
}

void Euclid::linDio(int num1, int num2, int num3){
	//an extended Euclidean algorithm.
	int a = num1;
	int b = num2;
	if (a < b){ //a must be greater than b in our algorithm, so we sort it should we need to.
		int temp = a;
		a = b;
		b = temp;
	}
	/*if (num3 < 0){
		a *= -1;
		b *= -1;
		num3 *= -1;
	}*/
	int gcdivisor = gcd(a, b);
	std::cout << "Finding solutions to " << a << "x + " << b << "y = " << num3 << "..." << std::endl;
	int acopy = a; //for later
	int bcopy = b; //for later
	if (num3%gcdivisor != 0){
		std::cout << "No solutions exist." << std::endl;
		std::cout << "The gcd of num1 and num2, " << gcdivisor << ", does not divide " << num3 << "." << std::endl;
		return;
	}
	if ((num3 == 0) && (num1 == 0) && (num2 == 0)){
		std::cout << "Any integer x, y will satisfy this equation." << std::endl;
		return;
	}
	if ((num3 == 0) && ((num1 != 0) && (num2 != 0))){
		std::cout << "A general solution to the system is" << std::endl;
		std::cout << "x = " << b << " + " << num2 / gcdivisor << "t" << std::endl;
		std::cout << "y = " << -a << " - " << num1 / gcdivisor << "t" << std::endl;
		return;
	}
	if (a == 0){
		std::cout << "The only solution is " << std::endl;
		std::cout << "y = " << num3 / b << std::endl;
		return;
	}
	if (b == 0){
		std::cout << "The only solution is " << std::endl;
		std::cout << "x = " << num3 / a << std::endl;
		return;
	}

	std::vector<int> rVec; //array of remainders
	std::vector<int> qVec; //array of q's
	std::vector<int> xVec;
	std::vector<int> yVec;
	int remainder = 0;
	rVec.push_back(0); //for convenience
	qVec.push_back(0); //for convenience
	xVec.push_back(0); //for convenience
	yVec.push_back(0); //for convenience
	int switcher = 0;

	for (int i = 1; remainder != gcdivisor; i++, switcher++){
		if (i != 1){
			a = b;
			b = rVec[i - 1];
		}
		int tester = a / b; //finds the floor of a/b, i.e. q, so that 0 <= r < b
		qVec.push_back(tester); //q_n = tester
		remainder = a - b*qVec[i]; //a - b*q_1 = r_1, but remainder = a%b*q_i works as well.
		if (remainder == 0){ //so that r_1 is zero.
			std::cout << "A general solution to the system is " << std::endl;
			std::cout << "x = " << 0 << " + " << b/gcdivisor << "t" << std::endl;
			std::cout << "y = " << num3/b << " - " << a/gcdivisor << "t" << std::endl;
			return;
		}
		rVec.push_back(remainder); //r_n = remainder
	}
	if (switcher == 1){ //r_2 is zero
		std::cout << "A general solution to the system is" << std::endl;
		std::cout << "x = " << 1 * (num3 / gcdivisor) << " + " << bcopy / gcdivisor << "t" << std::endl;
		std::cout << "y = " << (-1 * qVec[1])*(num3 / gcdivisor) << " - " << acopy / gcdivisor << "t" << std::endl;
		return;
	}
	if (switcher == 2){ //r_3 is zero
		std::cout << "A general solution to the system is" << std::endl;
		std::cout << "x = " << (-1 * qVec[2])*(num3 / gcdivisor) << " + " << bcopy / gcdivisor << "t" << std::endl;
		std::cout << "y = " << (1 + qVec[1] * qVec[2])*(num3 / gcdivisor) << " - " << acopy / gcdivisor << "t" << std::endl;
		return;
	}

	//now we have our values and we perform a recursive method to find a solution to the system.
	xVec.push_back(1); //x_1 = 1 always
	yVec.push_back(-1 * qVec[1]); //y_1 = -1*q_1 always
	xVec.push_back(-1 * qVec[2]); //x_2 = -1*q_2 always
	yVec.push_back(1 + qVec[1] * qVec[2]); //y_2 = 1 + q_1*q_2 always
	int rVecSize = rVec.size();
	for (int i = 3; i < rVecSize; i++){
		xVec.push_back(xVec[i - 2] - xVec[i - 1] * qVec[i]);
		yVec.push_back(yVec[i - 2] - yVec[i - 1] * qVec[i]);
	}
	std::cout << "A general solution to the system is" << std::endl;
	std::cout << "x = " << xVec[rVecSize - 1] * (num3 / gcdivisor) << " + " << bcopy / gcdivisor << "t" << std::endl;
	std::cout << "y = " << yVec[rVecSize - 1] * (num3 / gcdivisor) << " - " << acopy / gcdivisor << "t" << std::endl;
}