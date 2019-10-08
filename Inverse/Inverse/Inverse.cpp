/*

обратное преобразование

*/

#include <iostream>
#include <fstream>
#include <cmath>

#define M_PI 3.14159265358979323846

using namespace std;

static int N = 200; // количество отсчётов


int main() {
	double Re, Im;
	double sumOfRealPart = 0.0, sumOfImaginaryPart = 0.0;
	double arrayOfFuncResult[200];
	double t = 0.0;
	double argumentOfTrigFunc, argument;
	double arrayOfComplexNumbers[101][2];

	double resOfInputFuncAfterInverse[200];
	double ReInv, ImInv, Ck;
	double sumOfRealPartInv = 0.0, sumOfImaginaryPartInv = 0.0, summaCk = 0.0;
	double mag[101];
	double phi[101];
	double stepOfFreq;
	double tmp;

	stepOfFreq = 1.0 / (N + 2);

	//TO DO:

	fstream file;

	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfInputSignal.txt");
	for (int i = 0; i < N; i++) {
		file >> tmp >> arrayOfFuncResult[i];
	}
	file.close();

	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfOutputSignalRe.txt");
	for (int i = 0; i <= N/2; i++) {
		file >> tmp >> arrayOfComplexNumbers[i][0];
	}
	file.close();

	t = 0.0;

	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfOutputSignalIm.txt");
	for (int i = 0; i <= N/2; i++) {
		file >> tmp >> arrayOfComplexNumbers[i][1];
	}
	file.close();

	arrayOfComplexNumbers[0][0] = arrayOfComplexNumbers[0][0] / 2;
	arrayOfComplexNumbers[N / 2][0] = arrayOfComplexNumbers[N / 2][0] / 2;
	for (int n = 0; n < N-1; n++) {
		for (int k = 0; k <= N/2; k++) {

			argumentOfTrigFunc = 2 * M_PI * k * n / N;
			ReInv = 2 * arrayOfComplexNumbers[k][0] * cos(argumentOfTrigFunc) / N; //Ak
			ImInv = -2 * arrayOfComplexNumbers[k][1] * sin(argumentOfTrigFunc) / N;//Bk
			Ck = sqrt(ReInv*ReInv + ImInv*ImInv);
			phi[k] = atan2(arrayOfComplexNumbers[k][1], arrayOfComplexNumbers[k][0]);
			argument = (2*M_PI*n*(k+phi[k])) / N;
			sumOfRealPartInv += ReInv;
			sumOfImaginaryPartInv += ImInv;
			summaCk += Ck;
		

		}
		resOfInputFuncAfterInverse[n] = Ck*cos(argument)/*(sumOfRealPartInv + sumOfImaginaryPartInv)*/;
		sumOfRealPartInv = 0.0;
		sumOfImaginaryPartInv = 0.0;
		summaCk = 0.0;
	}

	file.open("C:\\Users\\Олег\\Desktop\\2019\\SuperResult.txt");
	for (int i = 0; i < N; i++) {
		file << arrayOfFuncResult[i] << "\t" << resOfInputFuncAfterInverse[i] << "\n";
	}
	file.close();

	return 0;
}




///*
//
//
//обратное преобразование
//
//*/
//
//#include <iostream>
//#include <fstream>
//#include <cmath>
//
//#define M_PI 3.14159265358979323846
//
//using namespace std;
//
//static int N = 200; // количество отсчётов
//
//
//int main() {
//	double Re, Im;
//	double sumOfRealPart = 0.0, sumOfImaginaryPart = 0.0;
//	double arrayOfFuncResult[200];
//	double t = 0.0;
//	double argumentOfTrigFunc;
//	double arrayOfComplexNumbers[101][2];
//
//	double resOfInputFuncAfterInverse[200];
//	double ReInv, ImInv;
//	double sumOfRealPartInv = 0.0, sumOfImaginaryPartInv = 0.0;
//	double mag[101];
//	double phase[101];
//	double stepOfFreq;
//	double tmp;
//
//	stepOfFreq = 1.0 / (N + 2);
//
//	//TO DO:
//
//	fstream file;
//
//	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfInputSignal.txt");
//	for (int i = 0; i < N; i++) {
//		file >> tmp >> arrayOfFuncResult[i];
//	}
//	file.close();
//
//	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfOutputSignalRe.txt");
//	for (int i = 0; i <= N / 2; i++) {
//		file >> tmp >> arrayOfComplexNumbers[i][0];
//	}
//	file.close();
//
//	t = 0.0;
//
//	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfOutputSignalIm.txt");
//	for (int i = 0; i <= N / 2; i++) {
//		file >> tmp >> arrayOfComplexNumbers[i][1];
//	}
//	file.close();
//
//	arrayOfComplexNumbers[0][0] = arrayOfComplexNumbers[0][0] / 2;
//	arrayOfComplexNumbers[N / 2][0] = arrayOfComplexNumbers[N / 2][0] / 2;
//	for (int n = 0; n < N; n++) {
//		for (int k = 0; k <= N / 2; k++) {
//
//			argumentOfTrigFunc = 2 * M_PI * k * n / N;
//			ReInv = 2 * arrayOfComplexNumbers[k][0] * cos(argumentOfTrigFunc) / N; //Ak
//			ImInv = -2 * arrayOfComplexNumbers[k][1] * sin(argumentOfTrigFunc) / N; //Bk
//			sumOfRealPartInv += ReInv;
//			sumOfImaginaryPartInv += ImInv;
//		}
//		resOfInputFuncAfterInverse[n] = (sumOfRealPartInv + sumOfImaginaryPartInv);
//		sumOfRealPartInv = 0.0;
//		sumOfImaginaryPartInv = 0.0;
//	}
//
//	file.open("C:\\Users\\Олег\\Desktop\\2019\\SuperResult.txt");
//	for (int i = 0; i < N; i++) {
//		file << arrayOfFuncResult[i] << "\t" << resOfInputFuncAfterInverse[i] << "\n";
//	}
//	file.close();
//
//	return 0;
//}