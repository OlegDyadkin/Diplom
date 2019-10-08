 /*

прямое преобразование

*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>

#define M_PI 3.14159265358979323846

using namespace std;

static int N = 200; // количество отсчётов
double LineInter(double w1, double w2, double c1, double c2, double R) {

	double res;

	res = c1 + (c2 - c1) / (w2 - w1)*(R - w1);



	return res;
}


int main() {
	double Re, Im;
	double sumOfRealPart = 0.0, sumOfImaginaryPart = 0.0;
	double arrayOfFuncResult[200];
	double t = 0.0;
	double argumentOfTrigFunc;
	double arrayOfComplexNumbers[101][2];
	double NaborCk[101];

	double resOfInputFuncAfterInverse[200];
	double ReInv, ImInv;
	double sumOfRealPartInv = 0.0, sumOfImaginaryPartInv = 0.0, summaCk = 0.0;
	double mag[101];
	double phase[101];
	double stepOfFreq;
	double tmp;
    double omega[15] = { 0, 77022.06, 115510.3815, 153971.203, 383590.5343, 755563.7828, 1089856.852, 1334024.614, 1493426.204, 1636126.952, 1789995.119, 1958672.888, 2331441.97, 2734717.393, 3155165.375 };
    double c[15] = { 4904.144463, 4903.376713, 4902.412853, 4901.055611, 4884.026372, 4810.068434, 4625.496139, 4246.332231, 3802.978599, 3471.969227, 3255.846529, 3117.32472, 2968.484112, 2901.625273,2869.486661 };
	double w1, w2, c1, c2;
    double Pn, Ck;
	double dT= 0.0000001 /** 10.0*/;
	double T0 =(N-1)*dT;
    double cres[101];
	double C[101];
	double phi[101];
    srand(time(NULL));

	stepOfFreq = 1.0 / (N * dT);

	//TO DO:
	fstream file;

	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfInputSignal.txt");
	for (int i = 0; i < N; i++) {
		file >> tmp >> arrayOfFuncResult[i];
	}
	file.close();

	t = 0.0;
	/*
	Прямое преобразование
	*/

	for (int k = 0; k <= N / 2; k++) {
		for (int n = 0; n < N; n++) {

			argumentOfTrigFunc = 2 * M_PI * k * n / N;
			Re = 2 * arrayOfFuncResult[n] * cos(argumentOfTrigFunc) / N;
			Im = 2 * arrayOfFuncResult[n] * sin(argumentOfTrigFunc) / N;
			Ck = sqrt(Re*Re + Im*Im);
			sumOfRealPart += Re;
			sumOfImaginaryPart += Im;
			summaCk += Ck;

		}
		NaborCk[k] =summaCk;
		arrayOfComplexNumbers[k][0] = sumOfRealPart; //Ak
		arrayOfComplexNumbers[k][1] = sumOfImaginaryPart; //Bk
		sumOfRealPart = 0.0;
		sumOfImaginaryPart = 0.0;
		summaCk = 0.0;

	}

	file.open("C:\\Users\\Олег\\Desktop\\2019\\Ck.txt");//Ck
	for (int i = 0; i <= N / 2; i++) {
		file << t << "\t" << NaborCk[i] << "\n";
		t += stepOfFreq;
	}
	file.close();


	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfOutputSignalRe.txt");//RE
	for (int i = 0; i <= N/2; i++) {
		file << t << "\t" << arrayOfComplexNumbers[i][0] << "\n";
		t += stepOfFreq;
	}
	file.close();

	t = 0.0;

	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfOutputSignalIm.txt");//Im
	for (int i = 0; i <= N/2; i++) {
		file << t << "\t" << arrayOfComplexNumbers[i][1] << "\n";
		t += stepOfFreq;
	}
	file.close();

	file.open("C:\\Users\\Олег\\Desktop\\2019\\PointsOfComplexCoord.txt"); //mag phase
	for (int k = 0; k <= N/2; k++) {
		mag[k] = sqrt(arrayOfComplexNumbers[k][0] * arrayOfComplexNumbers[k][0] + arrayOfComplexNumbers[k][1] * arrayOfComplexNumbers[k][1]);
		phase[k] = atan2(arrayOfComplexNumbers[k][1], arrayOfComplexNumbers[k][0]);

		double R = (((double)rand()) / RAND_MAX) * 3155165;

		for (int j = 0; j <= 13; j++) {
			if ((R > omega[j]) && (R < omega[j + 1])) {
				w1 = omega[j];
				w2 = omega[j + 1];
				c1 = c[j];
				c2 = c[j + 1];
				break;
			}
		}


	/*	printf("R=%lf", R);
		cout << endl;
		printf("w1=%lf", w1);
		cout << endl;
		printf("w2=%lf", w2);
		cout << endl;*/

		cres[k] = LineInter(w1, w2, c1, c2, R);
		/*cout << cres[k] << endl;*/
		file << mag[k] << "\t" << phase[k] << "\t" << R << "\t" << cres[k] << "\n";
	}
	file.close();

	//file.open("C:\\Users\\Олег\\Desktop\\2019\\PnNew.txt"); //C[k] phi
	//for (int j = 0; j <= N - 1; j++) {
	///*	C[j] = sqrt(arrayOfComplexNumbers[j][0] * arrayOfComplexNumbers[j][0] + arrayOfComplexNumbers[j][1] * arrayOfComplexNumbers[j][1]);
	//	phi[j] = atan2(arrayOfComplexNumbers[j][1], arrayOfComplexNumbers[j][0]);*/

	//	Pn = mag[k] * cos(2 * M_PI)*j*(n + phase[k]) / N;



	//	file <<  << "\t" <<  << "\t" << "\n";
	//}
	//file.close();


	//for (int k = 0; k <= N / 2; k++) {
	//	for (int n = 0; n < N; n++) {

	//		argumentOfTrigFunc = 2 * M_PI * k * n / N;
	//		Re = 2 * arrayOfFuncResult[n] * cos(argumentOfTrigFunc) / N;
	//		Im = 2 * arrayOfFuncResult[n] * sin(argumentOfTrigFunc) / N;
	//		Ck = sqrt(Re*Re + Im*Im);
	//		sumOfRealPart += Re;
	//		sumOfImaginaryPart += Im;
	//		summaCk += Ck;

	//	}
	//	arrayOfComplexNumbers[k][0] = sumOfRealPart; //Ak
	//	arrayOfComplexNumbers[k][1] = sumOfImaginaryPart; //Bk
	//	sumOfRealPart = 0.0;
	//	sumOfImaginaryPart = 0.0;
	//	summaCk = 0.0;

	//}





	file.open("C:\\Users\\Олег\\Desktop\\2019\\IsReAndImEqualAfterFunc.txt");
	for (int k = 0; k <= N / 2; k++) {
		if (mag[k] * cos(phase[k]) != arrayOfComplexNumbers[k][0])
			Re = -mag[k] * cos(phase[k]);
		else
			Re = mag[k] * cos(phase[k]);
		if (mag[k] * sin(phase[k]) != arrayOfComplexNumbers[k][1])
			Im = -mag[k] * sin(phase[k]);
		else
			Im = mag[k] * sin(phase[k]);
		file << Re << "\t" << arrayOfComplexNumbers[k][0] << "\t\t\t" << Im << "\t" << arrayOfComplexNumbers[k][1] << "\n";
	}
	file.close();

	return 0;
}


//#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <time.h>
//
//#define M_PI 3.14159265358979323846
//
//using namespace std;
//
//static int N = 200; // количество отсчётов
//
//double LineInter(double w1, double w2, double c1, double c2, double R) {
//
//	double res;
//
//	res = c1 + (c2 - c1) / (w2 - w1)*(R - w1);
//
//
//
//	return res;
//}
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
//	double omega[15] = { 0, 77022.06, 115510.3815, 153971.203, 383590.5343, 755563.7828, 1089856.852, 1334024.614, 1493426.204, 1636126.952, 1789995.119, 1958672.888, 2331441.97, 2734717.393, 3155165.375 };
//	double c[15] = { 4904.144463, 4903.376713, 4902.412853, 4901.055611, 4884.026372, 4810.068434, 4625.496139, 4246.332231, 3802.978599, 3471.969227, 3255.846529, 3117.32472, 2968.484112, 2901.625273,2869.486661 };
//	double w[512];
//	double w1;
//	double w2;
//	double c1;
//	double c2;
//	double cres[100];
//	srand(time(NULL));
//
//	
//	cout << "RAND_MAX = "<< RAND_MAX << endl;
//
//
//	for (int k=0; k <= N/2; k++) {
//		
//		double R = (((double)rand()) / RAND_MAX) * 3155165;
//		
//		for (int j = 0; j <= 13; j++) {
//			if ((R > omega[j]) && (R < omega[j + 1])) {
//				 w1 = omega[j];
//				 w2 = omega[j + 1];
//				 c1 = c[j];
//				 c2 = c[j + 1];
//				 break;
//			}
//		}
//
//		
//		printf("R=%lf", R);
//		cout << endl;
//		printf("w1=%lf",w1);
//		cout << endl;
//		printf("w2=%lf", w2);
//		cout << endl;
//	
//		cres[k] = LineInter(w1, w2, c1, c2, R);
//		cout << cres[k]<< endl;
//	
//	
//	} 
//	return 0;
//}
//

//12345679813245678123456789123456789
//
//
//
///*
//
//прямое преобразование
//
//*/
//
//#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <time.h>
//
//#define M_PI 3.14159265358979323846
//
//using namespace std;
//
//static int N = 200; // количество отсчётов
//double LineInter(double w1, double w2, double c1, double c2, double R) {
//
//	double res;
//
//	res = c1 + (c2 - c1) / (w2 - w1)*(R - w1);
//
//
//
//	return res;
//}
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
//	double omega[15] = { 0, 77022.06, 115510.3815, 153971.203, 383590.5343, 755563.7828, 1089856.852, 1334024.614, 1493426.204, 1636126.952, 1789995.119, 1958672.888, 2331441.97, 2734717.393, 3155165.375 };
//	double c[15] = { 4904.144463, 4903.376713, 4902.412853, 4901.055611, 4884.026372, 4810.068434, 4625.496139, 4246.332231, 3802.978599, 3471.969227, 3255.846529, 3117.32472, 2968.484112, 2901.625273,2869.486661 };
//	double w1, w2, c1, c2;
//	double Pn;
//	double n;
//	double cres[100];
//	double C[101];
//	double phi[101];
//	srand(time(NULL));
//
//	stepOfFreq = 1.0 / (N * 0.05);
//
//	//TO DO:
//	fstream file;
//
//	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfInputSignal.txt");
//	for (int i = 0; i < N; i++) {
//		file >> tmp >> arrayOfFuncResult[i];
//	}
//	file.close();
//
//	t = 0.0;
//	/*
//	Прямое преобразование
//	*/
//
//	for (int k = 0; k <= N / 2; k++) {
//		for (int n = 0; n < N; n++) {
//
//			argumentOfTrigFunc = 2 * M_PI * k * n / N;
//			Re = arrayOfFuncResult[n] * cos(argumentOfTrigFunc);
//			Im = -arrayOfFuncResult[n] * sin(argumentOfTrigFunc);
//			sumOfRealPart += Re;
//			sumOfImaginaryPart += Im;
//
//		}
//		arrayOfComplexNumbers[k][0] = sumOfRealPart;
//		arrayOfComplexNumbers[k][1] = sumOfImaginaryPart;
//		sumOfRealPart = 0.0;
//		sumOfImaginaryPart = 0.0;
//
//	}
//
//	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfOutputSignalRe.txt");//RE
//	for (int i = 0; i <= N / 2; i++) {
//		file << t << "\t" << arrayOfComplexNumbers[i][0] << "\n";
//		t += stepOfFreq;
//	}
//	file.close();
//
//	t = 0.0;
//
//	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfOutputSignalIm.txt");//Im
//	for (int i = 0; i <= N / 2; i++) {
//		file << t << "\t" << arrayOfComplexNumbers[i][1] << "\n";
//		t += stepOfFreq;
//	}
//	file.close();
//
//	file.open("C:\\Users\\Олег\\Desktop\\2019\\PointsOfComplexCoord.txt"); //mag phase
//	for (int k = 0; k <= N / 2; k++) {
//		mag[k] = sqrt(arrayOfComplexNumbers[k][0] * arrayOfComplexNumbers[k][0] + arrayOfComplexNumbers[k][1] * arrayOfComplexNumbers[k][1]);
//		phase[k] = atan2(arrayOfComplexNumbers[k][1], arrayOfComplexNumbers[k][0]);
//
//		double R = (((double)rand()) / RAND_MAX) * 3155165;
//
//		for (int j = 0; j <= 13; j++) {
//			if ((R > omega[j]) && (R < omega[j + 1])) {
//				w1 = omega[j];
//				w2 = omega[j + 1];
//				c1 = c[j];
//				c2 = c[j + 1];
//				break;
//			}
//		}
//
//
//		/*	printf("R=%lf", R);
//		cout << endl;
//		printf("w1=%lf", w1);
//		cout << endl;
//		printf("w2=%lf", w2);
//		cout << endl;*/
//
//		cres[k] = LineInter(w1, w2, c1, c2, R);
//		/*cout << cres[k] << endl;*/
//		file << mag[k] << "\t" << phase[k] << "\t" << R << "\t" << cres[k] << "\n";
//	}
//	file.close();
//
//	file.open("C:\\Users\\Олег\\Desktop\\2019\\PnNew.txt"); //C[k] phi
//	for (int j = 0; j <= N - 1; j++) {
//		/*	C[j] = sqrt(arrayOfComplexNumbers[j][0] * arrayOfComplexNumbers[j][0] + arrayOfComplexNumbers[j][1] * arrayOfComplexNumbers[j][1]);
//		phi[j] = atan2(arrayOfComplexNumbers[j][1], arrayOfComplexNumbers[j][0]);*/
//
//		Pn = mag[k] * cos(2 * M_PI)*j*(n + phase[k]) / N;
//
//
//
//		file << << "\t" << << "\t" << "\n";
//	}
//	file.close();
//
//
//	for (int k = 0; k <= N / 2; k++) {
//		for (int n = 0; n < N; n++) {
//
//			argumentOfTrigFunc = 2 * M_PI * k * n / N;
//			Re = arrayOfFuncResult[n] * cos(argumentOfTrigFunc);
//			Im = -arrayOfFuncResult[n] * sin(argumentOfTrigFunc);
//			sumOfRealPart += Re;
//			sumOfImaginaryPart += Im;
//
//		}
//		arrayOfComplexNumbers[k][0] = sumOfRealPart; //Ak
//		arrayOfComplexNumbers[k][1] = sumOfImaginaryPart; //Bk
//		sumOfRealPart = 0.0;
//		sumOfImaginaryPart = 0.0;
//
//	}
//
//
//
//
//
//	file.open("C:\\Users\\Олег\\Desktop\\2019\\IsReAndImEqualAfterFunc.txt");
//	for (int k = 0; k <= N / 2; k++) {
//		if (mag[k] * cos(phase[k]) != arrayOfComplexNumbers[k][0])
//			Re = -mag[k] * cos(phase[k]);
//		else
//			Re = mag[k] * cos(phase[k]);
//		if (mag[k] * sin(phase[k]) != arrayOfComplexNumbers[k][1])
//			Im = -mag[k] * sin(phase[k]);
//		else
//			Im = mag[k] * sin(phase[k]);
//		file << Re << "\t" << arrayOfComplexNumbers[k][0] << "\t\t\t" << Im << "\t" << arrayOfComplexNumbers[k][1] << "\n";
//	}
//	file.close();
//
//	return 0;
//}