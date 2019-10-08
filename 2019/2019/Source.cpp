/*

строит функцию

*/

#include <iostream>
#include <fstream>
#include <cmath>

#define M_PI 3.14159265358979323846

using namespace std;

const  int N = 200; // количество отсчётов
 const double eps = 1e-24;

 //double myFunc(double t) {  //прямоуг
	// t = t*1e5;

	// t = t - 2;
	//							
	//
	// if (t > -eps && t < 2)
	//
	//	 return 1.0;
	//	
	// else {
	//
	//	 return 0.0;
	//							}
	//						}


 //double myFunc(double t) {  //sin^2
 //	/*return sin(2.0 * M_PI * t);*/
	// t = t*1e5;
	// t = t - 2;
	// 	
	//	
	// 	if (t >= 0 && t <= 5 )
	//		return  (sin(2 * M_PI * t / 10 ) * (sin(2 * M_PI * t / 10)));
	// 	else {
	// 		return 0.0;
	// 	}
	//
	//	
 //} 



//double myFunc(double t) {  //sinc
//	/*return sin(2.0 * M_PI * t);*/
//	t = t - 0.0005;
//	t = t*1e4;
//	
//	if (fabs(t) < eps) {
//		return 1;
//	}
//	else {
//		
//		return  1 * ((sin(0.5 * M_PI * t)) / (0.5 * M_PI * t));
//	}
//} 

//double myFunc(double t) {
//t = t*1e4;
//	if (t >= 0)
//		return exp(-t) * 1.0 + 1.5 ;
//	else return exp(-t) * 0.0 + 1.5 ;
//}

//double myFunc(double t) {  //прямоуг
//	
//	/*t = t - 0.0001;*/
//	t = t*1e5;
//	t = t - 2;
//	if (t > -eps && t < 1.5)
//		return 1.0;
//	else {
//		return 0.0;
//	}
//}

double myFunc(double t) { //треуг
	t = t*1e5;
	t = t - 2;
	
	if (abs(t) <= 1)
		return 1 -  1 * abs(t);
	else return 0.0;
}

//double myFunc(double t) {
// 
//	if (t < 1) {
//		return 1.0;
//	}
//	else if (t = 1) {
//		return 1.0 / 2.0;
//	}
//	else if (t > 1) {
//		return 0;
//	}
//}


int main() {
	/*const  int N = 200;*/
	double sumOfRealPart = 0.0, sumOfImaginaryPart = 0.0;
	double arrayOfFuncResult[200];
	double t = 0.0;
    double const dT = 0.000005 /* * 10.0*/;

	


	fstream file;

	file.open("C:\\Users\\Олег\\Desktop\\2019\\DataOfInputSignal.txt");
	for (int i = 0; i < N; i++) {
		arrayOfFuncResult[i] = myFunc(t) /*2 * cos(2 * M_PI * t)*/;
		file << t << "\t" << arrayOfFuncResult[i] << "\n";
		t += dT;
	}
	file.close();

	return 0;
}