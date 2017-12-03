#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<math.h>
using namespace std;
double temperature[200][200][200];
double intermediateTemp[200][200];
double temp, hX, hY, hT, lX, lY, lT;
int maxX, maxY, maxT;
double Func(int numX, int numY, double numT) {
	return 2 * sin(numX * hX + numY * hY) + 2 * numT*hT / ((numT*hT)*(numT*hT) + 1);
}
double RealValue(int numX, int numY, double numT) {
	return sin(numX * hX + numY * hY) + log((numT*hT)*(numT*hT) + 1);
}
double L1(int numX, int numY, int numT){
	return (temperature[numT][numY][numX - 1] - 2 * temperature[numT][numY][numX] + temperature[numT][numY][numX + 1]) /( hX* hX);
}
double L1(double buttomValue, double middleValue, double topValue) {
	return (buttomValue - 2 * middleValue + topValue) / (hX* hX);
}

double L2(int numX, int numY, int numT) {
	return (temperature[numT][numY - 1][numX] - 2 * temperature[numT][numY][numX] + temperature[numT][numY + 1][numX]) / (hY* hY);
}
double mu1(int numY, int numT) {
	return sin(numY*hY) + log((numT*hT)*(numT*hT) + 1);
}
double mu2(int numY, int numT) {
	return sin(numY*hY) + log((numT*hT)*(numT*hT) + 1);
}
double _mu1(int numY, int numT) {
	double L2mu1Top = (mu1(numY - 1, numT + 1) - 2 * mu1(numY, numT + 1) + mu1(numY + 1, numT + 1)) / (hY*hY);
	double L2mu1 = (mu1(numY - 1, numT) - 2 * mu1(numY, numT) + mu1(numY + 1, numT)) / (hY*hY);
	return (mu1(numY, numT) + mu1(numY, numT + 1)) / 2 - (hT / 4)*(L2mu1Top - L2mu1);
}
double _mu2(int numY, int numT) {
	double L2mu2Top = (mu2(numY - 1, numT + 1) - 2 * mu2(numY, numT + 1) + mu2(numY + 1, numT + 1)) / (hY*hY);
	double L2mu2 = (mu2(numY - 1, numT) - 2 * mu2(numY, numT) + mu2(numY + 1, numT)) / (hY*hY);
	return (mu2(numY, numT) + mu2(numY, numT + 1)) / 2 - (hT / 4)* (L2mu2Top - L2mu2);
}
double mu3(int numX, int numT) {
	return sin(numX*hX) + log((numT*hT)*(numT*hT) + 1);
}
double mu4(int numX, int numT) {
	return sin(numX*hX) + log((numT*hT)*(numT*hT) + 1);
}
double startCondition(double x, double y) {
	return sin(x + y);
}
double alfa[200], beta[200], A, B, C, F;
int main() {
	freopen("out.txt", "w", stdout);
	lX = 2 * M_PI;
	lY = 2 * M_PI;
	lT = 10;
	hX = hY = M_PI / 5;
	hT = 0.25;
	maxX = int(lX / hX);
	maxY = int(lY / hY);
	maxT = int(lT / hT);
	double sigmaX = hT / (hX*hX);
	double sigmaY = hT / (hY*hY);
	for (int i = 0; i <= maxY; i++)
		for (int j = 0; j <= maxX; j++) {
		temperature[0][i][j] = startCondition(j*hX, i*hY);
		}
	for (int temp = 0; temp < maxT; temp++) {
		for (int m = 1; m < maxY; m++) {
			alfa[1] = 0; intermediateTemp[m][0] = beta[1] = _mu1(m, temp);
			intermediateTemp[m][maxX] = _mu2(m, temp);
			A = -sigmaX;
			C = 2 * (1 + sigmaX);
			B = -sigmaX;
			for (int n = 1; n < maxX; n++) {
				F = 2 * temperature[temp][m][n] + hT * L2(n, m, temp) + hT *Func(n, m, temp + 0.5);
				//F = sigmaY*(temperature[temp][m - 1][n] + temperature[temp][m + 1][n]) + 2 * (1 - sigmaY)*temperature[temp][m][n] + hT * Func(n, m, temp + 0.5);
				alfa[n + 1] = -B / (A*alfa[n] + C);
				beta[n + 1] = (F - A * beta[n]) / (A*alfa[n] + C);
			}
			for (int n = maxX - 1; n > 0; n--) {
				intermediateTemp[m][n] = alfa[n + 1] * intermediateTemp[m][n + 1] + beta[n + 1];
			}
		}
		//// считаем граничные условия по X 
		for (int m = 0; m <= maxY; m++) {
			temperature[temp + 1][m][0] = mu1(m, temp + 1);
			temperature[temp + 1][m][maxX] = mu2(m, temp + 1);
		}
		for (int n = 1; n < maxX; n++) {
			////переход на temp+1 уровень
			alfa[1] = 0; temperature[temp + 1][0][n] = beta[1] = mu3(n, temp + 1);
			temperature[temp + 1][maxY][n] = mu4(n, temp + 1);
			A = -sigmaY;
			B = -sigmaY;
			C = 2 * (1 + sigmaY);
			for (int m = 1; m < maxY; m++) {
				F = 2 * intermediateTemp[m][n] + hT * L1(intermediateTemp[m][n - 1], intermediateTemp[m][n], intermediateTemp[m][n + 1]) + hT *Func(n, m, temp + 0.5);
				//F = sigmaX*(intermediateTemp[m][n - 1] + intermediateTemp[m][n + 1]) + 2 * (1 - sigmaX)*intermediateTemp[m][n] + hT * Func(n, m, temp + 0.5);
				alfa[m + 1] = -B / (A*alfa[m] + C);
				beta[m + 1] = (F - A * beta[m]) / (A*alfa[m] + C);
			}
			for (int m = maxY - 1; m > 0; m--) {
				temperature[temp + 1][m][n] = alfa[m + 1] * temperature[temp + 1][m + 1][n] + beta[m + 1];
			}
		}
	}
	double err = 0;
	for (int t = 0; t <= maxT; t++) {
		for (int m = 0; m <= maxY; m++) {
			for (int n = 0; n <= maxX; n++) {
				//cout << temperature[maxT][m][n] << " " << RealValue(n, m, maxT) << endl;
				if (err < abs(temperature[t][m][n] - RealValue(n, m, t))) {
					err = abs(temperature[t][m][n] - RealValue(n, m, t));
				}

			}
		//	cout << endl;
		}
	}
	cout << err;
}