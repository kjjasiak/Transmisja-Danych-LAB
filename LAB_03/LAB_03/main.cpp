#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <complex>
#include "chart_utils.h"

using namespace chart_utils;
using std::cout;
using std::endl;
using std::string;
using std::to_string;
using std::stringstream;
using std::vector;
using std::complex;

const int _A = 8, _B = 4, _C = 3, f = _B;
const float A = 1.0, fi = _C * M_PI, n0 = 0, nN = 843;


float s(float t, float A, float f, float fi) {
	return A * sin(2 * M_PI*f*t + fi);
}

complex<float> x(int N) {
	int k = N - 1;
	complex<float> wNPow(0, (2 * M_PI) / N);
	complex<float> wN = exp(wNPow);

	complex<float> xk(0, 0);

	for (int n = 0; n <= N - 1; n++)
		xk += s(n, A, f, fi) * pow(wN, -k * n);

	//cout << real(xk) << " + i " << imag(xk) << endl;

	return xk;
}

vector<complex<float>> dft(int N) {
	vector<complex<float>> xk;

	//cout << real(wN) << " + i " << imag(wN) << endl;

	for (int k = 0; k <= N - 1; k++) {
		xk.push_back(x(k));
	}

	//complex<float> a(3.0, 4.0);
	//complex<float> b(5.2, 3.1);

	//xk.push_back(a);
	//xk.push_back(b);

	return xk;
}

void zad2() {

}

int main() {
	vector<complex<float>> xk = dft(nN);

	for (int i = 0; i < xk.size(); i++) {
		cout << real(xk[i]) << " + i " << imag(xk[i]) << endl;
	}

	getchar();
	return 0;
}