#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include "chart_utils.h"

using namespace chart_utils;
using std::string;
using std::to_string;
using std::stringstream;

int _A = 8, _B = 4, _C = 3, f = _B;
float A = 1.0, fi = _C * M_PI;

float s(float t, float A, float f, float fi) {
	return A * sin(2*M_PI*f*t + fi);
}

void zad1() {
	float t0 = 0, tN = _A, tn = t0, fs = 20, Ts = 1.0 / fs, dt = Ts;
	int n = 0;
	string title = "";
	stringstream st;

	while (tn <= tN) {
		st << tn << ";" << s(tn, A, fs, fi) << ";\n";
		n++;
		tn = t0 + (n*dt);
	}

	// s(t)
	string outS = st.str();

	title = "s(t), t <" + to_string(t0) + ";" + to_string(tN) + ">, Ts = " + to_string(Ts);

	dataToCsv("csv/LAB_02_ZAD1_s.csv", outS);
	drawChart(title, "t", "s(t)", "csv/LAB_02_ZAD1_s.csv", "charts/LAB_02_ZAD1_s.png");

	//return st.str();
}

int main() {
	zad1();
	getchar();
	return 0;
}