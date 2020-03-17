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

const int _A = 8, _B = 4, _C = 3, f = _B;
const float A = 1.0, fi = _C * M_PI, t0 = 0, tN = _A;

float s(float t, float A, float f, float fi) {
	return A * sin(2*M_PI*f*t + fi);
}

int quant(float st, unsigned int q) {
	int high = pow(2, q);
	int low = 0;

	return ((st + A) / (2 * A)) * (high - low) + low;
}

void zad1() {
	float tn = t0;
	float fs = 1000 * f;
	float Ts = 1.0 / fs;

	int n = 0;
	string title = "", filename = "";
	stringstream st;

	while (tn <= tN) {
		st << tn << ";" << s(tn, A, f, fi) << ";\n";
		n++;
		tn = t0 + (n*Ts);
	}

	// s(t)
	string outS = st.str();

	title = "s(t), t <" + to_string(t0) + ";" + to_string(tN) + ">, Ts = " + to_string(Ts);
	filename = "LAB_02_ZAD1_s";

	dataToCsv("csv/" + filename + ".csv", outS);
	drawChart(title, "t[s]", "A[v]", "csv/" + filename +".csv", "charts/" + filename + ".png");
}

void zad2() {
	float tn = t0;
	float fs = 1000 * f;
	float Ts = 1.0 / fs;

	int n = 0;
	string title = "", filename = "";
	stringstream st;

	while (tn <= tN) {
		st << tn << ";" << quant(s(tn, A, f, fi), 16) << ";\n";
		n++;
		tn = t0 + (n*Ts);
	}

	// s(t) skwantyzowane, q = 16
	string outS = st.str();

	title = "s(t) skwantyzowane, q = 16";
	filename = "LAB_02_ZAD2_s_q_16";

	dataToCsv("csv/" + filename + ".csv", outS);
	drawChartWithDots(title, "t[s]", "przedzialy kwantyzacji", "csv/" + filename + ".csv", "charts/" + filename + ".png");
}

void test() {
	int f = 1;
	float tn = t0;
	float fs = 10;
	float Ts = 1.0 / fs;

	int n = 0;
	string title = "", filename = "", filename2 = "";
	stringstream st, points;

	int qPrev = quant(s(tn, A, f, fi), 3);
	points << tn << ";" << qPrev << ";\n";

	while (tn <= tN) {
		if (quant(s(tn, A, f, fi), 3) != qPrev) {
			qPrev = quant(s(tn, A, f, fi), 3);
			points << tn << ";" << qPrev << ";\n";
		}

		st << tn << ";" << quant(s(tn, A, f, fi), 3) << ";\n";
		n++;
		tn = t0 + (n*Ts);
	}

	// s(t) skwantyzowane, q = 16
	string outS = st.str();
	string outSPoints = points.str();

	title = "s(t) skwantyzowane, q = 16 TEST";
	filename = "LAB_02_ZAD2_s_q_16_test";
	filename2 = "LAB_02_ZAD2_s_q_16_points_test";

	dataToCsv("csv/" + filename + ".csv", outS);
	dataToCsv("csv/" + filename2 + ".csv", outSPoints);
	drawChartWithStepsPoints(title, "t", "s(t)", "csv/" + filename + ".csv", "csv/" + filename2 + ".csv", "charts/" + filename + ".png");
}

void zad3() {
	float tn = t0;
	float fs = 500 * f;
	float Ts = 1.0 / fs;

	int n = 0;
	string title = "", filename = "";
	stringstream st;

	while (tn <= tN) {
		st << tn << ";" << quant(s(tn, A, f, fi), 8) << ";\n";
		n++;
		tn = t0 + (n*Ts);
	}

	// s(t) skwantyzowane, q = 8
	string outS = st.str();

	title = "s(t) skwantyzowane, q = 8";
	filename = "LAB_02_ZAD3_s_q_8";

	dataToCsv("csv/" + filename + ".csv", outS);
	drawChartWithDots(title, "t[s]", "przedzialy kwantyzacji", "csv/" + filename + ".csv", "charts/" + filename + ".png");
}

int main() {
	zad1();
	zad2();
	zad3();
	//test();
	getchar();
	return 0;
}