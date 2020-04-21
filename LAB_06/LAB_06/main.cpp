#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <complex>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include "chart_utils.h"
#include "gen_charts.h"
#include "gen_signals.h"
#include "dft.h"

using namespace chart_utils;
using std::cout;
using std::endl;
using std::string;
using std::to_string;
using std::stringstream;
using std::vector;
using std::complex;
using std::bitset;

// stale uzywane w calosci programu
float Tb = 0.1, fi = 0, a = 1.0;

float s_n(float a, float f, float tn, float fi) {
	return (a * sin(2 * M_PI*f*tn + fi));
}

float calcIntegral(float x, float a, float b, int n) {
	float step = (b - a) / n;
	float area = 0.0;

	for (int i = 0; i < n; i++) {
		area += x * step;
	}

	return area;
}

/*float calcIntegral(float x, float Tb, float fs) {
	float dt = 1 / fs;
	int num_samples = int(Tb / dt);
	float area = 0.0;

	float t0 = 0.0;
	int n = 0;
	float tn = t0 + (n*dt);

	while (tn <= Tb) {
		area += x * dt;
		n++;
		tn = t0 + (n*dt);
	}

	return area;
}*/

void demodulateAmplPhaseModSignal(vector<float> &m_t, vector<float> &m_x, vector<float> z_t, vector<float> z_x, float Tb, float h, float a, float f, float fs, float fi, float(*s_n)(float, float, float, float), string namePrefix) {
	float t0 = 0.0;
	int n = 0;
	float dt = 1 / fs;
	float tn = t0 + (n*dt);
	const int num_samples = int(Tb / dt);

	// etap 1: x(t)
	vector<float> t_x, x_x;

	for (int i = 0; i < z_x.size(); i++) {
		x_x.push_back(z_x[i] * s_n(a, f, tn, fi));
		t_x.push_back(tn);
		n++;
		tn = t0 + (n*dt);
		//cout << t_x[i] << "; " << x_x[i] << endl;
	}

	//cout << namePrefix << "xt" << endl;

	drawSignalChart(t_x, x_x, 10 * num_samples, "lines", "x(t)", namePrefix + "xt", "T_b[s]", "A", 1920, 800);

	// etap 2: p(t)
	// calka
	vector<float> t_p, x_p;

	n = 0;
	tn = t0 + (n*dt);
	
	
	float sum = 0.0f;
	float t = 0.0f;
	for (int i = 0; i < x_x.size(); i++) {
		sum += x_x[i] / dt;
		x_p.push_back(sum);
		t_p.push_back(tn);
		n++;
		tn = t0 + (n*dt);
		t += dt;
		if (t > Tb) {
			t = 0.0f;
			sum = 0.0f;
		}

		cout << t_p[i] << "; " << x_p[i] << endl;
	}

	drawSignalChart(t_p, x_p, 10 * num_samples, "lines", "p(t)", namePrefix + "xp", "T_b[s]", "A", 1920, 800);

	// etap 3: m'(t)
	n = 0;
	tn = t0 + (n*dt);

	for (int i = 0; i < x_p.size(); i++) {
		if (x_p[i] >= h)
			m_x.push_back(1);
		else
			m_x.push_back(0);

		m_t.push_back(tn);
		n++;
		tn = t0 + (n*dt);
	}
}

void demodulateFreqModSignal(vector<float> &m_t, vector<float> &m_x, vector<float> z_t, vector<float> z_x, float Tb, float h, float a, float f1, float f2, float fs, float fi, float(*s_n)(float, float, float, float), string namePrefix) {
	float t0 = 0.0;
	int n = 0;
	float dt = 1 / fs;
	float tn = t0 + (n*dt);
	const int num_samples = int(Tb / dt);

	// etap 1: x1(t) i x2(t)
	vector<float> z_t_copy, z_x_copy;

	z_t_copy = z_t;
	z_x_copy = z_x;


	vector<float> t_x, t_x_copy, x_x, x_x_copy;

	for (int i = 0; i < z_x.size(); i++) {
		x_x.push_back(z_x[i] * s_n(a, f1, tn, fi));
		t_x.push_back(tn);

		x_x_copy.push_back(z_x[i] * s_n(a, f2, tn, fi));
		t_x_copy.push_back(tn);

		n++;
		tn = t0 + (n*dt);
		//cout << t_x[i] << "; " << x_x[i] << endl;
	}

	drawSignalChart(t_x, x_x, 10 * num_samples, "lines", "x_1(t)", namePrefix + "x1_t", "T_b[s]", "A", 1920, 800);
	drawSignalChart(t_x_copy, x_x_copy, 10 * num_samples, "lines", "x_2(t)", namePrefix + "x2_t", "T_b[s]", "A", 1920, 800);

	// etap 2: p(t)
	// calka
	vector<float> t_p, t_p_copy, x_p, x_p_copy;

	n = 0;
	tn = t0 + (n*dt);

	for (int i = 0; i < x_x.size(); i++) {
		x_p.push_back(calcIntegral(x_x[i], 0, Tb, num_samples));
		t_p.push_back(tn);

		x_p_copy.push_back(calcIntegral(x_x_copy[i], 0, Tb, num_samples));
		t_p_copy.push_back(tn);

		n++;
		tn = t0 + (n*dt);
	}

	for (int i = 0; i < x_p.size(); i++) {
		x_p[i] = x_p[i] - x_p_copy[i];
	}

	drawSignalChart(t_p, x_p, 10 * num_samples, "lines", "p(t)", namePrefix + "xp", "T_b[s]", "A", 1920, 800);

	// etap 3: m'(t)
	n = 0;
	tn = t0 + (n*dt);

	for (int i = 0; i < x_p.size(); i++) {
		if (x_p[i] >= h)
			m_x.push_back(1);
		else
			m_x.push_back(0);

		m_t.push_back(tn);
		n++;
		tn = t0 + (n*dt);
	}
}

void reverseBitset(bitset<8> &bits) {
	for (int i = 0; i < 4; ++i) {
		bool tmp = bits[i];
		bits[i] = bits[7 - i];
		bits[7 - i] = tmp;
	}
}

vector<bitset<8>> strToBinStream(string input, string option) {
	vector<bitset<8>> output;

	for (int i = 0; i < input.size(); i++) {
		bitset<8> tmp(input[i]);

		if (option == "bigEndian")
			reverseBitset(tmp);

		output.push_back(tmp);
	}

	for (int i = 0; i < output.size(); i++)
		cout << output[i];

	cout << endl;

	return output;
}


void zad1_2(vector<float> t, vector<float> x) {
	int _N = 2;
	float f = _N * pow(Tb, -1), fs = 20 * f;

	float a1 = 0.0, a2 = 1.0; // ASK
	float f1 = (_N + 1) / Tb, f2 = (_N + 2) / Tb; // FSK
	float fi0 = 0, fi1 = M_PI; // PSK

	float dt = 1 / fs;
	const int num_samples = int(Tb / dt);

	string namePrefix = "LAB_06_ZAD_1_";

	// sygnal z_A(t), kluczowanie amplitudy
	// ----------------------------------------
	vector<float> x_A, t_A;
	generateAmplModSignal(t_A, x_A, x, a1, a2, Tb, f, fs, fi);

	// wykres sygnalu zmodulowanego z_A(t)
	stringstream title_A;

	title_A << std::fixed << std::setprecision(2)
		<< "z_A(t): sygnal zmodulowany amplituda (ASK)"
		<< ", A_1 = " << a1
		<< ", A_2 = " << a2
		<< ", T_b = " << Tb << " s"
		<< ", f = " << f << " Hz"
		<< ", N = " << _N;

	string filename_A = namePrefix + "sygnal_z_At";
	drawSignalChart(t_A, x_A, 10 * num_samples, "lines", title_A.str(), filename_A, "T_b[s]", "A", 1920, 800);

	// sygnal z_F(t), kluczowanie czestotliwosci
	// ---------------------------------------------
	vector<float> x_F, t_F;
	generateFreqModSignal(t_F, x_F, x, a, f1, f2, Tb, fs, fi);

	// wykres sygnalu zmodulowanego z_F(t)
	stringstream title_B;

	title_B << std::fixed << std::setprecision(2)
		<< "z_F(t): sygnal zmodulowany czestotliwoscia (FSK)"
		<< ", A = " << a
		<< ", f_1 = " << f1 << " Hz"
		<< ", f_2 = " << f2 << " Hz"
		<< ", N = " << _N;

	string filename_B = namePrefix + "sygnal_z_Ft";
	drawSignalChart(t_F, x_F, 10 * num_samples, "lines", title_B.str(), filename_B, "T_b[s]", "A", 1920, 800);

	// sygnal z_P(t), kluczowanie fazy
	// -----------------------------------
	vector<float> x_P, t_P;
	generatePhaseModSignal(t_P, x_P, x, a, fi0, fi1, Tb, f, fs);

	// wykres sygnalu zmodulowanego z_P(t)
	stringstream title_C;

	title_C << std::fixed << std::setprecision(2)
		<< "z_P(t): sygnal zmodulowany faza (PSK)"
		<< ", A = " << a
		<< ", fi_0 = " << fi0 << " Hz"
		<< ", fi_1 = " << fi1 << " Hz"
		<< ", f = " << f << " Hz"
		<< ", N = " << _N;

	string filename_C = namePrefix + "sygnal_z_Pt";
	drawSignalChart(t_P, x_P, 10 * num_samples, "lines", title_C.str(), filename_C, "T_b[s]", "A", 1920, 800);

	// wykres zbiorczy sygnalow
	// ----------------------------
	stringstream title_inf;

	title_inf << std::fixed << std::setprecision(2)
		<< "sygnal informacyjny"
		<< ", T_b = " << Tb << " s"
		<< ", f = " << f << " Hz"
		<< ", fs = " << fs << " Hz";

	drawSignalMultiChart(t, t_A, t_F, t_P, x, x_A, x_F, x_P, 10 * num_samples, "lines", "Sygnaly: inf., z_A(t), z_F(t), z_P(t), N = " + to_string(_N), title_inf.str(), title_A.str(), title_B.str(), title_C.str(), namePrefix + "sygnaly", "LAB_06_sygnal_inf", namePrefix + "sygnal_z_At", namePrefix + "sygnal_z_Ft", namePrefix + "sygnal_z_Pt", "T_b[s]", "A", 1920, 1200);

	// sygnaly zdemodulowane i ich wykresy
	// ---------------------------------------
		float h = 0.0001;
		vector<float> m_t, m_x;

		demodulateAmplPhaseModSignal(m_t, m_x, t_A, x_A, Tb, h, a2, f, fs, fi, s_n, namePrefix + "sygnal_z_At_");

		//for (int i = 0; i < m_x.size(); i++)
		//	cout << m_t[i] << "; " << m_x[i] << endl;

		float ymin = *min_element(m_x.begin(), m_x.end());
		float ymax = *max_element(m_x.begin(), m_x.end());

		drawSignalChart(m_t, m_x, 10 * num_samples, "lines", "m'(t)", namePrefix + "sygnal_mprim_t", "T_b[s]", "A", 1920, 800, 0, ymax * 0.1);

		vector<float> m2_t, m2_x;

		demodulateFreqModSignal(m2_t, m2_x, t_F, x_F, Tb, h, a, f1, f2, fs, fi, s_n, namePrefix + "sygnal_zFt_");

		//for (int i = 0; i < m2_x.size(); i++)
		//	cout << m2_t[i] << "; " << m2_x[i] << endl;

		ymin = *min_element(m2_x.begin(), m2_x.end());
		ymax = *max_element(m2_x.begin(), m2_x.end());

		drawSignalChart(m2_t, m2_x, 10 * num_samples, "lines", "m'(t)", namePrefix + "sygnal_mprim2_t", "T_b[s]", "A", 1920, 800, 0, ymax * 0.1);

		vector<float> m3_t, m3_x;

		demodulateAmplPhaseModSignal(m3_t, m3_x, t_P, x_P, Tb, h, a, f, fs, fi1, s_n, namePrefix + "sygnal_zPt_");

		ymin = *min_element(m3_x.begin(), m3_x.end());
		ymax = *max_element(m3_x.begin(), m3_x.end());

		drawSignalChart(m3_t, m3_x, 10 * num_samples, "lines", "m'(t)", namePrefix + "sygnal_mprim3_t", "T_b[s]", "A", 1920, 800, 0, ymax * 0.1);
}

int main() {
	int _N = 2;
	float f = _N * pow(Tb, -1), fs = 20 * f;

	float dt = 1 / fs;
	const int num_samples = int(Tb / dt);

	vector<bitset<8>> s1 = strToBinStream("KOT", "littleEndian");

	// sygnal informacyjny
	vector<float> t, x;
	generateInfSignal(t, x, s1, Tb, fs);

	float ymin = *min_element(x.begin(), x.end());
	float ymax = *max_element(x.begin(), x.end());

	// wykres sygnalu inf.
	stringstream title;

	title << std::fixed << std::setprecision(2)
		<< "sygnal informacyjny"
		<< ", T_b = " << Tb << " s"
		<< ", f = " << f << " Hz"
		<< ", fs = " << fs << " Hz";

	string filename = "LAB_06_sygnal_inf";
	drawSignalChart(t, x, 10 * num_samples, "lines", title.str(), filename, "T_b[s]", "A", 1920, 800, 0, ymax * 0.1);

	// zad. 1. i 2.
	// -----------
		zad1_2(t, x);


	getchar();
	return 0;
}
