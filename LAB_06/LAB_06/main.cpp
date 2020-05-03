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

void demodulateAmplPhaseModSignal(vector<float> &mprim_xaxis, vector<float> &mprim_yaxis, vector<float> z_xaxis, vector<float> z_yaxis, float Tb, float h, float a, float f, float fs, float fi, float(*s_n)(float, float, float, float), string namePrefix, string signalModType) {
	float t0 = 0.0;
	int n = 0;
	float dt = 1 / fs;
	float tn = t0 + (n*dt);
	const int num_samples = int(Tb / dt);

	// etap 1: x(t)
	vector<float> x_xaxis, x_yaxis;

	for (int i = 0; i < z_yaxis.size(); i++) {
		x_yaxis.push_back(z_yaxis[i] * s_n(a, f, tn, fi));
		x_xaxis.push_back(tn);
		n++;
		tn = t0 + (n*dt);
	}

	drawSignalChart(x_xaxis, x_yaxis, 10 * num_samples, "lines", signalModType + ": demodulacja, etap 1, x(t)", namePrefix + "demod_etap_1_xt", "T_b[s]", "A", 1920, 800);

	// etap 2: p(t), calka
	vector<float> p_xaxis, p_yaxis;

	n = 0;
	tn = t0 + (n*dt);
	
	float sum = 0.0, t = 0.0;

	for (int i = 0; i < x_yaxis.size(); i++) {
		sum += x_yaxis[i] / dt;
		p_yaxis.push_back(sum);
		p_xaxis.push_back(tn);
		n++;
		tn = t0 + (n*dt);
		t += dt;

		if (t > Tb) {
			t = 0.0;
			sum = 0.0;
		}
	}

	drawSignalChart(p_xaxis, p_yaxis, 10 * num_samples, "lines", signalModType + ": demodulacja, etap 2, p(t)", namePrefix + "demod_etap_2_pt", "T_b[s]", "x(t)dt", 1920, 800);

	// etap 3: m'(t)
	n = 0;
	tn = t0 + (n*dt);

	for (int i = 0; i < p_yaxis.size(); i++) {
		if (p_yaxis[i] >= h)
			mprim_yaxis.push_back(1);
		else
			mprim_yaxis.push_back(0);

		mprim_xaxis.push_back(tn);
		n++;
		tn = t0 + (n*dt);
	}

	float ymin = *min_element(mprim_yaxis.begin(), mprim_yaxis.end());
	float ymax = *max_element(mprim_yaxis.begin(), mprim_yaxis.end());

	drawSignalChart(mprim_xaxis, mprim_yaxis, 10 * num_samples, "lines", signalModType + ": demodulacja, etap 3, m'(t), h = " + to_string(int(h)), namePrefix + "demod_etap_3_mprimt", "T_b[s]", "A", 1920, 800, 0, ymax * 0.1);

	drawSignalMultiChart(
		z_xaxis, x_xaxis, p_xaxis, mprim_xaxis,
		z_yaxis, x_yaxis, p_yaxis, mprim_yaxis,
		10 * num_samples,
		"lines",
		signalModType + ": demodulacja",
		signalModType,
		signalModType + ": demodulacja, etap 1, x(t)",
		signalModType + ": demodulacja, etap 2, p(t)",
		signalModType + ": demodulacja, etap 3, m'(t)",
		namePrefix + "demod_sygnaly",
		namePrefix.substr(0, namePrefix.length() - 1),
		namePrefix + "demod_etap_1_xt",
		namePrefix + "demod_etap_2_pt",
		namePrefix + "demod_etap_3_mprimt",
		"T_b[s]", "A", "T_b[s]", "x(t)dt",
		1920, 1600
	);
}

void demodulateFreqModSignal(vector<float> &mprim_xaxis, vector<float> &mprim_yaxis, vector<float> z_xaxis, vector<float> z_yaxis, float Tb, float h, float a, float f1, float f2, float fs, float fi, float(*s_n)(float, float, float, float), string namePrefix, string signalModType) {
	float t0 = 0.0;
	int n = 0;
	float dt = 1 / fs;
	float tn = t0 + (n*dt);
	const int num_samples = int(Tb / dt);

	// etap 1: x1(t) i x2(t)
	vector<float> z_xaxis_copy, z_yaxis_copy;

	z_xaxis_copy = z_xaxis;
	z_yaxis_copy = z_yaxis;

	vector<float> x_xaxis, x_xaxis_copy, x_yaxis, x_yaxis_copy;

	for (int i = 0; i < z_yaxis.size(); i++) {
		x_yaxis.push_back(z_yaxis[i] * s_n(a, f1, tn, fi));
		x_xaxis.push_back(tn);

		x_yaxis_copy.push_back(z_yaxis[i] * s_n(a, f2, tn, fi));
		x_xaxis_copy.push_back(tn);

		n++;
		tn = t0 + (n*dt);
	}

	drawSignalChart(x_xaxis, x_yaxis, 10 * num_samples, "lines", signalModType + ": demodulacja, etap 1, x_1(t)", namePrefix + "demod_etap_1_x1t", "T_b[s]", "A", 1920, 800);
	drawSignalChart(x_xaxis_copy, x_yaxis_copy, 10 * num_samples, "lines", signalModType + ": demodulacja, etap 1, x_2(t)", namePrefix + "demod_etap_1_x2t", "T_b[s]", "A", 1920, 800);

	// etap 2: p(t), calka
	vector<float> p_xaxis, p_xaxis_copy, p_yaxis, p_yaxis_copy;

	n = 0;
	tn = t0 + (n*dt);

	float sum = 0.0, sum_copy = 0.0, t = 0.0;

	for (int i = 0; i < x_yaxis.size(); i++) {
		sum += x_yaxis[i] / dt;
		sum_copy += x_yaxis_copy[i] / dt;

		p_yaxis.push_back(sum);
		p_xaxis.push_back(tn);

		p_yaxis_copy.push_back(sum_copy);
		p_xaxis_copy.push_back(tn);

		n++;
		tn = t0 + (n*dt);
		t += dt;

		if (t > Tb) {
			t = 0.0;
			sum = 0.0;
			sum_copy = 0.0;
		}
	}

	for (int i = 0; i < p_yaxis.size(); i++) {
		p_yaxis[i] = - p_yaxis[i] + p_yaxis_copy[i];
	}

	drawSignalChart(p_xaxis, p_yaxis, 10 * num_samples, "lines", signalModType + ": demodulacja, etap 2, p(t)", namePrefix + "demod_etap_2_pt", "T_b[s]", "x(t)dt", 1920, 800);

	// etap 3: m'(t)
	n = 0;
	tn = t0 + (n*dt);

	for (int i = 0; i < p_yaxis.size(); i++) {
		if (p_yaxis[i] >= h)
			mprim_yaxis.push_back(1);
		else
			mprim_yaxis.push_back(0);

		mprim_xaxis.push_back(tn);
		n++;
		tn = t0 + (n*dt);
	}

	float ymin = *min_element(mprim_yaxis.begin(), mprim_yaxis.end());
	float ymax = *max_element(mprim_yaxis.begin(), mprim_yaxis.end());

	drawSignalChart(mprim_xaxis, mprim_yaxis, 10 * num_samples, "lines", signalModType + ": demodulacja, etap 3, m'(t), h = " + to_string(int(h)), namePrefix + "demod_etap_3_mprimt", "T_b[s]", "A", 1920, 800, 0, ymax * 0.1);

	drawSignalMultiChart(
		z_xaxis, x_xaxis, x_xaxis_copy, p_xaxis, mprim_xaxis,
		z_yaxis, x_yaxis, x_yaxis_copy, p_yaxis, mprim_yaxis,
		10 * num_samples,
		"lines",
		signalModType + ": demodulacja",
		signalModType,
		signalModType + ": demodulacja, etap 1, x_1(t)",
		signalModType + ": demodulacja, etap 1, x_2(t)",
		signalModType + ": demodulacja, etap 2, p(t)",
		signalModType + ": demodulacja, etap 3, m'(t)",
		namePrefix + "demod_sygnaly",
		namePrefix.substr(0, namePrefix.length() - 1),
		namePrefix + "demod_etap_1_x1t",
		namePrefix + "demod_etap_1_x2t",
		namePrefix + "demod_etap_2_pt",
		namePrefix + "demod_etap_3_mprimt",
		"T_b[s]", "A", "T_b[s]", "x(t)dt",
		1920, 1920
	);
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


void zad1_2(vector<float> inf_xaxis, vector<float> inf_yaxis, int _N, float f, float fs) {
	float a1 = 0.0, a2 = 1.0; // ASK
	float f1 = (_N + 1) / Tb, f2 = (_N + 2) / Tb; // FSK
	float fi0 = 0, fi1 = M_PI; // PSK

	float dt = 1 / fs;
	const int num_samples = int(Tb / dt);

	string namePrefix = "LAB_06_ZAD_1_";

	// sygnal z_A(t), kluczowanie amplitudy
	// ----------------------------------------
		vector<float> z_A_yaxis, z_A_xaxis;
		generateAmplModSignal(z_A_xaxis, z_A_yaxis, inf_yaxis, a1, a2, Tb, f, fs, fi);

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
		drawSignalChart(z_A_xaxis, z_A_yaxis, 10 * num_samples, "lines", title_A.str(), filename_A, "T_b[s]", "A", 1920, 800);

	// sygnal z_F(t), kluczowanie czestotliwosci
	// ---------------------------------------------
		vector<float> z_F_yaxis, z_F_xaxis;
		generateFreqModSignal(z_F_xaxis, z_F_yaxis, inf_yaxis, a, f1, f2, Tb, fs, fi);

		// wykres sygnalu zmodulowanego z_F(t)
		stringstream title_B;

		title_B << std::fixed << std::setprecision(2)
			<< "z_F(t): sygnal zmodulowany czestotliwoscia (FSK)"
			<< ", A = " << a
			<< ", f_1 = " << f1 << " Hz"
			<< ", f_2 = " << f2 << " Hz"
			<< ", N = " << _N;

		string filename_B = namePrefix + "sygnal_z_Ft";
		drawSignalChart(z_F_xaxis, z_F_yaxis, 10 * num_samples, "lines", title_B.str(), filename_B, "T_b[s]", "A", 1920, 800);

	// sygnal z_P(t), kluczowanie fazy
	// -----------------------------------
		vector<float> z_P_yaxis, z_P_xaxis;
		generatePhaseModSignal(z_P_xaxis, z_P_yaxis, inf_yaxis, a, fi0, fi1, Tb, f, fs);

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
		drawSignalChart(z_P_xaxis, z_P_yaxis, 10 * num_samples, "lines", title_C.str(), filename_C, "T_b[s]", "A", 1920, 800);

	// sygnaly zdemodulowane i ich wykresy
	// ---------------------------------------
		float h = 2;

		vector<float> z_A_mprim_xaxis, z_A_mprim_yaxis;
		vector<float> z_F_mprim_xaxis, z_F_mprim_yaxis;
		vector<float> z_P_mprim_xaxis, z_P_mprim_yaxis;

		demodulateAmplPhaseModSignal(z_A_mprim_xaxis, z_A_mprim_yaxis, z_A_xaxis, z_A_yaxis, Tb, h, a2, f, fs, fi, s_n, namePrefix + "sygnal_z_At_", "z_A(t)");
		demodulateFreqModSignal(z_F_mprim_xaxis, z_F_mprim_yaxis, z_F_xaxis, z_F_yaxis, Tb, h, a, f1, f2, fs, fi, s_n, namePrefix + "sygnal_zFt_", "z_F(t)");
		demodulateAmplPhaseModSignal(z_P_mprim_xaxis, z_P_mprim_yaxis, z_P_xaxis, z_P_yaxis, Tb, h, a, f, fs, fi1, s_n, namePrefix + "sygnal_zPt_", "z_P(t)");
}

int main() {
	int _N = 2;
	float f = _N * pow(Tb, -1), fs = 120 * f;

	float dt = 1 / fs;
	const int num_samples = int(Tb / dt);

	vector<bitset<8>> s1 = strToBinStream("KOT", "littleEndian");

	// sygnal informacyjny
	// -----------------------
		vector<float> inf_xaxis, inf_yaxis;
		generateInfSignal(inf_xaxis, inf_yaxis, s1, Tb, fs);

		float ymin = *min_element(inf_yaxis.begin(), inf_yaxis.end());
		float ymax = *max_element(inf_yaxis.begin(), inf_yaxis.end());

		// wykres sygnalu inf.
		stringstream title;

		title << std::fixed << std::setprecision(2)
			<< "sygnal informacyjny"
			<< ", T_b = " << Tb << " s"
			<< ", f = " << f << " Hz"
			<< ", fs = " << fs << " Hz";

		string filename = "LAB_06_sygnal_inf";
		drawSignalChart(inf_xaxis, inf_yaxis, 10 * num_samples, "lines", title.str(), filename, "T_b[s]", "A", 1920, 800, 0, ymax * 0.1);

	// zad. 1. i 2.
	// -----------
		zad1_2(inf_xaxis, inf_yaxis, _N, f, fs);


	getchar();
	return 0;
}
