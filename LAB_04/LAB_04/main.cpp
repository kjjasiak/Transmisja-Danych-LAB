#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <complex>
#include <algorithm>
#include "chart_utils.h"

using namespace chart_utils;
using std::cout;
using std::endl;
using std::string;
using std::to_string;
using std::stringstream;
using std::vector;
using std::complex;

// stale uzywane w calosci programu
const int _A = 8, _B = 4, _C = 3, _BA = 48, _AB = 84;
float ampl_m = 1.0;
int f = _B, f_m = _B, f_n = f_m * 10, fs = 1000;
float T = 1, Ts = 1.0 / fs, tN = T;


vector<float> getAmplSpectrum(vector<complex<float>> r, int N) {
	vector<float> M;

	for (int k = 0; k <= N - 1; k++) {
		float kTmp;

		kTmp = sqrt(real(r[k]) * real(r[k]) + imag(r[k]) * imag(r[k]));
		kTmp = (kTmp * 2) / N;

		M.push_back(kTmp);
	}

	return M;
}

vector<float> getAmplSpectrumBis(vector<float> M, float threshold, int N) {
	vector<float> MBis;

	for (int k = 0; k <= N - 1; k++) {
		if (M[k] < threshold)
			MBis.push_back(0);
		else
			MBis.push_back(10.0 * log10(M[k]));
	}

	return MBis;
}

// skala czestotliwosci
vector<float> getFreqScale(int fs, int N) {
	vector<float> f_k;

	for (int k = 0; k <= N - 1; k++) {
		f_k.push_back(k * (float(fs) / N));
	}

	return f_k;
}

// dft
vector<complex<float>> dft(vector<float> xn, int N) {
	vector<complex<float>> xk;

	for (int k = 0; k <= N - 1; k++) {
		complex<float> kTmp;

		for (int n = 0; n <= N - 1; n++) {
			complex<float> wNPow(cos((-2 * M_PI*k*n) / N), sin((-2 * M_PI*k*n) / N));

			kTmp += xn[n] * wNPow;
		}

		xk.push_back(kTmp);
	}

	return xk;
}

// idft
vector<float> idft(vector<complex<float>> xk, int N) {
	vector<float> xn;

	for (int n = 0; n <= N - 1; n++) {
		float nTmp = 0;

		for (int k = 0; k <= N - 1; k++) {
			complex<float> wNPow(cos((2 * M_PI*k*n) / N), sin((2 * M_PI*k*n) / N));

			nTmp += real(xk[k]) * real(wNPow) - imag(xk[k]) * imag(wNPow);
		}

		nTmp = nTmp / N;

		xn.push_back(nTmp);
	}

	return xn;
}

// szerokosc pasma sygnalu
float getBandWidth(vector<float> x_dB, vector<float> f, int N) {
	float fmin = f.back();
	float fmax = f[0];

	for (int i = 0; i < N / 2; i++) {
		if (x_dB[i] >= -3 && x_dB[i] != 0) {
			if (f[i] < fmin)
				fmin = f[i];

			if (f[i] > fmax)
				fmax = f[i];
		}
	}

	return (fmax - fmin >= 0) ? fmax - fmin : 0;
}


float fn_m(float t) {
	return ampl_m * sin(2 * M_PI * f_m * t);
}

float amplModulation(float t, float k_A, int f_n, float(*fn)(float)) {
	return (k_A * fn(t) + 1) * cos(2 * M_PI * f_n * t);
}

float phaseModulation(float t, float k_P, int f_n, float(*fn)(float)) {
	return cos(2 * M_PI * f_n * t + k_P * fn(t));
}

void generateSignal(vector<float> &t, vector<float> &x, float &t0, const float &tN, float &dt, const int &N, float(*fn)(float)) {
	int n = 0;
	float tn = t0 + (n*dt);

	while (tn < tN) {
		float xn = fn(tn);

		x.push_back(xn);
		t.push_back(tn);

		n++;
		tn = t0 + (n*dt);
	}
}

void generateAmplModSignal(vector<float> &t, vector<float> &x, float &t0, const float &tN, float &dt, const int &N, float &k_A, int &f_n, float(*fn)(float, float, int, float(*fn_m)(float))) {
	int n = 0;
	float tn = t0 + (n*dt);

	while (tn < tN) {
		float xn = fn(tn, k_A, f_n, fn_m);

		x.push_back(xn);
		t.push_back(tn);

		n++;
		tn = t0 + (n*dt);
	}
}

void generatePhaseModSignal(vector<float> &t, vector<float> &x, float &t0, const float &tN, float &dt, const int &N, float &k_P, int &f_n, float(*fn)(float, float, int, float(*fn_m)(float))) {
	int n = 0;
	float tn = t0 + (n*dt);

	while (tn < tN) {
		float xn = fn(tn, k_P, f_n, fn_m);

		x.push_back(xn);
		t.push_back(tn);

		n++;
		tn = t0 + (n*dt);
	}
}
void drawSignalChart(vector<float> &t, vector<float> &x, int N, string plotType, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << t[i] << ";" << x[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());

	if (plotType == "dots")
		drawChartWithDots(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);
	else
		drawChart(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSpectrumChart(vector<float> &t, vector<float> &x, int N, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << t[i] << ";" << x[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSpectrumPartChart(vector<float> &t, vector<float> &x, int N, float threshold, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;

	for (int i = 0; i < N / 2; i++) {
		if (x[i] < threshold) {
			st << t[i] << ";" << 0 << ";\n";
		}
		else {
			st << t[i] << ";" << x[i] << ";\n";
		}
	}

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSpectrumPartChart(vector<float> &t, vector<float> &x, int N, float threshold, string title, string filename, string xlabel, string ylabel, int width, int height, int freqCap) {
	stringstream st;

	for (int i = 0; i < N / 2; i++) {
		if (x[i] < threshold) {
			st << t[i] << ";" << 0 << ";\n";
		}
		else {
			st << t[i] << ";" << x[i] << ";\n";
		}

		if (t[i] > freqCap)
			break;
	}

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSpectrumPartdBChart(vector<float> &t, vector<float> &x, int N, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;
	int counter = 0;

	for (int i = 0; i < N / 2; i++)
		st << t[i] << ";" << x[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void sygnalInf() {
	// sygnal informacyjny m(t) (ton prosty)
	int n = 0;
	float t0 = 0.0, dt = Ts, tn = t0 + (n*dt);
	const int N = ceil(T*fs);

	vector<float> x, t;
	stringstream st;

	generateSignal(t, x, t0, tN, dt, N, fn_m);

	// wykres sygnalu tonu prostego
	string title = "m(t): sygnal informacyjny, t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f_m = " + to_string(int(f_m)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	string filename = "LAB_04_sygnal_inf_mt";
	drawSignalChart(t, x, N, "dots", title, filename, "t[s]", "A", 1920, 1080);
}

void sygnalA() {
	// zad. 1
	// -------------------------

	const int N = ceil(T*fs);

	int n = 0;
	float t0 = 0.0, dt = Ts, tN = 1, tn = t0 + (n*dt);

	float k_A = 0.5;
	float k_P = 1.5;

	vector<float> x_A, t_A, x_P, t_P;

	generateAmplModSignal(t_A, x_A, t0, tN, dt, N, k_A, f_n, amplModulation);
	generatePhaseModSignal(t_P, x_P, t0, tN, dt, N, k_P, f_n, phaseModulation);

	// wykres sygnalu zmodulowanego z_A(t)
	string title = "Sygnal a), z_A(t): sygnal zmodulowany amplituda, k_A = " + to_string(k_A) + ", t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f_n = " + to_string(int(f_n)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	string filename = "LAB_04_sygnal_a_z_At";
	drawSignalChart(t_A, x_A, N, "lines", title, filename, "t[s]", "A", 1920, 1080);

	// wykres sygnalu zmodulowanego z_P(t)
	title = "Sygnal a), z_P(t): sygnal zmodulowany faza, k_P = " + to_string(k_P) + ", t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f_n = " + to_string(int(f_n)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	filename = "LAB_04_sygnal_a_z_Pt";
	drawSignalChart(t_P, x_P, N, "lines", title, filename, "t[s]", "A", 1920, 1080);

	// zad. 2
	// -------------------------

	// sygnal z_A(t)

		// dft
		vector<complex<float>> x_A_AfterDFT = dft(x_A, N);

		// obliczanie widma amplitudowego
		vector<float> x_A_AmplSpectrum = getAmplSpectrum(x_A_AfterDFT, N);

		// skala decybelowa
		float threshold = *std::max_element(x_A_AmplSpectrum.begin(), x_A_AmplSpectrum.end()) / 10000;
		vector<float> x_A_AmplSpectrumBis = getAmplSpectrumBis(x_A_AmplSpectrum, threshold, N);

		// skala czestotliwosci
		vector<float> f_k = getFreqScale(fs, N);

		// wykres widma amplitudowego (fragment)
		title = "Sygnal a), z_A(t): fragment widma ampl. w dziedzinie czest.";
		filename = "LAB_04_sygnal_a_z_At_widmo";
		drawSpectrumPartChart(f_k, x_A_AmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080);

		// wykres widma amplitudowego - skala decybelowa (fragment)
		title = "Sygnal a), z_A(t): fragment widma ampl. w dziedzinie czest. - skala dB";
		filename = "LAB_04_sygnal_a_z_At_widmo_dB";
		drawSpectrumPartdBChart(f_k, x_A_AmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);

		// sygnal a, z_A(t) => szerokosc pasma: 0 Hz
		float bandwidth = getBandWidth(x_A_AmplSpectrumBis, f_k, N);
		cout << "sygnal a, z_A(t) => szerokosc pasma: " << bandwidth << " Hz" << endl;

	// sygnal z_P(t)

		// dft
		vector<complex<float>> x_P_AfterDFT = dft(x_P, N);

		// obliczanie widma amplitudowego
		vector<float> x_P_AmplSpectrum = getAmplSpectrum(x_P_AfterDFT, N);

		// skala decybelowa
		threshold = *std::max_element(x_P_AmplSpectrum.begin(), x_P_AmplSpectrum.end()) / 10000;
		vector<float> x_P_AmplSpectrumBis = getAmplSpectrumBis(x_P_AmplSpectrum, threshold, N);

		// skala czestotliwosci
		f_k = getFreqScale(fs, N);

		// wykres widma amplitudowego (fragment)
		title = "Sygnal a), z_P(t): fragment widma ampl. w dziedzinie czest.";
		filename = "LAB_04_sygnal_a_z_Pt_widmo";
		drawSpectrumPartChart(f_k, x_P_AmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080);

		// wykres widma amplitudowego - skala decybelowa (fragment)
		title = "Sygnal a), z_P(t): fragment widma ampl. w dziedzinie czest. - skala dB";
		filename = "LAB_04_sygnal_a_z_Pt_widmo_dB";
		drawSpectrumPartdBChart(f_k, x_P_AmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);

		// sygnal a, z_P(t) => szerokosc pasma: 8 Hz
		bandwidth = getBandWidth(x_P_AmplSpectrumBis, f_k, N);
		cout << "sygnal a, z_P(t) => szerokosc pasma: " << bandwidth << " Hz" << endl;
}

void sygnalB() {
	// zad. 1
	// -------------------------

	const int N = ceil(T*fs);

	int n = 0;
	float t0 = 0.0, dt = Ts, tN = 1, tn = t0 + (n*dt);

	float k_A = 8;
	float k_P = 2.8;

	vector<float> x_A, t_A, x_P, t_P;

	generateAmplModSignal(t_A, x_A, t0, tN, dt, N, k_A, f_n, amplModulation);
	generatePhaseModSignal(t_P, x_P, t0, tN, dt, N, k_P, f_n, phaseModulation);

	// wykres sygnalu zmodulowanego z_A(t)
	string title = "Sygnal b), z_A(t): sygnal zmodulowany amplituda, k_A = " + to_string(k_A) + ", t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f_n = " + to_string(int(f_n)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	string filename = "LAB_04_sygnal_b_z_At";
	drawSignalChart(t_A, x_A, N, "lines", title, filename, "t[s]", "A", 1920, 1080);

	// wykres sygnalu zmodulowanego z_P(t)
	title = "Sygnal b), z_P(t): sygnal zmodulowany faza, k_P = " + to_string(k_P) + ", t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f_n = " + to_string(int(f_n)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	filename = "LAB_04_sygnal_b_z_Pt";
	drawSignalChart(t_P, x_P, N, "lines", title, filename, "t[s]", "A", 1920, 1080);

	// zad. 2
	// -------------------------

	// sygnal z_A(t)

		// dft
		vector<complex<float>> x_A_AfterDFT = dft(x_A, N);

		// obliczanie widma amplitudowego
		vector<float> x_A_AmplSpectrum = getAmplSpectrum(x_A_AfterDFT, N);

		// skala decybelowa
		float threshold = *std::max_element(x_A_AmplSpectrum.begin(), x_A_AmplSpectrum.end()) / 10000;
		vector<float> x_A_AmplSpectrumBis = getAmplSpectrumBis(x_A_AmplSpectrum, threshold, N);

		// skala czestotliwosci
		vector<float> f_k = getFreqScale(fs, N);

		// wykres widma amplitudowego (fragment)
		title = "Sygnal b), z_A(t): fragment widma ampl. w dziedzinie czest.";
		filename = "LAB_04_sygnal_b_z_At_widmo";
		drawSpectrumPartChart(f_k, x_A_AmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080);

		// wykres widma amplitudowego - skala decybelowa (fragment)
		title = "Sygnal b), z_A(t): fragment widma ampl. w dziedzinie czest. - skala dB";
		filename = "LAB_04_sygnal_b_z_At_widmo_dB";
		drawSpectrumPartdBChart(f_k, x_A_AmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);

		// sygnal b, z_A(t) => szerokosc pasma: 8 Hz
		float bandwidth = getBandWidth(x_A_AmplSpectrumBis, f_k, N);
		cout << "sygnal b, z_A(t) => szerokosc pasma: " << bandwidth << " Hz" << endl;

	// sygnal z_P(t)

		// dft
		vector<complex<float>> x_P_AfterDFT = dft(x_P, N);

		// obliczanie widma amplitudowego
		vector<float> x_P_AmplSpectrum = getAmplSpectrum(x_P_AfterDFT, N);

		// skala decybelowa
		threshold = *std::max_element(x_P_AmplSpectrum.begin(), x_P_AmplSpectrum.end()) / 10000;
		vector<float> x_P_AmplSpectrumBis = getAmplSpectrumBis(x_P_AmplSpectrum, threshold, N);

		// skala czestotliwosci
		f_k = getFreqScale(fs, N);

		// wykres widma amplitudowego (fragment)
		title = "Sygnal b), z_P(t): fragment widma ampl. w dziedzinie czest.";
		filename = "LAB_04_sygnal_b_z_Pt_widmo";
		drawSpectrumPartChart(f_k, x_P_AmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080);

		// wykres widma amplitudowego - skala decybelowa (fragment)
		title = "Sygnal b), z_P(t): fragment widma ampl. w dziedzinie czest. - skala dB";
		filename = "LAB_04_sygnal_b_z_Pt_widmo_dB";
		drawSpectrumPartdBChart(f_k, x_P_AmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);

		// sygnal b, z_P(t) => szerokosc pasma: 0 Hz
		bandwidth = getBandWidth(x_P_AmplSpectrumBis, f_k, N);
		cout << "sygnal b, z_P(t) => szerokosc pasma: " << bandwidth << " Hz" << endl;
}

void sygnalC() {
	// zad. 1
	// -------------------------

	const int N = ceil(T*fs);

	int n = 0;
	float t0 = 0.0, dt = Ts, tN = 1, tn = t0 + (n*dt);

	float k_A = _BA + 2;
	float k_P = _AB + 1;

	vector<float> x_A, t_A, x_P, t_P;

	generateAmplModSignal(t_A, x_A, t0, tN, dt, N, k_A, f_n, amplModulation);
	generatePhaseModSignal(t_P, x_P, t0, tN, dt, N, k_P, f_n, phaseModulation);

	// wykres sygnalu zmodulowanego z_A(t)
	string title = "Sygnal c), z_A(t): sygnal zmodulowany amplituda, k_A = " + to_string(k_A) + ", t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f_n = " + to_string(int(f_n)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	string filename = "LAB_04_sygnal_c_z_At";
	drawSignalChart(t_A, x_A, N, "lines", title, filename, "t[s]", "A", 1920, 1080);

	// wykres sygnalu zmodulowanego z_P(t)
	title = "Sygnal c), z_P(t): sygnal zmodulowany faza, k_P = " + to_string(k_P) + ", t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f_n = " + to_string(int(f_n)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	filename = "LAB_04_sygnal_c_z_Pt";
	drawSignalChart(t_P, x_P, N, "lines", title, filename, "t[s]", "A", 1920, 1080);

	// zad. 2
	// -------------------------

	// sygnal z_A(t)

		// dft
		vector<complex<float>> x_A_AfterDFT = dft(x_A, N);

		// obliczanie widma amplitudowego
		vector<float> x_A_AmplSpectrum = getAmplSpectrum(x_A_AfterDFT, N);

		// skala decybelowa
		float threshold = *std::max_element(x_A_AmplSpectrum.begin(), x_A_AmplSpectrum.end()) / 10000;
		vector<float> x_A_AmplSpectrumBis = getAmplSpectrumBis(x_A_AmplSpectrum, threshold, N);

		// skala czestotliwosci
		vector<float> f_k = getFreqScale(fs, N);

		// wykres widma amplitudowego (fragment)
		title = "Sygnal c), z_A(t): fragment widma ampl. w dziedzinie czest.";
		filename = "LAB_04_sygnal_c_z_At_widmo";
		drawSpectrumPartChart(f_k, x_A_AmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080);

		// wykres widma amplitudowego - skala decybelowa (fragment)
		title = "Sygnal c), z_A(t): fragment widma ampl. w dziedzinie czest. - skala dB";
		filename = "LAB_04_sygnal_c_z_At_widmo_dB";
		drawSpectrumPartdBChart(f_k, x_A_AmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);

		// sygnal c, z_A(t) => szerokosc pasma: 8 Hz
		float bandwidth = getBandWidth(x_A_AmplSpectrumBis, f_k, N);
		cout << "sygnal c, z_A(t) => szerokosc pasma: " << bandwidth << " Hz" << endl;

	// sygnal z_P(t)

		// dft
		vector<complex<float>> x_P_AfterDFT = dft(x_P, N);

		// obliczanie widma amplitudowego
		vector<float> x_P_AmplSpectrum = getAmplSpectrum(x_P_AfterDFT, N);

		// skala decybelowa
		threshold = *std::max_element(x_P_AmplSpectrum.begin(), x_P_AmplSpectrum.end()) / 10000;
		vector<float> x_P_AmplSpectrumBis = getAmplSpectrumBis(x_P_AmplSpectrum, threshold, N);

		// skala czestotliwosci
		f_k = getFreqScale(fs, N);

		// wykres widma amplitudowego (fragment)
		title = "Sygnal c), z_P(t): fragment widma ampl. w dziedzinie czest.";
		filename = "LAB_04_sygnal_c_z_Pt_widmo";
		drawSpectrumPartChart(f_k, x_P_AmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080);

		// wykres widma amplitudowego - skala decybelowa (fragment)
		title = "Sygnal c), z_P(t): fragment widma ampl. w dziedzinie czest. - skala dB";
		filename = "LAB_04_sygnal_c_z_Pt_widmo_dB";
		drawSpectrumPartdBChart(f_k, x_P_AmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);

		// sygnal c, z_P(t) => szerokosc pasma: 0 Hz
		bandwidth = getBandWidth(x_P_AmplSpectrumBis, f_k, N);
		cout << "sygnal c, z_P(t) => szerokosc pasma: " << bandwidth << " Hz" << endl;
}

int main() {
	sygnalInf();
	sygnalA();
	sygnalB();
	sygnalC();

	getchar();
	return 0;
}