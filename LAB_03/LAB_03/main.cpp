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
const int _A = 8, _B = 4, _C = 3;
const float ampl = 1.0;
const float fi = _C * M_PI;

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

// generowanie i rysowanie sygnalow/widm
void generateSignal(vector<float> &t, vector<float> &x, const int &f, float &t0, const float &tN, float &dt, const int &N, float(*fn)(float, float, float, float)) {
	int n = 0;
	float tn = t0 + (n*dt);

	while (tn < tN) {
		float xn = fn(tn, ampl, f, fi);

		x.push_back(xn);
		t.push_back(tn);

		n++;
		tn = t0 + (n*dt);
	}
}

void generateSignalLab01(vector<float> &t, vector<float> &x, const int &f, float &t0, const float &tN, float &dt, const int &N, float(*fn)(float)) {
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

void generateTestSignal(vector<float> &t, vector<float> &x, const int &f, float &t0, const float &tN, float &dt, const int &M) {
	int n = 0;
	float tn = t0 + (n*dt);

	while (tn <= tN) {
		int m = 1;

		float xn = 0;

		while (m < M) {
			xn = xn + sin(m * 2 * M_PI*f*tn) / m;
			m = m + 2;
		}

		xn = xn * ((4 * ampl) / M_PI);
		x.push_back(xn);
		t.push_back(tn);

		n++;
		tn = t0 + (n*dt);
	}
}

void drawSignalChart(vector<float> &t, vector<float> &x, int N, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << t[i] << ";" << x[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithDots(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

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
	int counter = 0;

	for (int i = 0; i < N / 2; i++) {
		if (x[i] < threshold) {
			counter++;
			st << t[i] << ";" << 0 << ";\n";
		}
		else {
			counter = 0;
			st << t[i] << ";" << x[i] << ";\n";
		}

		if (counter == 5)
			break;
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

	for (int i = 0; i < N / 2; i++) {
		if (x[i] == 0)
			counter++;
		else
			counter = 0;

		if (counter == 5)
			break;

		st << t[i] << ";" << x[i] << ";\n";
	}

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}


// funkcje do generowania sygnalow z LAB_01
float fn_xx(float t) {
	return _A * pow(t, 2) + _B * t + _C;
}

float fn_y(float t) {
	return 2 * pow(fn_xx(t), 2) + 12 * cos(t);
}

float fn_z(float t) {
	return sin(2 * M_PI * 7 * t) * fn_xx(t) - 0.2 * log10(abs(fn_y(t)) + M_PI);
}

float fn_u(float t) {
	return sqrt(abs(fn_y(t) * fn_y(t) * fn_z(t))) - 1.8 * sin(0.4 * t * fn_z(t) * fn_xx(t));
}

float fn_v(float t) {
	if ((t >= 0) && (t < 0.22))
		return (1 - 7 * t) * sin((2 * M_PI * t * 10) / (t + 0.04));

	if ((t < 0.7) && (t >= 0.22))
		return 0.63 * t * sin(125 * t);

	if ((t >= 0.7) && (t <= 1.0))
		return pow(t, -0.662) + 0.77 * sin(8 * t);
}

float fn_p(float t) {
	float pt = 0;
	int N = 84;

	for (int n = 1; n <= N; n++)
		pt += (cos(12 * t * pow(n, 2)) + cos(16 * t * n)) / pow(n, 2);

	return pt;
}


void zad2() {
	// generowanie sygnalu
	const int f = 100, fs = 100000, M = 49;
	const float T = 0.01, Ts = 1.0 / fs;
	const int N = ceil(T*fs);

	int n = 0;
	float t0 = 0.0, dt = Ts, tN = T - dt, tn = t0 + (n*dt);

	vector<float> x, t;
	stringstream st;

	generateTestSignal(t, x, f, t0, tN, dt, M);

	// wykres sygnalu
	string title = "x(t): sygnal testowy, t <" + to_string(int(t0)) + ";" + to_string(tN+dt) + ">, f = " + to_string(int(f)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	string filename = "LAB_03_ZAD2_test_xt_01";
	drawSignalChart(t, x, N, title, filename, "t[s]", "A", 1920, 1080);

	// dft/idft
	vector<complex<float>> xAfterDFT = dft(x, N);
	vector<float> xAfterIDFT = idft(xAfterDFT, N);

	// wykres sygnalu po IDFT
	title = "x(t): sygnal X(k) => x(t) po IDFT, t <" + to_string(int(t0)) + ";" + to_string(T) + ">, f = " + to_string(int(f)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	filename = "LAB_03_ZAD2_test_xt_02_idft";
	drawSignalChart(t, xAfterIDFT, N, title, filename, "t[s]", "A", 1920, 1080);

	// obliczanie widma amplitudowego
	vector<float> xAmplSpectrum = getAmplSpectrum(xAfterDFT, N);

	// skala decybelowa
	float threshold = *std::max_element(xAmplSpectrum.begin(), xAmplSpectrum.end()) / 10000;
	vector<float> xAmplSpectrumBis = getAmplSpectrumBis(xAmplSpectrum, threshold, N);

	// skala czestotliwosci
	vector<float> f_k = getFreqScale(fs, N);

	// wykres widma amplitudowego (fragment)
	title = "x(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci";
	filename = "LAB_03_ZAD2_test_xt_04_widmo_frag";
	drawSpectrumPartChart(f_k, xAmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080);

	// wykres widma amplitudowego - skala decybelowa (fragment)
	title = "x(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci - skala dB";
	filename = "LAB_03_ZAD2_test_xt_06_widmo_dB_frag";
	drawSpectrumPartdBChart(f_k, xAmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);
}

void zad3_xt() {
	const int f = _B;
	const float fs = 42.15, T = 20, Ts = 1.0 / fs, tN = 10;
	const int N = ceil(T*fs);
	float t0 = -10, dt = Ts;
	
	int n = 0;
	float tn = t0 + (n*dt);

	vector<float> x, t;
	stringstream st;

	generateSignalLab01(t, x, f, t0, tN, dt, N, fn_xx);

	// wykres x(t)
	string title = "x(t): sygnal, t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f = " + to_string(int(f)) + " Hz, f_s = " + to_string(fs) + " Hz, n <0;" + to_string(N) + ">";
	string filename = "LAB_03_ZAD3_xt_01";
	drawSignalChart(t, x, N, title, filename, "t[s]", "A", 1920, 1080);

	// dft
	vector<complex<float>> xAfterDFT = dft(x, N);

	// obliczanie widma amplitudowego
	vector<float> xAmplSpectrum = getAmplSpectrum(xAfterDFT, N);

	// skala decybelowa
	float threshold = *std::max_element(xAmplSpectrum.begin(), xAmplSpectrum.end()) / 10000;
	vector<float> xAmplSpectrumBis = getAmplSpectrumBis(xAmplSpectrum, threshold, N);

	// skala czestotliwosci
	vector<float> f_k = getFreqScale(fs, N);

	// wykres widma amplitudowego (fragment)
	title = "x(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci";
	filename = "LAB_03_ZAD3_xt_03_widmo_frag";
	drawSpectrumPartChart(f_k, xAmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080, 5);

	// wykres widma amplitudowego - skala decybelowa (fragment)
	title = "x(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci - skala dB";
	filename = "LAB_03_ZAD3_xt_05_widmo_dB_frag";
	drawSpectrumPartdBChart(f_k, xAmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);
}

void zad3_yt() {
	const int f = _B, fs = 843;
	const float T = 1, Ts = 1.0 / fs, tN = 1;
	const int N = ceil(T*fs);
	float t0 = 0, dt = Ts;

	int n = 0;
	float tn = t0 + (n*dt);

	vector<float> x, t;
	stringstream st;

	generateSignalLab01(t, x, f, t0, tN, dt, N, fn_y);

	// wykres y(t)
	string title = "y(t): sygnal, t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f = " + to_string(int(f)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	string filename = "LAB_03_ZAD3_yt_01";
	drawSignalChart(t, x, N, title, filename, "t[s]", "A", 1920, 1080);

	// dft
	vector<complex<float>> xAfterDFT = dft(x, N);

	// obliczanie widma amplitudowego
	vector<float> xAmplSpectrum = getAmplSpectrum(xAfterDFT, N);

	// skala decybelowa
	float threshold = *std::max_element(xAmplSpectrum.begin(), xAmplSpectrum.end()) / 10000;
	vector<float> xAmplSpectrumBis = getAmplSpectrumBis(xAmplSpectrum, threshold, N);

	// skala czestotliwosci
	vector<float> f_k = getFreqScale(fs, N);

	// wykres widma amplitudowego (fragment)
	title = "y(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci";
	filename = "LAB_03_ZAD3_yt_03_widmo_frag";
	drawSpectrumPartChart(f_k, xAmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080, 100);

	// wykres widma amplitudowego - skala decybelowa (fragment)
	title = "y(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci - skala dB";
	filename = "LAB_03_ZAD3_yt_05_widmo_dB_frag";
	drawSpectrumPartdBChart(f_k, xAmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);
}

void zad3_zt() {
	const int f = _B, fs = 843;
	const float T = 1, Ts = 1.0 / fs, tN = 1;
	const int N = ceil(T*fs);
	float t0 = 0, dt = Ts;

	int n = 0;
	float tn = t0 + (n*dt);

	vector<float> x, t;
	stringstream st;

	generateSignalLab01(t, x, f, t0, tN, dt, N, fn_z);

	// wykres z(t)
	string title = "z(t): sygnal, t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f = " + to_string(int(f)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	string filename = "LAB_03_ZAD3_zt_01";
	drawSignalChart(t, x, N, title, filename, "t[s]", "A", 1920, 1080);

	// dft
	vector<complex<float>> xAfterDFT = dft(x, N);

	// obliczanie widma amplitudowego
	vector<float> xAmplSpectrum = getAmplSpectrum(xAfterDFT, N);

	// skala decybelowa
	float threshold = *std::max_element(xAmplSpectrum.begin(), xAmplSpectrum.end()) / 10000;
	vector<float> xAmplSpectrumBis = getAmplSpectrumBis(xAmplSpectrum, threshold, N);

	// skala czestotliwosci
	vector<float> f_k = getFreqScale(fs, N);

	// wykres widma amplitudowego (fragment)
	title = "z(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci";
	filename = "LAB_03_ZAD3_zt_03_widmo_frag";
	drawSpectrumPartChart(f_k, xAmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080, 100);

	// wykres widma amplitudowego - skala decybelowa (fragment)
	title = "z(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci - skala dB";
	filename = "LAB_03_ZAD3_zt_05_widmo_dB_frag";
	drawSpectrumPartdBChart(f_k, xAmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);
}

void zad3_ut() {
	const int f = _B, fs = 843;
	const float T = 1, Ts = 1.0 / fs, tN = 1;
	const int N = ceil(T*fs);
	float t0 = 0, dt = Ts;

	int n = 0;
	float tn = t0 + (n*dt);

	vector<float> x, t;
	stringstream st;

	generateSignalLab01(t, x, f, t0, tN, dt, N, fn_u);

	// wykres u(t)
	string title = "u(t): sygnal, t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f = " + to_string(int(f)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	string filename = "LAB_03_ZAD3_ut_01";
	drawSignalChart(t, x, N, title, filename, "t[s]", "A", 1920, 1080);

	// dft
	vector<complex<float>> xAfterDFT = dft(x, N);

	// obliczanie widma amplitudowego
	vector<float> xAmplSpectrum = getAmplSpectrum(xAfterDFT, N);

	// skala decybelowa
	float threshold = *std::max_element(xAmplSpectrum.begin(), xAmplSpectrum.end()) / 10000;
	vector<float> xAmplSpectrumBis = getAmplSpectrumBis(xAmplSpectrum, threshold, N);

	// skala czestotliwosci
	vector<float> f_k = getFreqScale(fs, N);

	// wykres widma amplitudowego (fragment)
	title = "u(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci";
	filename = "LAB_03_ZAD3_ut_03_widmo_frag";
	drawSpectrumPartChart(f_k, xAmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080, 100);

	// wykres widma amplitudowego - skala decybelowa (fragment)
	title = "u(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci - skala dB";
	filename = "LAB_03_ZAD3_ut_05_widmo_dB_frag";
	drawSpectrumPartdBChart(f_k, xAmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);
}

void zad3_vt() {
	const int f = _B, fs = 843;
	const float T = 1, Ts = 1.0 / fs, tN = 1;
	const int N = ceil(T*fs);
	float t0 = 0, dt = Ts;

	int n = 0;
	float tn = t0 + (n*dt);

	vector<float> x, t;
	stringstream st;

	generateSignalLab01(t, x, f, t0, tN, dt, N, fn_v);

	// wykres v(t)
	string title = "v(t): sygnal, t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f = " + to_string(int(f)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	string filename = "LAB_03_ZAD3_vt_01";
	drawSignalChart(t, x, N, title, filename, "t[s]", "A", 1920, 1080);

	// dft
	vector<complex<float>> xAfterDFT = dft(x, N);

	// obliczanie widma amplitudowego
	vector<float> xAmplSpectrum = getAmplSpectrum(xAfterDFT, N);

	// skala decybelowa
	float threshold = *std::max_element(xAmplSpectrum.begin(), xAmplSpectrum.end()) / 10000;
	vector<float> xAmplSpectrumBis = getAmplSpectrumBis(xAmplSpectrum, threshold, N);

	// skala czestotliwosci
	vector<float> f_k = getFreqScale(fs, N);

	// wykres widma amplitudowego (fragment)
	title = "v(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci";
	filename = "LAB_03_ZAD3_vt_03_widmo_frag";
	drawSpectrumPartChart(f_k, xAmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080, 200);

	// wykres widma amplitudowego - skala decybelowa (fragment)
	title = "v(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci - skala dB";
	filename = "LAB_03_ZAD3_vt_05_widmo_dB_frag";
	drawSpectrumPartdBChart(f_k, xAmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);
}

void zad3_pt() {
	const int f = _B, fs = 843;
	const float T = 1, Ts = 1.0 / fs, tN = 1;
	const int N = ceil(T*fs);
	float t0 = 0, dt = Ts;

	int n = 0;
	float tn = t0 + (n*dt);

	vector<float> x, t;
	stringstream st;

	generateSignalLab01(t, x, f, t0, tN, dt, N, fn_p);

	// wykres p(t)
	string title = "p(t): sygnal, N = 84, t <" + to_string(int(t0)) + ";" + to_string(int(tN)) + ">, f = " + to_string(int(f)) + " Hz, f_s = " + to_string(int(fs)) + " Hz, n <0;" + to_string(N) + ">";
	string filename = "LAB_03_ZAD3_pt_01";
	drawSignalChart(t, x, N, title, filename, "t[s]", "A", 1920, 1080);

	// dft
	vector<complex<float>> xAfterDFT = dft(x, N);

	// obliczanie widma amplitudowego
	vector<float> xAmplSpectrum = getAmplSpectrum(xAfterDFT, N);

	// skala decybelowa
	float threshold = *std::max_element(xAmplSpectrum.begin(), xAmplSpectrum.end()) / 10000;
	vector<float> xAmplSpectrumBis = getAmplSpectrumBis(xAmplSpectrum, threshold, N);

	// skala czestotliwosci
	vector<float> f_k = getFreqScale(fs, N);

	// wykres widma amplitudowego (fragment)
	title = "p(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci";
	filename = "LAB_03_ZAD3_pt_03_widmo_frag";
	drawSpectrumPartChart(f_k, xAmplSpectrum, N, threshold, title, filename, "f[Hz]", "A", 1920, 1080, 200);

	// wykres widma amplitudowego - skala decybelowa (fragment)
	title = "p(t): fragment widma amplitudowego sygnalu w dziedzinie czestotliwosci - skala dB";
	filename = "LAB_03_ZAD3_pt_05_widmo_dB_frag";
	drawSpectrumPartdBChart(f_k, xAmplSpectrumBis, N, title, filename, "f[Hz]", "A", 1920, 1080);
}


int main() {
	zad2();
	zad3_xt();
	zad3_yt();
	zad3_zt();
	zad3_ut();
	zad3_vt();
	zad3_pt();

	getchar();
	return 0;
}