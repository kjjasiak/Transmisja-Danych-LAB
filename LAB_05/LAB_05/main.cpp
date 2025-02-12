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
float getBandWidth(vector<float> x_dB, vector<float> f, int N, float dBlevel) {
	float fmin = f.back();
	float fmax = f[0];

	for (int i = 0; i < N / 2; i++) {
		if (x_dB[i] >= dBlevel && x_dB[i] != 0) {
			if (f[i] < fmin)
				fmin = f[i];

			if (f[i] > fmax)
				fmax = f[i];
		}
	}

	return (fmax - fmin >= 0) ? fmax - fmin : 0;
}

void generateInfSignal(vector<float> &t, vector<float> &x, vector<bitset<8>> &s, float &Tb, float fs) {
	float t0 = 0.0;
	int n = 0;
	float dt = 1 / fs;
	float tn = t0 + (n*dt);

	const int num_samples = int(Tb / dt);

	for (int i = 0; i < s.size(); i++) {
		for (int j = 0; j < 8; j++) {
			for (int k = 0; k < num_samples; k++) {
				x.push_back(s[i][j]);
				t.push_back(tn);
				n++;
				tn = t0 + (n*dt);
			}
		}
	}
}

void generateAmplModSignal(vector<float> &t, vector<float> &x, vector<float> &m, float &a1, float &a2, float &Tb, float f, float fs) {
	float t0 = 0.0;
	int n = 0;
	float dt = 1 / fs;
	float tn = t0 + (n*dt);

	for (int i = 0; i < m.size(); i++) {
		if (m[i] == 1) {
			x.push_back(a2 * sin(2 * M_PI*f*tn + fi));
		}
		else {
			x.push_back(a1 * sin(2 * M_PI*f*tn + fi));
		}

		t.push_back(tn);
		n++;
		tn = t0 + (n*dt);
	}
}

void generateFreqModSignal(vector<float> &t, vector<float> &x, vector<float> &m, float &a, float &f1, float &f2, float &Tb, float fs) {
	float t0 = 0.0;
	int n = 0;
	float dt = 1 / fs;
	float tn = t0 + (n*dt);

	for (int i = 0; i < m.size(); i++) {
		if (m[i] == 1) {
			x.push_back(a * sin(2 * M_PI*f2*tn + fi));
		}
		else {
			x.push_back(a * sin(2 * M_PI*f1*tn + fi));
		}

		t.push_back(tn);
		n++;
		tn = t0 + (n*dt);
	}
}

void generatePhaseModSignal(vector<float> &t, vector<float> &x, vector<float> &m, float &a, float &fi0, float &fi1, float &Tb, float f, float fs) {
	float t0 = 0.0;
	int n = 0;
	float dt = 1 / fs;
	float tn = t0 + (n*dt);

	for (int i = 0; i < m.size(); i++) {
		if (m[i] == 1) {
			x.push_back(a * sin(2 * M_PI*f*tn + fi1));
		}
		else {
			x.push_back(a * sin(2 * M_PI*f*tn + fi0));
		}

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

void drawSignalChart(vector<float> &t, vector<float> &x, int N, string plotType, string title, string filename, string xlabel, string ylabel, int width, int height, float xOffset, float yOffset) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << t[i] << ";" << x[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());

	if (plotType == "dots")
		drawChartWithDots(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);
	else
		drawChart(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height, xOffset, yOffset);

	st.str("");
}

void drawSignalMultiChart(vector<float> &t1, vector<float> &t2, vector<float> &t3, vector<float> &t4, vector<float> &x1, vector<float> &x2, vector<float> &x3, vector<float> &x4, int N, string plotType, string title, string title1, string title2, string title3, string title4, string outputFile, string filename1, string filename2, string filename3, string filename4, string xlabel, string ylabel, int width, int height) {
	stringstream st;
	float ymin = *min_element(x1.begin(), x1.end());
	float ymax = *max_element(x1.begin(), x1.end());

	// x1
	for (int i = 0; i < N; i++)
		st << t1[i] << ";" << x1[i] << ";\n";

	dataToCsv("csv/" + filename1 + ".csv", st.str());
	st.str("");

	// x2
	for (int i = 0; i < N; i++)
		st << t2[i] << ";" << x2[i] << ";\n";

	dataToCsv("csv/" + filename2 + ".csv", st.str());
	st.str("");

	// x3
	for (int i = 0; i < N; i++)
		st << t3[i] << ";" << x3[i] << ";\n";

	dataToCsv("csv/" + filename3 + ".csv", st.str());
	st.str("");

	// x4
	for (int i = 0; i < N; i++)
		st << t4[i] << ";" << x4[i] << ";\n";

	dataToCsv("csv/" + filename4 + ".csv", st.str());
	st.str("");

	drawChartMultiPlot(title, title1, title2, title3, title4, xlabel, ylabel, "charts/" + outputFile + ".png", "csv/" + filename1 + ".csv", "csv/" + filename2 + ".csv", "csv/" + filename3 + ".csv", "csv/" + filename4 + ".csv", width, height, 0, ymax * 0.1);

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

void zad2_4(vector<bitset<8>> &s1) {
	int _N = 4;
	float f = _N * pow(Tb, -1), fs = 20 * f;

	float a1 = 0.0, a2 = 1.0; // ASK
	float f1 = (_N + 1) / Tb, f2 = (_N + 2) / Tb; // FSK
	float fi0 = 0, fi1 = M_PI; // PSK

	string namePrefix = "LAB_05_ZAD_2_4_";

	// sygnal informacyjny
	// -----------------------
		vector<float> t, x;
		generateInfSignal(t, x, s1, Tb, fs);

		float ymin = *min_element(x.begin(), x.end());
		float ymax = *max_element(x.begin(), x.end());

		// wykres sygnalu inf.
		stringstream title_I;

		title_I << std::fixed << std::setprecision(2)
			<< "sygnal informacyjny"
			<< ", T_b = " << Tb << " s"
			<< ", f = " << f << " Hz"
			<< ", fs = " << fs << " Hz";

		string filename_I = namePrefix + "sygnal_inf";
		drawSignalChart(t, x, x.size(), "lines", title_I.str(), filename_I, "T_b[s]", "A", 1920, 800, 0, ymax * 0.1);

	// sygnal z_A(t), kluczowanie amplitudy
	// ----------------------------------------
		vector<float> x_A, t_A;
		generateAmplModSignal(t_A, x_A, x, a1, a2, Tb, f, fs);

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
		drawSignalChart(t_A, x_A, x_A.size(), "lines", title_A.str(), filename_A, "T_b[s]", "A", 1920, 800);

		// dft
		vector<complex<float>> x_A_AfterDFT = dft(x_A, x_A.size());

		// obliczanie widma amplitudowego
		vector<float> x_A_AmplSpectrum = getAmplSpectrum(x_A_AfterDFT, x_A_AfterDFT.size());

		// skala decybelowa
		float threshold = *std::max_element(x_A_AmplSpectrum.begin(), x_A_AmplSpectrum.end()) / 10000;
		vector<float> x_A_AmplSpectrumBis = getAmplSpectrumBis(x_A_AmplSpectrum, threshold, x_A_AmplSpectrum.size());

		// skala czestotliwosci
		vector<float> f_k = getFreqScale(fs, x_A_AmplSpectrum.size());

		// wykres widma amplitudowego (fragment)
		string title = "z_A(t): fragment widma ampl. w dziedzinie czest.";
		string filename = namePrefix + "sygnal_z_At_widmo_frag";
		drawSpectrumPartChart(f_k, x_A_AmplSpectrum, x_A_AmplSpectrum.size(), threshold, title, filename, "f[Hz]", "A", 1920, 800, 120);

		// wykres widma amplitudowego - skala decybelowa (fragment)
		title = "z_A(t): fragment widma ampl. w dziedzinie czest. - skala dB";
		filename = namePrefix + "sygnal_z_At_widmo_dB_frag";
		drawSpectrumPartdBChart(f_k, x_A_AmplSpectrumBis, x_A_AmplSpectrum.size(), title, filename, "f[Hz]", "A", 1920, 800);

		// sygnal z_A(t) => szerokosc pasma:
		// 0 Hz (poziom -3 dB)
		// 10 Hz (poziom -10 dB)
		float bandwidth = getBandWidth(x_A_AmplSpectrumBis, f_k, x_A_AmplSpectrumBis.size(), -3);
		cout << "sygnal z_A(t) => szerokosc pasma: " << bandwidth << " Hz" << endl;

	// sygnal z_F(t), kluczowanie czestotliwosci
	// ---------------------------------------------
		vector<float> x_F, t_F;
		generateFreqModSignal(t_F, x_F, x, a, f1, f2, Tb, fs);

		// wykres sygnalu zmodulowanego z_F(t)
		stringstream title_B;

		title_B << std::fixed << std::setprecision(2)
			<< "z_F(t): sygnal zmodulowany czestotliwoscia (FSK)"
			<< ", A = " << a
			<< ", f_1 = " << f1 << " Hz"
			<< ", f_2 = " << f2 << " Hz"
			<< ", N = " << _N;

		string filename_B = namePrefix + "sygnal_z_Ft";
		drawSignalChart(t_F, x_F, x_F.size(), "lines", title_B.str(), filename_B, "T_b[s]", "A", 1920, 800);

		// dft
		vector<complex<float>> x_F_AfterDFT = dft(x_F, x_F.size());

		// obliczanie widma amplitudowego
		vector<float> x_F_AmplSpectrum = getAmplSpectrum(x_F_AfterDFT, x_F_AfterDFT.size());

		// skala decybelowa
		threshold = *std::max_element(x_F_AmplSpectrum.begin(), x_F_AmplSpectrum.end()) / 10000;
		vector<float> x_F_AmplSpectrumBis = getAmplSpectrumBis(x_F_AmplSpectrum, threshold, x_F_AmplSpectrum.size());

		// skala czestotliwosci
		f_k = getFreqScale(fs, x_F_AmplSpectrum.size());

		// wykres widma amplitudowego (fragment)
		title = "z_F(t): fragment widma ampl. w dziedzinie czest.";
		filename = namePrefix + "sygnal_z_Ft_widmo_frag";
		drawSpectrumPartChart(f_k, x_F_AmplSpectrum, x_F_AmplSpectrum.size(), threshold, title, filename, "f[Hz]", "A", 1920, 800, 120);

		// wykres widma amplitudowego - skala decybelowa (fragment)
		title = "z_F(t): fragment widma ampl. w dziedzinie czest. - skala dB";
		filename = namePrefix + "sygnal_z_Ft_widmo_dB_frag";
		drawSpectrumPartdBChart(f_k, x_F_AmplSpectrumBis, x_F_AmplSpectrumBis.size(), title, filename, "f[Hz]", "A", 1920, 800);

		// sygnal z_F(t) => szerokosc pasma:
		// 0 Hz (poziom -3 dB)
		// 11.6667 Hz (poziom -10 dB)
		bandwidth = getBandWidth(x_F_AmplSpectrumBis, f_k, x_F_AmplSpectrumBis.size(), -3);
		cout << "sygnal z_F(t) => szerokosc pasma: " << bandwidth << " Hz" << endl;

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
		drawSignalChart(t_P, x_P, x_P.size(), "lines", title_C.str(), filename_C, "T_b[s]", "A", 1920, 800);

		// dft
		vector<complex<float>> x_P_AfterDFT = dft(x_P, x_P.size());

		// obliczanie widma amplitudowego
		vector<float> x_P_AmplSpectrum = getAmplSpectrum(x_P_AfterDFT, x_P_AfterDFT.size());

		// skala decybelowa
		threshold = *std::max_element(x_P_AmplSpectrum.begin(), x_P_AmplSpectrum.end()) / 10000;
		vector<float> x_P_AmplSpectrumBis = getAmplSpectrumBis(x_P_AmplSpectrum, threshold, x_P_AmplSpectrum.size());

		// skala czestotliwosci
		f_k = getFreqScale(fs, x_P_AmplSpectrum.size());

		// wykres widma amplitudowego (fragment)
		title = "z_P(t): fragment widma ampl. w dziedzinie czest.";
		filename = namePrefix + "sygnal_z_Pt_widmo_frag";
		drawSpectrumPartChart(f_k, x_P_AmplSpectrum, x_P_AmplSpectrum.size(), threshold, title, filename, "f[Hz]", "A", 1920, 800, 180);

		// wykres widma amplitudowego - skala decybelowa (fragment)
		title = "z_P(t): fragment widma ampl. w dziedzinie czest. - skala dB";
		filename = namePrefix + "sygnal_z_Pt_widmo_dB_frag";
		drawSpectrumPartdBChart(f_k, x_P_AmplSpectrumBis, x_P_AmplSpectrumBis.size(), title, filename, "f[Hz]", "A", 1920, 800);

		// sygnal z_P(t) => szerokosc pasma:
		// 0 Hz (poziom -3 dB)
		// 12.5 Hz (poziom -10 dB)
		bandwidth = getBandWidth(x_P_AmplSpectrumBis, f_k, x_P_AmplSpectrumBis.size(), -3);
		cout << "sygnal z_P(t) => szerokosc pasma: " << bandwidth << " Hz" << endl;

	// wykres zbiorczy sygnalow
	// ----------------------------
		stringstream title_inf;

		title_inf << std::fixed << std::setprecision(2)
			<< "sygnal informacyjny"
			<< ", T_b = " << Tb << " s"
			<< ", f = " << f << " Hz"
			<< ", fs = " << fs << " Hz";

		drawSignalMultiChart(t, t_A, t_F, t_P, x, x_A, x_F, x_P, x.size(), "lines", "Sygnaly: inf., z_A(t), z_F(t), z_P(t), N = " + to_string(_N), title_inf.str(), title_A.str(), title_B.str(), title_C.str(), namePrefix + "sygnaly", "LAB_05_sygnal_inf_ogr", namePrefix + "sygnal_z_At_ogr", namePrefix + "sygnal_z_Ft_ogr", namePrefix + "sygnal_z_Pt_ogr", "T_b[s]", "A", 1920, 1600);
}

void zad3(vector<bitset<8>> &s1) {
	int _N = 2;
	float f = _N * pow(Tb, -1), fs = 40 * f;

	float a1 = 0.0, a2 = 1.0; // ASK
	float f1 = (_N + 1) / Tb, f2 = (_N + 2) / Tb; // FSK
	float fi0 = 0, fi1 = M_PI; // PSK

	float dt = 1 / fs;
	const int num_samples = int(Tb / dt);

	string namePrefix = "LAB_05_ZAD_3_";

	// sygnal informacyjny
	// -----------------------
		vector<float> t, x;
		generateInfSignal(t, x, s1, Tb, fs);

		float ymin = *min_element(x.begin(), x.end());
		float ymax = *max_element(x.begin(), x.end());

		// wykres sygnalu inf.
			stringstream title_I;

			title_I << std::fixed << std::setprecision(2)
				<< "sygnal informacyjny"
				<< ", T_b = " << Tb << " s"
				<< ", f = " << f << " Hz"
				<< ", fs = " << fs << " Hz";

			string filename_I = namePrefix + "sygnal_inf";
			drawSignalChart(t, x, 10 * num_samples, "lines", title_I.str(), filename_I, "T_b[s]", "A", 1920, 800, 0, ymax * 0.1);

	// sygnal z_A(t), kluczowanie amplitudy
	// ----------------------------------------
		vector<float> x_A, t_A;
		generateAmplModSignal(t_A, x_A, x, a1, a2, Tb, f, fs);

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
		generateFreqModSignal(t_F, x_F, x, a, f1, f2, Tb, fs);

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

		drawSignalMultiChart(t, t_A, t_F, t_P, x, x_A, x_F, x_P, 10 * num_samples, "lines", "Sygnaly: inf., z_A(t), z_F(t), z_P(t), N = " + to_string(_N), title_I.str(), title_A.str(), title_B.str(), title_C.str(), namePrefix + "sygnaly", namePrefix + "sygnal_inf", namePrefix + "sygnal_z_At", namePrefix + "sygnal_z_Ft", namePrefix + "sygnal_z_Pt", "T_b[s]", "A", 1920, 1600);
}

int main() {
	// zad. 1.
	// -----------

		// wynik: 010010110100111101010100
		vector<bitset<8>> s1 = strToBinStream("KOT", "littleEndian");

		// wynik: 110100101111001000101010
		vector<bitset<8>> s2 = strToBinStream("KOT", "bigEndian");

	// zad. 2. i 4.
	// ----------------

		// sygnal inf., sygnaly z_A(t), z_F(t), z_P(t) i ich widma
		zad2_4(s1);

	// zad. 3.
	// -----------
		zad3(s1);

	getchar();
	return 0;
}
