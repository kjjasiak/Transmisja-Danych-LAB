#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <bitset>
#include <cmath>
#include <complex>

using std::vector;
using std::bitset;
using std::complex;

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

void generateAmplModSignal(vector<float> &t, vector<float> &x, vector<float> &m, float &a1, float &a2, float &Tb, float f, float fs, float fi) {
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

void generateFreqModSignal(vector<float> &t, vector<float> &x, vector<float> &m, float &a, float &f1, float &f2, float &Tb, float fs, float fi) {
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
