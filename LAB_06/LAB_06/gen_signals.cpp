#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <bitset>
#include <cmath>

using std::vector;
using std::bitset;

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
