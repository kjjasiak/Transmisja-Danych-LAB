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
float f = 10.0, Tb = 1.0 / f, fs = 20 * f;
const int num_samples = 100 * f;

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

		if (option == "littleEndian")
			reverseBitset(tmp);

		output.push_back(tmp);
	}

	for (int i = 0; i < output.size(); i++)
		cout << output[i] << " ";

	cout << endl;

	return output;
}

// generatory sygnalow
void generateClockSignal(vector<float> &t, vector<float> &x, float &f, float &fi, float &Tb, float &tN, int num_samples) {
	float t0 = 0.0;
	int n = 0;
	float dt = Tb / num_samples;
	float tn = t0 + (n*dt);
	float a = 1.0;

	cout << "dt = " << dt << endl;

	while (tn < tN) {
		for (int k = 0; k <= num_samples; k++) {
			float sample = a * sin(2 * M_PI * f * tn + fi);

			if (sample > 0)
				x.push_back(1);
			else
				x.push_back(0);

			t.push_back(tn);

			n++;
			tn = t0 + (n*dt);
		}
	}
}

void generateClockSignal2(vector<float> &t, vector<float> &x, float &Tb, float &tN, int num_samples) { // poprzedni generator CLK
	float t0 = 0.0;
	int n = 0;
	float dt = Tb / num_samples;
	float tn = t0 + (n*dt);

	cout << "dt = " << dt << endl;

	while (tn < tN) {
		for (int k = 0; k <= num_samples; k++) {
			if (k <= num_samples / 2) {
				x.push_back(1);
				t.push_back(tn);
			}
			else {
				x.push_back(0);
				t.push_back(tn);
			}

			n++;
			tn = t0 + (n*dt);
		}
	}
}

bool isClkDown(float &clk_curr, float &clk_next) {
	return ((clk_curr != clk_next) && (clk_next == 0));
}

bool isClkUp(float &clk_curr, float &clk_next) {
	return ((clk_curr != clk_next) && (clk_next == 1));
}

bool isClkChanging(float &clk_curr, float &clk_next) {
	return (clk_curr != clk_next);
}

void generateTTLSignal(vector<float> &ttl_xaxis, vector<float> &ttl_yaxis, vector<bitset<8>> &s, float &Tb, int num_samples) {
	float t0 = 0.0;
	int n = 0;
	float dt = Tb / num_samples;
	float tn = t0 + (n*dt);

	for (int i = 0; i < s.size(); i++) { // bajty
		for (int j = 7; j >= 0; j--) { // bity
			for (int k = 0; k < num_samples; k++) {
				ttl_yaxis.push_back(s[i][j]);
				ttl_xaxis.push_back(tn);
				n++;
				tn = t0 + (n*dt);
			}
		}
	}
}

void generateNRZISignal(vector<float> &nrzi_xaxis, vector<float> &nrzi_yaxis, vector<float> &clk_xaxis, vector<float> &clk_yaxis, vector<float> &ttl_xaxis, vector<float> &ttl_yaxis) {
	int state = 0;

	for (int i = 0; i < ttl_yaxis.size(); i++) {
		if (isClkDown(clk_yaxis[i], clk_yaxis[i + 1])) {
			nrzi_yaxis.push_back(state);
			nrzi_xaxis.push_back(clk_xaxis[i]);

			if (ttl_yaxis[i] == 1) {
				if (state == 0)
					state = -1;
				else
					state = -state;
			}
		}
		else {
			nrzi_yaxis.push_back(state);
			nrzi_xaxis.push_back(clk_xaxis[i]);
		}
	}
}

void generateBAMISignal(vector<float> &bami_xaxis, vector<float> &bami_yaxis, vector<float> &clk_xaxis, vector<float> &clk_yaxis, vector<float> &ttl_xaxis, vector<float> &ttl_yaxis) {
	int state = 0, prevState = 0;

	for (int i = 0; i < ttl_yaxis.size(); i++) {
		if (isClkUp(clk_yaxis[i], clk_yaxis[i + 1])) {
			if (i + 1 < ttl_yaxis.size()) {
				if (ttl_yaxis[i + 1] == 1) {
					if (state == 0 && prevState != -1)
						state = -1;
					else if (state == 0 && prevState == -1)
						state = 1;
					else
						state = -state;

					prevState = state;
				}
				else {
					state = 0;
				}
			}
		}

		bami_yaxis.push_back(state);
		bami_xaxis.push_back(clk_xaxis[i]);
	}
}

void generateManchesterSignal(vector<float> &manch_xaxis, vector<float> &manch_yaxis, vector<float> &clk_xaxis, vector<float> &clk_yaxis, vector<float> &ttl_xaxis, vector<float> &ttl_yaxis) {
	int state = 0, prevState = 0;

	int lastIndex = 0;

	for (int i = 0; i < ttl_yaxis.size(); i++) {
		if (isClkDown(clk_yaxis[i], clk_yaxis[i + 1])) {
			if (ttl_yaxis[i] == 1) {
				state = -1;
			}
			else {
				state = 1;
			}

			cout << "clkDown, state = " << state << endl;
		}

		if (((i - 1 > 0) && (i + 1 < ttl_yaxis.size())) && (isClkUp(clk_yaxis[i], clk_yaxis[i + 1]))) {
			cout << i << " " << i + 1 << endl;
			if (ttl_yaxis[i - 1] == ttl_yaxis[i + 1]) {
				state = -state;
			}

			cout << "clkUp, state = " << state << endl;
		}

		manch_yaxis.push_back(state);
		manch_xaxis.push_back(clk_xaxis[i]);

		lastIndex = i;
	}

	cout << "lastIndex = " << lastIndex << endl;
	cout << "clk_xaxis[lastIndex] = " << clk_xaxis[lastIndex] << endl;
	cout << "ttl_xaxis.size() = " << ttl_xaxis.size() << endl;
	cout << "ttl_yaxis.size() = " << ttl_yaxis.size() << endl;

	cout << "manch_xaxis.size() = " << manch_xaxis.size() << endl;
	cout << "manch_yaxis.size() = " << manch_yaxis.size() << endl;
}

void decodeTTLSignal(vector<float> &d_ttl_xaxis, vector<float> &d_ttl_yaxis, vector<float> &ttl_xaxis, vector<float> &ttl_yaxis) {
	for (int i = 0; i < ttl_yaxis.size(); i++) {
		d_ttl_yaxis.push_back(ttl_yaxis[i]);
		d_ttl_xaxis.push_back(ttl_xaxis[i]);
	}
}

void decodeNRZISignal(vector<float> &d_nrzi_xaxis, vector<float> &d_nrzi_yaxis, vector<float> &nrzi_xaxis, vector<float> &nrzi_yaxis, int num_samples) {
	int startIndex = num_samples;

	for (int i = 0; i < num_samples; i++) {
		d_nrzi_yaxis.push_back(0);
		d_nrzi_xaxis.push_back(nrzi_xaxis[i]);
	}

	for (int i = startIndex; i < nrzi_yaxis.size(); i++) {
		if (nrzi_yaxis[i] == -1)
			nrzi_yaxis[i] = 0;
	}

	int lastBitIndex = 0;
	int xorTmp = 0;

	for (int i = startIndex + (num_samples + 1); i < nrzi_yaxis.size(); i = i + (num_samples + 1)) {
		xorTmp = int(nrzi_yaxis[i - (num_samples + 1)]) ^ int(nrzi_yaxis[i]);

		for (int j = i - (num_samples + 1); j < i; j++) {
			d_nrzi_xaxis.push_back(nrzi_xaxis[j]);
			d_nrzi_yaxis.push_back(xorTmp);
			lastBitIndex = j;
		}
	}

	int toXOR = nrzi_yaxis[lastBitIndex];

	for (int l = lastBitIndex + 1; l < nrzi_yaxis.size(); l++) {
		d_nrzi_yaxis.push_back(toXOR ^ int(nrzi_yaxis[l]));
		d_nrzi_xaxis.push_back(nrzi_xaxis[l]);
	}
}

void decodeBAMISignal(vector<float> &d_bami_xaxis, vector<float> &d_bami_yaxis, vector<float> &bami_xaxis, vector<float> &bami_yaxis) {
	for (int i = 0; i < bami_yaxis.size(); i++) {
		d_bami_yaxis.push_back(abs(bami_yaxis[i]));
		d_bami_xaxis.push_back(bami_xaxis[i]);
	}
}

void decodeManchesterSignal(vector<float> &d_manch_xaxis, vector<float> &d_manch_yaxis, vector<float> &manch_xaxis, vector<float> &manch_yaxis, float &f, float &Tb, float &tN, float num_samples) {
	vector<float> s_clk_xaxis, s_clk_yaxis;
	float fi = M_PI / 2;
	int state = 0;

	generateClockSignal(s_clk_xaxis, s_clk_yaxis, f, fi, Tb, tN, num_samples);

	for (int i = 0; i < manch_yaxis.size(); i++) {
		if (isClkDown(s_clk_yaxis[i], s_clk_yaxis[i + 1])) {
			cout << "clkDown: i = " << i << ", state = " << state << endl;
			state = manch_yaxis[i];

			if (state == -1)
				state = 0;
		}

		d_manch_yaxis.push_back(~state);
		//d_manch_yaxis.push_back(state);
		d_manch_xaxis.push_back(s_clk_xaxis[i]);
	}
}

int main() {
	string namePrefix = "LAB_07_", filename = "";
	stringstream title;

	float tN = Tb * 30.0, fi = 0;

	// sygnal zegarowy (CLK)
	// -------------------------
		vector<float> clk_xaxis, clk_yaxis;
		generateClockSignal(clk_xaxis, clk_yaxis, f, fi, Tb, tN, num_samples);

		title << std::fixed << std::setprecision(2)
			<< "sygnal zegarowy (CLK)"
			<< ", f = " << f << " Hz";

		filename = namePrefix + "clk";

		drawSignalChart(clk_xaxis, clk_yaxis, 16 * num_samples, "lines", title.str(), filename, "t[s]", "A", 1920, 800, 0.1, 0, 0.1);
		title.str("");

	// sygnal TTL
	// -----------------------------
		vector<bitset<8>> s1 = strToBinStream("KOT", "bigEndian");
		//vector<bitset<8>> s1 = strToBinStream("KOT", "littleEndian");
		//vector<bitset<8>> s1 = strToBinStream("ALA", "bigEndian");

		vector<float> ttl_xaxis, ttl_yaxis;
		generateTTLSignal(ttl_xaxis, ttl_yaxis, s1, Tb, num_samples);

		//float xTics = (ttl_xaxis.back() - ttl_xaxis.front()) / 10;

		float ymin = *min_element(ttl_yaxis.begin(), ttl_yaxis.end());
		float ymax = *max_element(ttl_yaxis.begin(), ttl_yaxis.end());

		// wykres sygnalu TTL
		title << std::fixed << std::setprecision(2)
			<< "sygnal TTL"
			<< ", T_b = " << Tb << " s"
			<< ", f = " << f << " Hz"
			<< ", fs = " << fs << " Hz";

		filename = namePrefix + "ttl";
		drawSignalChart(ttl_xaxis, ttl_yaxis, 16 * num_samples, "steps", title.str(), filename, "T_b[s]", "A", 1920, 800, 0.1, 0, ymax * 0.1);
		title.str("");

	// sygnal NRZI
	// ---------------
		vector<float> nrzi_xaxis, nrzi_yaxis;

		generateNRZISignal(nrzi_xaxis, nrzi_yaxis, clk_xaxis, clk_yaxis, ttl_xaxis, ttl_yaxis);

		// wykres sygnalu NRZI
		title << std::fixed << std::setprecision(2)
			<< "sygnal NRZI"
			<< ", T_b = " << Tb << " s"
			<< ", f = " << f << " Hz"
			<< ", fs = " << fs << " Hz";

		filename = namePrefix + "nrzi";
		drawSignalChart(nrzi_xaxis, nrzi_yaxis, 16 * num_samples, "lines", title.str(), filename, "T_b[s]", "A", 1920, 800, 0.1, 0, ymax * 0.1);
		title.str("");

	// sygnal BAMI
	// ---------------
		vector<float> bami_xaxis, bami_yaxis;

		generateBAMISignal(bami_xaxis, bami_yaxis, clk_xaxis, clk_yaxis, ttl_xaxis, ttl_yaxis);

		// wykres sygnalu BAMI
		title << std::fixed << std::setprecision(2)
			<< "sygnal BAMI"
			<< ", T_b = " << Tb << " s"
			<< ", f = " << f << " Hz"
			<< ", fs = " << fs << " Hz";

		filename = namePrefix + "bami";
		drawSignalChart(bami_xaxis, bami_yaxis, 16 * num_samples, "lines", title.str(), filename, "T_b[s]", "A", 1920, 800, 0.1, 0, ymax * 0.1);
		title.str("");

		cout << "clk size: " << clk_yaxis.size() << endl;
		cout << "ttl size: " << ttl_yaxis.size() << endl;
		cout << "nrzi size: " << nrzi_yaxis.size() << endl;
		cout << "bami size: " << bami_yaxis.size() << endl;

	// sygnal Manchester
	// ---------------------
		vector<float> manch_xaxis, manch_yaxis;

		generateManchesterSignal(manch_xaxis, manch_yaxis, clk_xaxis, clk_yaxis, ttl_xaxis, ttl_yaxis);

		// wykres sygnalu Manchester
		title << std::fixed << std::setprecision(2)
			<< "sygnal Manchester"
			<< ", T_b = " << Tb << " s"
			<< ", f = " << f << " Hz"
			<< ", fs = " << fs << " Hz";

		filename = namePrefix + "manch";
		drawSignalChart(manch_xaxis, manch_yaxis, 16 * num_samples, "lines", title.str(), filename, "T_b[s]", "A", 1920, 800, 0.1, 0, ymax * 0.1);
		title.str("");

	// dekoder TTL
	// ---------------
		vector<float> d_ttl_xaxis, d_ttl_yaxis;

		decodeTTLSignal(d_ttl_xaxis, d_ttl_yaxis, ttl_xaxis, ttl_yaxis);

		// wykres sygnalu TTL
		title << std::fixed << std::setprecision(2)
			<< "zdekodowany sygnal TTL"
			<< ", T_b = " << Tb << " s"
			<< ", f = " << f << " Hz"
			<< ", fs = " << fs << " Hz";

		filename = namePrefix + "ttl_decoded";
		drawSignalChart(d_ttl_xaxis, d_ttl_yaxis, 16 * num_samples, "lines", title.str(), filename, "T_b[s]", "A", 1920, 800, 0.1, 0, ymax * 0.1);
		title.str("");

	// dekoder NRZI
	// ----------------
		vector<float> d_nrzi_xaxis, d_nrzi_yaxis;

		decodeNRZISignal(d_nrzi_xaxis, d_nrzi_yaxis, nrzi_xaxis, nrzi_yaxis, num_samples);

		// wykres sygnalu NRZI
		title << std::fixed << std::setprecision(2)
			<< "zdekodowany sygnal NRZI"
			<< ", T_b = " << Tb << " s"
			<< ", f = " << f << " Hz"
			<< ", fs = " << fs << " Hz";

		filename = namePrefix + "nrzi_decoded";
		drawSignalChart(d_nrzi_xaxis, d_nrzi_yaxis, 16 * num_samples, "lines", title.str(), filename, "T_b[s]", "A", 1920, 800, 0.1, 0, ymax * 0.1);
		title.str("");

	// dekoder BAMI
	// ----------------
		vector<float> d_bami_xaxis, d_bami_yaxis;

		decodeBAMISignal(d_bami_xaxis, d_bami_yaxis, bami_xaxis, bami_yaxis);

		// wykres sygnalu BAMI
		title << std::fixed << std::setprecision(2)
			<< "zdekodowany sygnal BAMI"
			<< ", T_b = " << Tb << " s"
			<< ", f = " << f << " Hz"
			<< ", fs = " << fs << " Hz";

		filename = namePrefix + "bami_decoded";
		drawSignalChart(d_bami_xaxis, d_bami_yaxis, 16 * num_samples, "lines", title.str(), filename, "T_b[s]", "A", 1920, 800, 0.1, 0, ymax * 0.1);
		title.str("");

	// dekoder Manchester
	// ----------------------
		vector<float> d_manch_xaxis, d_manch_yaxis;

		decodeManchesterSignal(d_manch_xaxis, d_manch_yaxis, manch_xaxis, manch_yaxis, f, Tb, tN, num_samples);

		// wykres sygnalu Manchester
		title << std::fixed << std::setprecision(2)
			<< "zdekodowany sygnal Manchester"
			<< ", T_b = " << Tb << " s"
			<< ", f = " << f << " Hz"
			<< ", fs = " << fs << " Hz";

		filename = namePrefix + "manch_decoded";
		drawSignalChart(d_manch_xaxis, d_manch_yaxis, 16 * num_samples, "lines", title.str(), filename, "T_b[s]", "A", 1920, 800, 0.1, 0, ymax * 0.1);
		title.str("");

	getchar();
	return 0;
}