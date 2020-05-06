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
#include "structs.h"
#include "chart_utils.h"
#include "gen_charts.h"


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
float f = 10.0;

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
void generateClockSignal(vector<float> &t, vector<float> &x, float &f, float &fi, float &tN) {
	float a = 1.0;
	int n = 0;
	const int num_samples = 100 * f;
	float Tb = 1.0 / f, dt = Tb / num_samples;
	float t0 = 0.0, tn = t0 + (n*dt);
	
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

bool isClkDown(float &clk_curr, float &clk_next) {
	return ((clk_curr != clk_next) && (clk_next == 0));
}

bool isClkUp(float &clk_curr, float &clk_next) {
	return ((clk_curr != clk_next) && (clk_next == 1));
}

void generateTTLSignal(vector<float> &ttl_xaxis, vector<float> &ttl_yaxis, vector<bitset<8>> &s, float &f) {
	int n = 0;
	const int num_samples = 100 * f;
	float Tb = 1.0 / f, dt = Tb / num_samples;
	float t0 = 0.0, tn = t0 + (n*dt);

	for (int i = 0; i < s.size(); i++) {
		for (int j = 7; j >= 0; j--) {
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
	int state = 0;

	int lastIndex = 0;

	for (int i = 0; i < ttl_yaxis.size(); i++) {
		if (isClkDown(clk_yaxis[i], clk_yaxis[i + 1])) {
			if (ttl_yaxis[i] == 1)
				state = -1;
			else
				state = 1;
		}

		if (((i - 1 > 0) && (i + 1 < ttl_yaxis.size())) && (isClkUp(clk_yaxis[i], clk_yaxis[i + 1]))) {
			if (ttl_yaxis[i - 1] == ttl_yaxis[i + 1])
				state = -state;
		}

		manch_yaxis.push_back(state);
		manch_xaxis.push_back(clk_xaxis[i]);

		lastIndex = i;
	}
}

// dekodery
void decodeTTLSignal(vector<float> &d_ttl_xaxis, vector<float> &d_ttl_yaxis, vector<float> ttl_xaxis, vector<float> ttl_yaxis) {
	for (int i = 0; i < ttl_yaxis.size(); i++) {
		d_ttl_yaxis.push_back(ttl_yaxis[i]);
		d_ttl_xaxis.push_back(ttl_xaxis[i]);
	}
}

void decodeNRZISignal(vector<float> &d_nrzi_xaxis, vector<float> &d_nrzi_yaxis, vector<float> nrzi_xaxis, vector<float> nrzi_yaxis, float &f) {
	const int num_samples = 100 * f;
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

void decodeBAMISignal(vector<float> &d_bami_xaxis, vector<float> &d_bami_yaxis, vector<float> bami_xaxis, vector<float> bami_yaxis) {
	for (int i = 0; i < bami_yaxis.size(); i++) {
		d_bami_yaxis.push_back(abs(bami_yaxis[i]));
		d_bami_xaxis.push_back(bami_xaxis[i]);
	}
}

void decodeManchesterSignal(vector<float> &d_manch_xaxis, vector<float> &d_manch_yaxis, vector<float> manch_xaxis, vector<float> manch_yaxis, float &f, float &tN) {
	vector<float> s_clk_xaxis, s_clk_yaxis;
	const int num_samples = 100 * f;
	float Tb = 1.0 / f;
	float fi = -M_PI / 2;
	int state = 0;

	generateClockSignal(s_clk_xaxis, s_clk_yaxis, f, fi, tN);

	for (int i = 0; i < manch_yaxis.size(); i++) {
		if (isClkDown(s_clk_yaxis[i], s_clk_yaxis[i + 1])) {
			state = manch_yaxis[i];

			if (state == -1)
				state = 0;

			state = !state;
		}

		d_manch_yaxis.push_back(state);
		d_manch_xaxis.push_back(s_clk_xaxis[i]);
	}
}

int main() {
	string namePrefix = "LAB_07_", filename = "";
	stringstream title;

	vector<bitset<8>> s1 = strToBinStream("KOT", "bigEndian");

	int samplesToPlot = 16 * (100 * f);

	float Tb = 1.0 / f, tN = Tb * (s1.size() * 8);
	float fi = 0;

	// sygnal zegarowy (CLK)
	// -------------------------
		vector<float> clk_xaxis, clk_yaxis;
		generateClockSignal(clk_xaxis, clk_yaxis, f, fi, tN);

		title << std::fixed << std::setprecision(2)
			<< "sygnal zegarowy (CLK)"
			<< ", f = " << f << " Hz";

		filename = namePrefix + "clk";
		signal clk(clk_xaxis, clk_yaxis, title.str(), filename, "t[s]", "A");
		title.str("");

	// sygnal TTL
	// -----------------------------
		vector<float> ttl_xaxis, ttl_yaxis;
		generateTTLSignal(ttl_xaxis, ttl_yaxis, s1, f);

		float ymax = *max_element(ttl_yaxis.begin(), ttl_yaxis.end());

		// wykres sygnalu TTL
		title << std::fixed << std::setprecision(2)
			<< "sygnal TTL"
			<< ", f = " << f << " Hz"
			<< ", T_b = " << Tb << " s";

		filename = namePrefix + "ttl";
		signal ttl(ttl_xaxis, ttl_yaxis, title.str(), filename, "t[s]", "A");
		title.str("");

	// sygnal NRZI
	// ---------------
		vector<float> nrzi_xaxis, nrzi_yaxis;

		generateNRZISignal(nrzi_xaxis, nrzi_yaxis, clk_xaxis, clk_yaxis, ttl_xaxis, ttl_yaxis);

		// wykres sygnalu NRZI
		title << std::fixed << std::setprecision(2)
			<< "sygnal NRZI"
			<< ", f = " << f << " Hz"
			<< ", T_b = " << Tb << " s";

		filename = namePrefix + "nrzi";
		signal nrzi(nrzi_xaxis, nrzi_yaxis, title.str(), filename, "t[s]", "A");
		title.str("");

	// sygnal BAMI
	// ---------------
		vector<float> bami_xaxis, bami_yaxis;

		generateBAMISignal(bami_xaxis, bami_yaxis, clk_xaxis, clk_yaxis, ttl_xaxis, ttl_yaxis);

		// wykres sygnalu BAMI
		title << std::fixed << std::setprecision(2)
			<< "sygnal BAMI"
			<< ", f = " << f << " Hz"
			<< ", T_b = " << Tb << " s";

		filename = namePrefix + "bami";
		signal bami(bami_xaxis, bami_yaxis, title.str(), filename, "t[s]", "A");

		title.str("");

	// sygnal Manchester (G.E. Thomas)
	// -----------------------------------
		vector<float> manch_xaxis, manch_yaxis;

		generateManchesterSignal(manch_xaxis, manch_yaxis, clk_xaxis, clk_yaxis, ttl_xaxis, ttl_yaxis);

		// wykres sygnalu Manchester
		title << std::fixed << std::setprecision(2)
			<< "sygnal Manchester (G.E. Thomas)"
			<< ", f = " << f << " Hz"
			<< ", T_b = " << Tb << " s";

		filename = namePrefix + "manch";
		signal manch(manch_xaxis, manch_yaxis, title.str(), filename, "t[s]", "A");

		title.str("");

	// dekoder TTL
	// ---------------
		vector<float> d_ttl_xaxis, d_ttl_yaxis;

		decodeTTLSignal(d_ttl_xaxis, d_ttl_yaxis, ttl_xaxis, ttl_yaxis);

		// wykres sygnalu TTL
		title << std::fixed << std::setprecision(2)
			<< "zdekodowany sygnal TTL"
			<< ", f = " << f << " Hz"
			<< ", T_b = " << Tb << " s";

		filename = namePrefix + "ttl_decoded";
		signal d_ttl(d_ttl_xaxis, d_ttl_yaxis, title.str(), filename, "t[s]", "A");
		title.str("");

	// dekoder NRZI
	// ----------------
		vector<float> d_nrzi_xaxis, d_nrzi_yaxis;

		decodeNRZISignal(d_nrzi_xaxis, d_nrzi_yaxis, nrzi_xaxis, nrzi_yaxis, f);

		// wykres sygnalu NRZI
		title << std::fixed << std::setprecision(2)
			<< "zdekodowany sygnal NRZI"
			<< ", f = " << f << " Hz"
			<< ", T_b = " << Tb << " s";

		filename = namePrefix + "nrzi_decoded";
		signal d_nrzi(d_nrzi_xaxis, d_nrzi_yaxis, title.str(), filename, "t[s]", "A");
		title.str("");

	// dekoder BAMI
	// ----------------
		vector<float> d_bami_xaxis, d_bami_yaxis;

		decodeBAMISignal(d_bami_xaxis, d_bami_yaxis, bami_xaxis, bami_yaxis);

		// wykres sygnalu BAMI
		title << std::fixed << std::setprecision(2)
			<< "zdekodowany sygnal BAMI"
			<< ", f = " << f << " Hz"
			<< ", T_b = " << Tb << " s";

		filename = namePrefix + "bami_decoded";
		signal d_bami(d_bami_xaxis, d_bami_yaxis, title.str(), filename, "t[s]", "A");
		title.str("");

	// dekoder Manchester (G.E. Thomas)
	// ------------------------------------
		vector<float> d_manch_xaxis, d_manch_yaxis;

		decodeManchesterSignal(d_manch_xaxis, d_manch_yaxis, manch_xaxis, manch_yaxis, f, tN);

		// wykres sygnalu Manchester
		title << std::fixed << std::setprecision(2)
			<< "zdekodowany sygnal Manchester (G.E. Thomas)"
			<< ", f = " << f << " Hz"
			<< ", T_b = " << Tb << " s";

		filename = namePrefix + "manch_decoded";
		signal d_manch(d_manch_xaxis, d_manch_yaxis, title.str(), filename, "t[s]", "A");
		title.str("");

	// PDF z pojedynczymi wykresami sygnalow
	// -----------------------------------------
		vector<signal> signalsToPlot;

		signalsToPlot.push_back(clk);
		signalsToPlot.push_back(ttl);
		signalsToPlot.push_back(nrzi);
		signalsToPlot.push_back(manch);
		signalsToPlot.push_back(bami);

		signalsToPlot.push_back(ttl);
		signalsToPlot.push_back(d_ttl);
		signalsToPlot.push_back(d_nrzi);
		signalsToPlot.push_back(d_manch);
		signalsToPlot.push_back(d_bami);

		cout << endl << "Generowanie pliku PDF z pojedynczymi wykresami..." << endl;

		drawSignalsCharts(signalsToPlot, samplesToPlot, namePrefix + "wykresy_kod_dekod", 0.1, 0, ymax * 0.1);

	// PDF ze zbiorczymi wykresami sygnalow
	// ----------------------------------------
		vector<int> plotsPerMPlot;
		plotsPerMPlot.push_back(5);
		plotsPerMPlot.push_back(5);

		cout << "Generowanie pliku PDF ze zbiorczymi wykresami..." << endl;

		drawSignalsMultiCharts(signalsToPlot, plotsPerMPlot, samplesToPlot, "", namePrefix + "wykresy_zbiorcze", 0.1, 0, ymax * 0.1);

		cout << endl << "Nacisnij dowolny klawisz, aby zakonczyc..." << endl;

	getchar();
	return 0;
}
