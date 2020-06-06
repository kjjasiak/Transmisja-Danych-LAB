#define _USE_MATH_DEFINES

#include <iostream>
#include <bitset>
#include <vector>
#include <cmath>
#include <sstream>
#include "helper_functions.h"

using std::cout;
using std::endl;
using std::bitset;
using std::vector;
using std::stringstream;


// kodowania/dekodowania
vector<bitset<8>> secded(vector<bitset<8>> binStream) {
	vector<bitset<4>> dVec;

	for (int i = 0; i < binStream.size(); i++) {
		bitset<8> streamByte = binStream[i];
		stringstream st;

		streamByte = reverseBitset(streamByte);

		for (int j = 0; j < 4; j++)
			st << streamByte[j];

		bitset<4> bits(st.str());
		dVec.push_back(bits);
		st.str("");

		for (int j = 4; j < 8; j++)
			st << streamByte[j];

		bitset<4> bits2(st.str());
		dVec.push_back(bits2);
		st.str("");
	}

	vector<bitset<8>> hBitsVec;

	for (int i = 0; i < dVec.size(); i++) {
		bitset<4> d = dVec[i];

		d = reverseBitset(d);

		bool p1 = (d[1 - 1] + d[2 - 1] + d[4 - 1]) % 2,
			p2 = (d[1 - 1] + d[3 - 1] + d[4 - 1]) % 2,
			p3 = (d[2 - 1] + d[3 - 1] + d[4 - 1]) % 2,
			p4 = (p1 + p2 + d[1 - 1] + p3 + d[2 - 1] + d[3 - 1] + d[4 - 1]) % 2;

		stringstream st;
		st << p1 << p2 << d[1 - 1] << p3 << d[2 - 1] << d[3 - 1] << d[4 - 1] << p4;

		bitset<8> h(st.str());
		hBitsVec.push_back(h);

		st.str("");
	}

	return hBitsVec;
}

vector<bitset<4>> decodeSECDED(vector<bitset<8>> scd) {
	vector<bitset<4>> dVec;
	stringstream st, out;

	for (int i = 0; i < scd.size(); i++) {
		bitset<8> h = scd[i];
		h = reverseBitset(h);

		bool p4 = (h[1 - 1] + h[2 - 1] + h[3 - 1] + h[4 - 1] + h[5 - 1] + h[6 - 1] + h[7 - 1]) % 2;

		if (p4 != h[8 - 1]) {
			out << "#" << i + 1 << ":" << endl;
			out << "Wystapil blad, nastapi proba naprawy..." << endl;
		}

		bool p1 = (h[1 - 1] + h[3 - 1] + h[5 - 1] + h[7 - 1]) % 2,
			p2 = (h[2 - 1] + h[3 - 1] + h[6 - 1] + h[7 - 1]) % 2,
			p3 = (h[4 - 1] + h[5 - 1] + h[6 - 1] + h[7 - 1]) % 2;

		int n = p1 * pow(2, 0) + p2 * pow(2, 1) + p3 * pow(2, 2);

		if (n > 0) {
			h[n - 1] = !h[n - 1];

			out << "#" << i + 1 << ":" << endl;
			out << "Korekta na bicie nr " << n;

			if (p4 == h[8 - 1])
				out << " (XOR numerow bitow z bledami)";

			out << endl;

			p4 = (h[1 - 1] + h[2 - 1] + h[3 - 1] + h[4 - 1] + h[5 - 1] + h[6 - 1] + h[7 - 1]) % 2;

			if (p4 == h[8 - 1]) {
				out << "Wystapil 1 blad w transmisji, ktory naprawiono" << endl << endl;
			}
			else {
				out << "Wykryto co najmniej 2 bledy w transmisji" << endl;
				out << "Pakiet odrzucono, zdekodowany strumien jest niepoprawny - konieczna ponowna transmisja" << endl << endl;
			}
		}

		st << h[3 - 1] << h[5 - 1] << h[6 - 1] << h[7 - 1];

		bitset<4> d(st.str());
		dVec.push_back(d);
		st.str("");
	}

	return dVec;
}


// modulacje/demodulacje
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

void desampleSignal(vector<float> signal_xaxis, vector<float> signal_yaxis, vector<float> &desampled_xaxis, vector<float> &desampled_yaxis, int num_samples, int start_pos = 0) {
	for (int i = start_pos; i < signal_yaxis.size(); i = i + (num_samples)) {
		desampled_xaxis.push_back(signal_xaxis[i]);
		desampled_yaxis.push_back(signal_yaxis[i]);
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
		if ((i > 0) && isClkUp(clk_yaxis[i], clk_yaxis[i + 1])) {
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