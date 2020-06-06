#ifndef _CODERS_MODULATORS
#define _CODERS_MODULATORS

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
vector<bitset<8>> secded(vector<bitset<8>> binStream);
vector<bitset<4>> decodeSECDED(vector<bitset<8>> scd);

// modulacje/demodulacje
void generateClockSignal(vector<float> &t, vector<float> &x, float &f, float &fi, float &tN);
bool isClkDown(float &clk_curr, float &clk_next);
bool isClkUp(float &clk_curr, float &clk_next);
void generateTTLSignal(vector<float> &ttl_xaxis, vector<float> &ttl_yaxis, vector<bitset<8>> &s, float &f);
void desampleSignal(vector<float> signal_xaxis, vector<float> signal_yaxis, vector<float> &desampled_xaxis, vector<float> &desampled_yaxis, int num_samples, int start_pos = 0);
void generateNRZISignal(vector<float> &nrzi_xaxis, vector<float> &nrzi_yaxis, vector<float> &clk_xaxis, vector<float> &clk_yaxis, vector<float> &ttl_xaxis, vector<float> &ttl_yaxis);
void generateBAMISignal(vector<float> &bami_xaxis, vector<float> &bami_yaxis, vector<float> &clk_xaxis, vector<float> &clk_yaxis, vector<float> &ttl_xaxis, vector<float> &ttl_yaxis);
void generateManchesterSignal(vector<float> &manch_xaxis, vector<float> &manch_yaxis, vector<float> &clk_xaxis, vector<float> &clk_yaxis, vector<float> &ttl_xaxis, vector<float> &ttl_yaxis);

void decodeTTLSignal(vector<float> &d_ttl_xaxis, vector<float> &d_ttl_yaxis, vector<float> ttl_xaxis, vector<float> ttl_yaxis);
void decodeNRZISignal(vector<float> &d_nrzi_xaxis, vector<float> &d_nrzi_yaxis, vector<float> nrzi_xaxis, vector<float> nrzi_yaxis, float &f);
void decodeBAMISignal(vector<float> &d_bami_xaxis, vector<float> &d_bami_yaxis, vector<float> bami_xaxis, vector<float> bami_yaxis);
void decodeManchesterSignal(vector<float> &d_manch_xaxis, vector<float> &d_manch_yaxis, vector<float> manch_xaxis, vector<float> manch_yaxis, float &f, float &tN);

#endif