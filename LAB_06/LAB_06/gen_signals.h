#ifndef _GEN_SIGNALS
#define _GEN_SIGNALS

#include <vector>
#include <bitset>

using std::vector;
using std::bitset;

void generateInfSignal(vector<float> &, vector<float> &, vector<bitset<8>> &, float &, float);

void generateAmplModSignal(vector<float> &, vector<float> &, vector<float> &, float &, float &, float &, float, float, float);

void generateFreqModSignal(vector<float> &, vector<float> &, vector<float> &, float &, float &, float &, float &, float, float);

void generatePhaseModSignal(vector<float> &, vector<float> &, vector<float> &, float &, float &, float &, float &, float, float);

#endif