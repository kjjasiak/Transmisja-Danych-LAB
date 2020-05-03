#ifndef _GEN_SIGNALS
#define _GEN_SIGNALS

#include <vector>
#include <bitset>

using std::vector;
using std::bitset;
using std::complex;


void generateInfSignal(vector<float> &, vector<float> &, vector<bitset<8>> &, float &, float);

void generateAmplModSignal(vector<float> &, vector<float> &, vector<float> &, float &, float &, float &, float, float, float);

void generateFreqModSignal(vector<float> &, vector<float> &, vector<float> &, float &, float &, float &, float &, float, float);

void generatePhaseModSignal(vector<float> &, vector<float> &, vector<float> &, float &, float &, float &, float &, float, float);

vector<float> getAmplSpectrum(vector<complex<float>>, int);
vector<float> getAmplSpectrumBis(vector<float>, float, int);
vector<float> getFreqScale(int, int);
vector<complex<float>> dft(vector<float>, int);
vector<float> idft(vector<complex<float>>, int);
float getBandWidth(vector<float>, vector<float>, int, float);

#endif