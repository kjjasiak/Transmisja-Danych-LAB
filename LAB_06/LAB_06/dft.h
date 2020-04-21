#ifndef _DFT
#define _DFT

#include <vector>
#include <complex>

using std::vector;
using std::complex;

vector<float> getAmplSpectrum(vector<complex<float>>, int);
vector<float> getAmplSpectrumBis(vector<float>, float, int);
vector<float> getFreqScale(int, int);
vector<complex<float>> dft(vector<float>, int);
vector<float> idft(vector<complex<float>>, int);
float getBandWidth(vector<float>, vector<float>, int, float);

#endif