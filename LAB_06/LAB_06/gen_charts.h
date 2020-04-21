#ifndef _GEN_CHARTS
#define _GEN_CHARTS

#include <vector>
#include <string>

using std::vector;
using std::string;

void drawSignalChart(vector<float> &, vector<float> &, int, string, string, string, string, string, int, int);

void drawSignalChart(vector<float> &, vector<float> &, int, string, string, string, string, string, int, int, float, float);

void drawSignalMultiChart(vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, int, string, string, string, string, string, string, string, string, string, string, string, string, string, int, int);

void drawSpectrumChart(vector<float> &, vector<float> &, int, string, string, string, string, int, int);

void drawSpectrumPartChart(vector<float> &, vector<float> &, int, float, string, string, string, string, int, int);

void drawSpectrumPartChart(vector<float> &, vector<float> &, int, float, string, string, string, string, int, int, int);

void drawSpectrumPartdBChart(vector<float> &, vector<float> &, int, string, string, string, string, int, int);

#endif