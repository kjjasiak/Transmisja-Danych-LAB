#ifndef _GEN_CHARTS
#define _GEN_CHARTS

#include <vector>
#include <string>

using std::vector;
using std::string;

void drawSignalChart(vector<float> &, vector<float> &, int, string, string, string, string, string, int, int, float xTics = 0, float xOffset = 0, float yOffset = 0);

void drawSignalMultiChart(vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, int, string, string, string, string, string, string, string, string, string, string, string, string, string, string, string, int, int);

void drawSignalMultiChart(vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, vector<float> &, int, string, string, string, string, string, string, string, string, string, string, string, string, string, string, string, string, string, int, int);

void drawSpectrumChart(vector<float> &, vector<float> &, int, string, string, string, string, int, int);

void drawSpectrumPartChart(vector<float> &, vector<float> &, int, float, string, string, string, string, int, int);

void drawSpectrumPartChart(vector<float> &, vector<float> &, int, float, string, string, string, string, int, int, int);

void drawSpectrumPartdBChart(vector<float> &, vector<float> &, int, string, string, string, string, int, int);

#endif