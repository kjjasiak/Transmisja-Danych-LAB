#ifndef _GEN_CHARTS
#define _GEN_CHARTS

#include <vector>
#include <string>

using std::vector;
using std::string;

void drawSignalChart(vector<float> &, vector<float> &, int, string, string, string, string, string, int, int, float xTics = 0, float xOffset = 0, float yOffset = 0);

void drawSignalsCharts(vector<signal>, int, string, float xTics = 0, float xOffset = 0, float yOffset = 0);

void drawSignalsMultiCharts(vector<signal>, vector<int>, int, string, string, float xTics = 0, float xOffset = 0, float yOffset = 0);

#endif