#ifndef _HELPER_FUNCTIONS
#define _HELPER_FUNCTIONS

#include <iostream>
#include <string>
#include <vector>
#include <bitset>
#include "structs.h"

using std::cout;
using std::endl;
using std::string;
using std::bitset;
using std::vector;

void dataToFile(string outputFile, string data);
bitset<8> reverseBitset(bitset<8> bits);
bitset<7> reverseBitset(bitset<7> bits);
bitset<4> reverseBitset(bitset<4> bits);
vector<bitset<8>> strToBinStream(string input, string option);
string binStreamToStr(vector<bitset<8>> input, string option);
vector<bitset<8>> signalToBitsets(vector<float> signal);
string bitsetVectorToStr(vector<bitset<8>> input, bool splitHalf = false);

namespace chart_utils {
	void dataToCsv(string, string);
	void drawChartToPDF(string, string, string, string, string, int, int, float, float, float);
	void drawChartsToPDF(string, vector<signal>, float, float, float);
	void drawChartsMultiplotToPDF(string, string, vector<signal>, vector<int>, float, float, float);
	void drawChart(string, string, string, string, string, int, int, float, float, float);
}

void drawSignalChart(vector<float> &, vector<float> &, int, string, string, string, string, string, int, int, float xTics = 0, float xOffset = 0, float yOffset = 0);

void drawSignalsCharts(vector<signal>, int, string, float xTics = 0, float xOffset = 0, float yOffset = 0);

void drawSignalsMultiCharts(vector<signal>, vector<int>, int, string, string, float xTics = 0, float xOffset = 0, float yOffset = 0);

#endif