#ifndef _CHART_UTILS
#define _CHART_UTILS

#include <string>
#include <vector>
#include "structs.h"

using std::string;
using std::vector;


namespace chart_utils {
	void dataToCsv(string, string);
	void drawChartToPDF(string, string, string, string, string, int, int, float, float, float);
	void drawChartsToPDF(string, vector<signal>, float, float, float);
	void drawChartsMultiplotToPDF(string, string, vector<signal>, vector<int>, float, float, float);
	void drawChart(string, string, string, string, string, int, int, float, float, float);
}

#endif