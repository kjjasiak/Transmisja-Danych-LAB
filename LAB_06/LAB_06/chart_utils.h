#ifndef _CHART_UTILS
#define _CHART_UTILS

#include <string>

using std::string;

namespace chart_utils {
	void dataToCsv(string, string);
	void drawChart(string, string, string, string, string, int, int);
	void drawChart(string, string, string, string, string, int, int, float, float);
	void drawChartWithDots(string, string, string, string, string, int, int);
	void drawChartWithSteps(string, string, string, string, string, int, int);
	void drawChartWithImpulses(string, string, string, string, string, int, int);
	void drawChartMultiPlot(string, string, string, string, string, string, string, string, string, string, string, string, string, string, int, int, float, float);
	void drawChartMultiPlot(string, string, string, string, string, string, string, string, string, string, string, string, string, string, string, string, int, int, float, float);
}

#endif