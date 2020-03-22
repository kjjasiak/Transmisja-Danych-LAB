#ifndef _CHART_UTILS
#define _CHART_UTILS

#include <string>

using std::string;

namespace chart_utils {
	string generateData(float, float, float, float, float(*fn)(float));
	string generateData(float, float, float, float, int, float(*fn)(float, int));
	void dataToCsv(string, string);
	void drawChart(string, string, string, string, string);
	void drawChartWithDots(string, string, string, string, string);
	void drawChartWithSteps(string, string, string, string, string);
}

#endif