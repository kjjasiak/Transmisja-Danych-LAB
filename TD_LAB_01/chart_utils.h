#ifndef _CHART_UTILS
#define _CHART_UTILS

#include <string>

using std::string;

namespace chart_utils {
	string generateData(float, float, float, float, float(*fn)(float));
	void dataToCsv(string, string);
	void drawChart(string, string);
}

#endif