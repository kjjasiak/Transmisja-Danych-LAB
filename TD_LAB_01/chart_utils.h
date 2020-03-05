#ifndef _CHART_UTILS
#define _CHART_UTILS

#include <string>
#include "./gnuplot-cpp/gnuplot_i.hpp"

using std::string;

void drawChart(string inputFile, string outputFile) {
	Gnuplot gplot;

	gplot.cmd("set terminal png");
	gplot.cmd("set output \"./" + outputFile + "\"");
	gplot.cmd("set datafile separator ';'");
	gplot.cmd("plot \"./" + inputFile + "\"");
}

#endif