#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include "structs.h"
#include "gnuplot_i.hpp"

using std::string;
using std::to_string;
using std::stringstream;
using std::ofstream;
using std::ios;
using std::vector;

namespace chart_utils {
	void dataToCsv(string outputFile, string data) {
		ofstream of { outputFile };

		if (of.is_open())
		{
			of << data;
			of.close();
		}
		else {
			std::cout << "Nie udalo sie otworzyc: " << outputFile << std::endl;
		}
	}

	void drawChartToPDF(string title, string xlabel, string ylabel, string inputFile, string outputFile, int width, int height, float xTics, float xOffset, float yOffset) {
		Gnuplot gplot;

		gplot.cmd("set terminal pdf size 29.7 cm, 21 cm enhanced font 'Verdana,8'");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("set encoding utf8");
		gplot.cmd("set size ratio 0.4");
		gplot.cmd("set title \"" + title + "\"");
		gplot.cmd("set xlabel \"" + xlabel + "\"");
		gplot.cmd("set ylabel \"" + ylabel + "\"");

		if (xTics) {
			gplot.cmd("set xtics " + to_string(xTics));
			gplot.cmd("set mxtics 2");
		}

		gplot.cmd("set offsets " + to_string(xOffset) + ", " + to_string(xOffset) + ", " + to_string(yOffset) + ", " + to_string(yOffset));
		gplot.cmd("set style line 100 lt 2 dt 3");
		gplot.cmd("set style line 101 lt 2 dt 2");
		gplot.cmd("set grid xtics mxtics ls 101, ls 100");
		gplot.cmd("set xzeroaxis ls 101");
		gplot.cmd("set autoscale xfixmin");
		gplot.cmd("set autoscale xfixmax");
		gplot.cmd("set style line 1 linetype 1 linewidth 1");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("plot \"./csv/" + inputFile + ".csv\" u 1:2 with lines linestyle 1 lc '#4169E1");
	}

	void drawChartsToPDF(string outputFile, vector<signal> signals, float xTics, float xOffset, float yOffset) {
		Gnuplot gplot;
		int count = signals.size();

		gplot.cmd("set terminal pdf size 29.7 cm, 21 cm enhanced font 'Verdana,8'");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("set encoding utf8");
		gplot.cmd("set size ratio 0.4");
		gplot.cmd("set offsets " + to_string(xOffset) + ", " + to_string(xOffset) + ", " + to_string(yOffset) + ", " + to_string(yOffset));
		gplot.cmd("set style line 100 lt 2 dt 3");
		gplot.cmd("set style line 101 lt 2 dt 2");
		gplot.cmd("set grid xtics mxtics ls 101, ls 100");
		gplot.cmd("set xzeroaxis ls 101");
		gplot.cmd("set autoscale xfixmin");
		gplot.cmd("set autoscale xfixmax");
		gplot.cmd("set style line 1 linetype 1 linewidth 1");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");

		if (xTics) {
			gplot.cmd("set xtics " + to_string(xTics));
			gplot.cmd("set mxtics 2");
		}

		for (int i = 0; i < count; i++) {
			gplot.cmd("set title \"" + signals[i].title + "\"");
			gplot.cmd("set xlabel \"" + signals[i].xlabel + "\"");
			gplot.cmd("set ylabel \"" + signals[i].ylabel + "\"");
			gplot.cmd("plot \"./csv/" + signals[i].filename + ".csv\" u 1:2 with lines linestyle 1 lc '#4169E1");
		}
	}

	void drawChartsMultiplotToPDF(string title, string outputFile, vector<signal> signals, vector<int> plotsPerMPlot, float xTics, float xOffset, float yOffset) {
		Gnuplot gplot;

		gplot.cmd("set terminal pdf size 21 cm, 29.7 cm enhanced font 'Verdana,6'");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("set encoding utf8");
		gplot.cmd("set offsets " + to_string(xOffset) + ", " + to_string(xOffset) + ", " + to_string(yOffset) + ", " + to_string(yOffset));
		gplot.cmd("set style line 100 lt 2 dt 3");
		gplot.cmd("set style line 101 lt 2 dt 2");
		gplot.cmd("set grid xtics mxtics ls 101, ls 100");
		gplot.cmd("set xzeroaxis ls 101");
		gplot.cmd("set autoscale xfixmin");
		gplot.cmd("set autoscale xfixmax");
		gplot.cmd("set style line 1 linetype 1 linewidth 1");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");

		if (xTics) {
			gplot.cmd("set xtics " + to_string(xTics));
			gplot.cmd("set mxtics 2");
		}

		int plotCount = 0;
		int tmp = 0;

		for (int i = 0; i < plotsPerMPlot.size(); i++) {
			gplot.cmd("set multiplot layout " + to_string(plotsPerMPlot[i]) + ", 1 title \"" + title + "\" font \", 14\"");

			tmp = plotCount;

			for (int j = tmp; j < tmp + plotsPerMPlot[i]; j++) {
				gplot.cmd("set title \"" + signals[j].title + "\"");
				gplot.cmd("set xlabel \"" + signals[j].xlabel + "\"");
				gplot.cmd("set ylabel \"" + signals[j].ylabel + "\"");
				gplot.cmd("plot \"./csv/" + signals[j].filename + ".csv\" u 1:2 with lines linestyle 1 lc '#4169E1");
				plotCount++;
			}
		}
	}

	void drawChart(string title, string xlabel, string ylabel, string inputFile, string outputFile, int width, int height, float xTics, float xOffset, float yOffset) {
		Gnuplot gplot;

		gplot.cmd("set terminal pngcairo size " + to_string(width) + ", " + to_string(height) + " enhanced font 'Verdana,8'");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("set encoding utf8");
		gplot.cmd("set title \"" + title + "\"");
		gplot.cmd("set xlabel \"" + xlabel + "\"");
		gplot.cmd("set ylabel \"" + ylabel + "\"");

		if (xTics) {
			gplot.cmd("set xtics " + to_string(xTics));
			gplot.cmd("set mxtics 2");
		}

		gplot.cmd("set offsets " + to_string(xOffset) + ", " + to_string(xOffset) + ", " + to_string(yOffset) + ", " + to_string(yOffset));
		gplot.cmd("set style line 100 lt 2 dt 3");
		gplot.cmd("set style line 101 lt 2 dt 2");
		gplot.cmd("set grid xtics mxtics ls 101, ls 100");
		gplot.cmd("set xzeroaxis ls 101");
		gplot.cmd("set autoscale xfixmin");
		gplot.cmd("set autoscale xfixmax");
		gplot.cmd("set style line 1 linetype 1 linewidth 1");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("plot \"./" + inputFile + "\" u 1:2 with lines linestyle 1 lc '#4169E1");
	}
}
