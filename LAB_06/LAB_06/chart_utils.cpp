#include <string>
#include <sstream>
#include <fstream>
#include "gnuplot_i.hpp"

using std::string;
using std::to_string;
using std::stringstream;
using std::ofstream;
using std::ios;

namespace chart_utils {
	void dataToCsv(string outputFile, string data) {
		ofstream of;

		of.open(outputFile, ios::out | ios::trunc);
		of << data;
		of.close();
	}

	void drawChart(string title, string xlabel, string ylabel, string inputFile, string outputFile, int width, int height) {
		Gnuplot gplot;

		gplot.cmd("set terminal pngcairo size " + to_string(width) + ", " + to_string(height) + " enhanced font 'Verdana,8'");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("set encoding utf8");
		gplot.cmd("set title \"" + title + "\"");
		gplot.cmd("set xlabel \"" + xlabel + "\"");
		gplot.cmd("set ylabel \"" + ylabel + "\"");
		gplot.cmd("set mxtics");
		gplot.cmd("set autoscale xfixmin");
		gplot.cmd("set autoscale xfixmax");
		gplot.cmd("set style line 1 linetype 1 linewidth 1");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("plot \"./" + inputFile + "\" u 1:2 with lines linestyle 1 lc '#4169E1");
	}

	void drawChart(string title, string xlabel, string ylabel, string inputFile, string outputFile, int width, int height, float xOffset, float yOffset) {
		Gnuplot gplot;

		gplot.cmd("set terminal pngcairo size " + to_string(width) + ", " + to_string(height) + " enhanced font 'Verdana,8'");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("set encoding utf8");
		gplot.cmd("set title \"" + title + "\"");
		gplot.cmd("set xlabel \"" + xlabel + "\"");
		gplot.cmd("set ylabel \"" + ylabel + "\"");
		gplot.cmd("set mxtics");
		gplot.cmd("set offsets " + to_string(xOffset) + ", " + to_string(xOffset) + ", " + to_string(yOffset) + ", " + to_string(yOffset));
		gplot.cmd("set autoscale xfixmin");
		gplot.cmd("set autoscale xfixmax");
		gplot.cmd("set style line 1 linetype 1 linewidth 1");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("plot \"./" + inputFile + "\" u 1:2 with lines linestyle 1 lc '#4169E1");
	}

	void drawChartWithDots(string title, string xlabel, string ylabel, string inputFile, string outputFile, int width, int height) {
		Gnuplot gplot;

		gplot.cmd("set terminal pngcairo size " + to_string(width) + ", " + to_string(height) + " enhanced font 'Verdana,8'");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("set encoding utf8");
		gplot.cmd("set title \"" + title + "\"");
		gplot.cmd("set xlabel \"" + xlabel + "\"");
		gplot.cmd("set ylabel \"" + ylabel + "\"");
		gplot.cmd("set mxtics");
		gplot.cmd("set autoscale xfixmin");
		gplot.cmd("set autoscale xfixmax");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("plot \"./" + inputFile + "\" u 1:2 with points pointtype 7 pointsize 0.4 lc '#4169E1");
	}

	void drawChartWithSteps(string title, string xlabel, string ylabel, string inputFile, string outputFile, int width, int height) {
		Gnuplot gplot;

		gplot.cmd("set terminal pngcairo size " + to_string(width) + ", " + to_string(height) + " enhanced font 'Verdana,8'");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("set encoding utf8");
		gplot.cmd("set title \"" + title + "\"");
		gplot.cmd("set xlabel \"" + xlabel + "\"");
		gplot.cmd("set ylabel \"" + ylabel + "\"");
		gplot.cmd("set mxtics");
		gplot.cmd("set autoscale xfixmin");
		gplot.cmd("set autoscale xfixmax");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("plot \"./" + inputFile + "\" u 1:2 with steps lc '#4169E1");
	}

	void drawChartWithImpulses(string title, string xlabel, string ylabel, string inputFile, string outputFile, int width, int height) {
		Gnuplot gplot;

		gplot.cmd("set terminal pngcairo size " + to_string(width) + ", " + to_string(height) + " enhanced font 'Verdana,8'");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("set encoding utf8");
		gplot.cmd("set title \"" + title + "\"");
		gplot.cmd("set xlabel \"" + xlabel + "\"");
		gplot.cmd("set ylabel \"" + ylabel + "\"");
		gplot.cmd("set mxtics");
		gplot.cmd("set autoscale xfixmin");
		gplot.cmd("set autoscale xfixmax");
		gplot.cmd("set xzeroaxis");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("set style line 2 lc rgb '#4169E1' ps 1.4 pt 6");
		gplot.cmd("plot \"./" + inputFile + "\" u 1:2 with impulses lc '#4169E1', \"./" + inputFile + "\" u 1:2 with points ls 2");
	}

	void drawChartMultiPlot(string title, string title1, string title2, string title3, string title4, string xlabel, string ylabel, string outputFile, string inputFile1, string inputFile2, string inputFile3, string inputFile4, int width, int height, float xOffset = 0, float yOffset = 0) {
		Gnuplot gplot;

		gplot.cmd("set terminal pngcairo size " + to_string(width) + ", " + to_string(height) + " enhanced font 'Verdana,8'");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("set encoding utf8");
		gplot.cmd("set multiplot layout 4, 1 title \"" + title + "\" font \", 14\"");
		gplot.cmd("set title \"" + title + "\"");
		gplot.cmd("set xlabel \"" + xlabel + "\"");
		gplot.cmd("set ylabel \"" + ylabel + "\"");
		gplot.cmd("set tics font \", 6\"");
		gplot.cmd("set mxtics");
		gplot.cmd("set autoscale xfixmin");
		gplot.cmd("set autoscale xfixmax");
		gplot.cmd("set offsets " + to_string(xOffset) + ", " + to_string(xOffset) + ", " + to_string(yOffset) + ", " + to_string(yOffset));
		gplot.cmd("set autoscale xfixmin");
		gplot.cmd("set autoscale xfixmax");
		gplot.cmd("set style line 1 linetype 1 linewidth 1");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("set title \"" + title1 + "\"");
		gplot.cmd("plot \"./" + inputFile1 + "\" u 1:2 with lines linestyle 1 lc '#4169E1");
		gplot.cmd("set offsets 0, 0, 0, 0");
		gplot.cmd("set title \"" + title2 + "\"");
		gplot.cmd("plot \"./" + inputFile2 + "\" u 1:2 with lines linestyle 1 lc '#4169E1");
		gplot.cmd("set title \"" + title3 + "\"");
		gplot.cmd("plot \"./" + inputFile3 + "\" u 1:2 with lines linestyle 1 lc '#4169E1");
		gplot.cmd("set title \"" + title4 + "\"");
		gplot.cmd("plot \"./" + inputFile4 + "\" u 1:2 with lines linestyle 1 lc '#4169E1");
	}	
}
