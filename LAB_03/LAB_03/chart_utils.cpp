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
		gplot.cmd("set style line 1 linetype 1 linewidth 1");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("plot \"./" + inputFile + "\" u 1:2 with lines linestyle 1");
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
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("plot \"./" + inputFile + "\" u 1:2 with dots lc '#4169E1");
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
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("plot \"./" + inputFile + "\" u 1:2 with steps");
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
		gplot.cmd("set xzeroaxis");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("set style line 2 lc rgb '#4169E1' ps 1.4 pt 6");
		gplot.cmd("plot \"./" + inputFile + "\" u 1:2 with impulses lc '#4169E1', \"./" + inputFile + "\" u 1:2 with points ls 2");
	}
}