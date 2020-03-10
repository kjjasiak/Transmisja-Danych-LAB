#include <string>
#include <sstream>
#include <fstream>
#include "./gnuplot-cpp/gnuplot_i.hpp"

using std::string;
using std::stringstream;
using std::ofstream;
using std::ios;

namespace chart_utils {
	string generateData(float t0, float tN, float tn, float dt, float(*fn)(float)) {
		int n = 0;
		stringstream s;

		while (tn <= tN) {
			s << tn << ";" << fn(tn) << ";\n";
			n++;
			tn = t0 + (n*dt);
		}

		return s.str();
	}

	string generateData(float t0, float tN, float tn, float dt, int N, float(*fn)(float, int)) {
		int n = 0;
		stringstream s;

		while (tn <= tN) {
			s << tn << ";" << fn(tn, N) << ";\n";
			n++;
			tn = t0 + (n*dt);
		}

		return s.str();
	}

	void dataToCsv(string outputFile, string data) {
		ofstream of;

		of.open(outputFile, ios::out | ios::trunc);
		of << data;
		of.close();
	}

	void drawChart(string title, string xlabel, string ylabel, string inputFile, string outputFile) {
		Gnuplot gplot;

		gplot.cmd("set terminal pngcairo size 800, 600 enhanced font 'Verdana,9'");
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

	void drawChartWithDots(string title, string xlabel, string ylabel, string inputFile, string outputFile) {
		Gnuplot gplot;

		gplot.cmd("set terminal pngcairo size 800, 600 enhanced font 'Verdana,9'");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("set encoding utf8");
		gplot.cmd("set title \"" + title + "\"");
		gplot.cmd("set xlabel \"" + xlabel + "\"");
		gplot.cmd("set ylabel \"" + ylabel + "\"");
		//gplot.cmd("set style line 1 linetype 1 linewidth 1");
		gplot.cmd("unset key");
		gplot.cmd("set samples 1000");
		gplot.cmd("plot \"./" + inputFile + "\" u 1:2 with dots");
	}
}