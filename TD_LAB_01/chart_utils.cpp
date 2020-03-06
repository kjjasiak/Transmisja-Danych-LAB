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

	void dataToCsv(string outputFile, string data) {
		ofstream of;

		of.open(outputFile, ios::out | ios::trunc);
		of << data;
		of.close();
	}

	void drawChart(string inputFile, string outputFile) {
		Gnuplot gplot;

		gplot.cmd("set terminal png");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("plot \"./" + inputFile + "\" u 1:2");
	}
}