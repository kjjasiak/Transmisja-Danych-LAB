#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <bitset>
#include <algorithm>
#include "structs.h"
#include "gnuplot_i.hpp"
#include <complex>

using std::cout;
using std::endl;
using std::string;
using std::bitset;
using std::vector;
using std::stringstream;
using std::ofstream;
using std::to_string;
using std::complex;

vector<float> getAmplSpectrum(vector<complex<float>> r, int N) {
	vector<float> M;

	for (int k = 0; k <= N - 1; k++) {
		float kTmp;

		kTmp = sqrt(real(r[k]) * real(r[k]) + imag(r[k]) * imag(r[k]));
		kTmp = (kTmp * 2) / N;

		M.push_back(kTmp);
	}

	return M;
}

vector<float> getAmplSpectrumBis(vector<float> M, float threshold, int N) {
	vector<float> MBis;

	for (int k = 0; k <= N - 1; k++) {
		if (M[k] < threshold)
			MBis.push_back(0);
		else
			MBis.push_back(10.0 * log10(M[k]));
	}

	return MBis;
}

// skala czestotliwosci
vector<float> getFreqScale(int fs, int N) {
	vector<float> f_k;

	for (int k = 0; k <= N - 1; k++) {
		f_k.push_back(k * (float(fs) / N));
	}

	return f_k;
}

// dft
vector<complex<float>> dft(vector<float> xn, int N) {
	vector<complex<float>> xk;

	for (int k = 0; k <= N - 1; k++) {
		complex<float> kTmp;

		for (int n = 0; n <= N - 1; n++) {
			complex<float> wNPow(cos((-2 * M_PI*k*n) / N), sin((-2 * M_PI*k*n) / N));

			kTmp += xn[n] * wNPow;
		}

		xk.push_back(kTmp);
	}

	return xk;
}

// idft
vector<float> idft(vector<complex<float>> xk, int N) {
	vector<float> xn;

	for (int n = 0; n <= N - 1; n++) {
		float nTmp = 0;

		for (int k = 0; k <= N - 1; k++) {
			complex<float> wNPow(cos((2 * M_PI*k*n) / N), sin((2 * M_PI*k*n) / N));

			nTmp += real(xk[k]) * real(wNPow) - imag(xk[k]) * imag(wNPow);
		}

		nTmp = nTmp / N;

		xn.push_back(nTmp);
	}

	return xn;
}


namespace chart_utils {
	void dataToCsv(string outputFile, string data) {
		ofstream of{ outputFile };

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

		gplot.cmd("set terminal pdf size 40 cm, 15 cm enhanced font 'Verdana,8'");
		gplot.cmd("set output \"./" + outputFile + "\"");
		gplot.cmd("set datafile separator \";\"");
		gplot.cmd("set encoding utf8");
		gplot.cmd("set size ratio 0.2");
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
			gplot.cmd("set xtics " + to_string(xTics) + " font 'Verdana,7'");
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

		gplot.cmd("set terminal pdf size 50 cm, 40 cm enhanced font 'Verdana,8'");
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
			gplot.cmd("set xtics " + to_string(xTics) + " font 'Verdana,6'");
		}

		int plotCount = 0;
		int tmp = 0;

		for (int i = 0; i < plotsPerMPlot.size(); i++) {
			gplot.cmd("set multiplot layout " + to_string(plotsPerMPlot[i]) + ", 1 title \"" + title + "\" font \", 12\"");

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


void drawSpectrumPartChart(vector<float> &t, vector<float> &x, int N, float threshold, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;

	for (int i = 0; i < N / 2; i++) {
		if (x[i] < threshold) {
			st << t[i] << ";" << 0 << ";\n";
		}
		else {
			st << t[i] << ";" << x[i] << ";\n";
		}
	}

	chart_utils::dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSpectrumPartChart(vector<float> &t, vector<float> &x, int N, float threshold, string title, string filename, string xlabel, string ylabel, int width, int height, int freqCap) {
	stringstream st;

	for (int i = 0; i < N / 2; i++) {
		if (x[i] < threshold) {
			st << t[i] << ";" << 0 << ";\n";
		}
		else {
			st << t[i] << ";" << x[i] << ";\n";
		}

		if (t[i] > freqCap)
			break;
	}

	chart_utils::dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSpectrumChart(vector<float> &t, vector<float> &x, int N, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << t[i] << ";" << x[i] << ";\n";

	chart_utils::dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSpectrumPartdBChart(vector<float> &t, vector<float> &x, int N, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;
	int counter = 0;

	for (int i = 0; i < N / 2; i++)
		st << t[i] << ";" << x[i] << ";\n";

	chart_utils::dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}


void dataToFile(string outputFile, string data) {
	ofstream of{ outputFile };

	if (of.is_open())
	{
		of << data;
		of.close();
	}
	else
		cout << "Nie udalo sie otworzyc: " << outputFile << endl;
}

bitset<8> reverseBitset(bitset<8> bits) {
	for (int i = 0; i < 4; ++i) {
		bool tmp = bits[i];
		bits[i] = bits[7 - i];
		bits[7 - i] = tmp;
	}

	return bits;
}

bitset<7> reverseBitset(bitset<7> bits) {
	for (int i = 0; i < 3; ++i) {
		bool tmp = bits[i];
		bits[i] = bits[6 - i];
		bits[6 - i] = tmp;
	}

	return bits;
}

bitset<4> reverseBitset(bitset<4> bits) {
	for (int i = 0; i < 2; ++i) {
		bool tmp = bits[i];
		bits[i] = bits[3 - i];
		bits[3 - i] = tmp;
	}

	return bits;
}

vector<bitset<8>> strToBinStream(string input, string option) {
	vector<bitset<8>> output;

	for (int i = 0; i < input.size(); i++) {
		bitset<8> tmp(input[i]);

		if (option == "littleEndian")
			reverseBitset(tmp);

		output.push_back(tmp);
	}

	return output;
}

string binStreamToStr(vector<bitset<8>> input, string option) {
	string output = "";

	for (int i = 0; i < input.size(); i++) {
		bitset<8> tmp(input[i]);

		char tmpChar = tmp.to_ulong();

		output += tmpChar;
	}

	return output;
}

vector<bitset<8>> signalToBitsets(vector<float> signal) {
	vector<bitset<8>> output;
	stringstream st;

	for (int i = 0; i < signal.size(); i = i + 8) {
		for (int j = i; j < i + 8; j++) {
			st << signal[j];
		}

		bitset<8> bits(st.str());
		output.push_back(bits);
		st.str("");
	}

	return output;
	return output;
}

string bitsetVectorToStr(vector<bitset<8>> input, bool splitHalf = false) {
	stringstream st;

	for (int i = 0; i < input.size(); i++) {
		for (int j = 7; j >= 0; j--) {
			st << input[i][j];
			
			if (splitHalf && (j % 4 == 0))
				st << " ";
		}
		
		if (!splitHalf)
			st << " ";
	}

	return st.str();
}



void drawSignalChart(vector<float> &xaxis, vector<float> &yaxis, int N, string plotType, string title, string filename, string xlabel, string ylabel, int width, int height, float xTics = 0, float xOffset = 0, float yOffset = 0) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << xaxis[i] << ";" << yaxis[i] << ";\n";

	chart_utils::dataToCsv("csv/" + filename + ".csv", st.str());

	chart_utils::drawChart(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height, xTics, xOffset, yOffset);

	st.str("");
}

void drawSignalsCharts(vector<signal> signals, int N, string filename, float xTics = 0, float xOffset = 0, float yOffset = 0) {
	stringstream st;

	for (int i = 0; i < signals.size(); i++) {
		for (int j = 0; j < N; j++) {
			st << signals[i].xaxis[j] << ";" << signals[i].yaxis[j] << ";\n";
		}

		chart_utils::dataToCsv("csv/" + signals[i].filename + ".csv", st.str());
		st.str("");
	}

	chart_utils::drawChartsToPDF("charts/" + filename + ".pdf", signals, xTics, xOffset, yOffset);
	st.str("");
}

void drawSignalsMultiCharts(vector<signal> signals, vector<int> plotsPerMPlot, int N, string title, string filename, float xTics = 0, float xOffset = 0, float yOffset = 0) {
	stringstream st;

	for (int i = 0; i < signals.size(); i++) {
		for (int j = 0; j < N; j++) {
			st << signals[i].xaxis[j] << ";" << signals[i].yaxis[j] << ";\n";
		}

		chart_utils::dataToCsv("csv/" + signals[i].filename + ".csv", st.str());
		st.str("");
	}

	chart_utils::drawChartsMultiplotToPDF(title, "charts/" + filename + ".pdf", signals, plotsPerMPlot, xTics, xOffset, yOffset);
	st.str("");
}
