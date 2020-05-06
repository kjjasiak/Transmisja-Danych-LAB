#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include "chart_utils.h"

using namespace chart_utils;
using std::vector;
using std::string;
using std::stringstream;


void drawSignalChart(vector<float> &xaxis, vector<float> &yaxis, int N, string plotType, string title, string filename, string xlabel, string ylabel, int width, int height, float xTics = 0, float xOffset = 0, float yOffset = 0) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << xaxis[i] << ";" << yaxis[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());

	drawChart(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height, xTics, xOffset, yOffset);

	st.str("");
}

void drawSignalsCharts(vector<signal> signals, int N, string filename, float xTics = 0, float xOffset = 0, float yOffset = 0) {
	stringstream st;
	
	for (int i = 0; i < signals.size(); i++) {
		for (int j = 0; j < N; j++) {
			st << signals[i].xaxis[j] << ";" << signals[i].yaxis[j] << ";\n";
		}

		dataToCsv("csv/" + signals[i].filename + ".csv", st.str());
		st.str("");
	}

	drawChartsToPDF("charts/" + filename + ".pdf", signals, xTics, xOffset, yOffset);
	st.str("");
}

void drawSignalsMultiCharts(vector<signal> signals, vector<int> plotsPerMPlot, int N, string title, string filename, float xTics = 0, float xOffset = 0, float yOffset = 0) {
	stringstream st;

	for (int i = 0; i < signals.size(); i++) {
		for (int j = 0; j < N; j++) {
			st << signals[i].xaxis[j] << ";" << signals[i].yaxis[j] << ";\n";
		}

		dataToCsv("csv/" + signals[i].filename + ".csv", st.str());
		st.str("");
	}

	drawChartsMultiplotToPDF(title, "charts/" + filename + ".pdf", signals, plotsPerMPlot, xTics, xOffset, yOffset);
	st.str("");
}
