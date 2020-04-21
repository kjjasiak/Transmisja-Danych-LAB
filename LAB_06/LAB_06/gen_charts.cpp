#pragma once

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

void drawSignalChart(vector<float> &t, vector<float> &x, int N, string plotType, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << t[i] << ";" << x[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());

	if (plotType == "dots")
		drawChartWithDots(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);
	else
		drawChart(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSignalChart(vector<float> &t, vector<float> &x, int N, string plotType, string title, string filename, string xlabel, string ylabel, int width, int height, float xOffset, float yOffset) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << t[i] << ";" << x[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());

	if (plotType == "dots")
		drawChartWithDots(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);
	else
		drawChart(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height, xOffset, yOffset);

	st.str("");
}

void drawSignalMultiChart(vector<float> &t1, vector<float> &t2, vector<float> &t3, vector<float> &t4, vector<float> &x1, vector<float> &x2, vector<float> &x3, vector<float> &x4, int N, string plotType, string title, string title1, string title2, string title3, string title4, string outputFile, string filename1, string filename2, string filename3, string filename4, string xlabel, string ylabel, int width, int height) {
	stringstream st;
	float ymin = *min_element(x1.begin(), x1.end());
	float ymax = *max_element(x1.begin(), x1.end());

	// x1
	for (int i = 0; i < N; i++)
		st << t1[i] << ";" << x1[i] << ";\n";

	dataToCsv("csv/" + filename1 + ".csv", st.str());
	st.str("");

	// x2
	for (int i = 0; i < N; i++)
		st << t2[i] << ";" << x2[i] << ";\n";

	dataToCsv("csv/" + filename2 + ".csv", st.str());
	st.str("");

	// x3
	for (int i = 0; i < N; i++)
		st << t3[i] << ";" << x3[i] << ";\n";

	dataToCsv("csv/" + filename3 + ".csv", st.str());
	st.str("");

	// x4
	for (int i = 0; i < N; i++)
		st << t4[i] << ";" << x4[i] << ";\n";

	dataToCsv("csv/" + filename4 + ".csv", st.str());
	st.str("");

	drawChartMultiPlot(title, title1, title2, title3, title4, xlabel, ylabel, "charts/" + outputFile + ".png", "csv/" + filename1 + ".csv", "csv/" + filename2 + ".csv", "csv/" + filename3 + ".csv", "csv/" + filename4 + ".csv", 1920, 1200, 0, ymax * 0.1);

	st.str("");
}

void drawSpectrumChart(vector<float> &t, vector<float> &x, int N, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << t[i] << ";" << x[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
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

	dataToCsv("csv/" + filename + ".csv", st.str());
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

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSpectrumPartdBChart(vector<float> &t, vector<float> &x, int N, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;
	int counter = 0;

	for (int i = 0; i < N / 2; i++)
		st << t[i] << ";" << x[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}
