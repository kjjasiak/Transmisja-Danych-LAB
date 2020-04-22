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


void drawSignalChart(vector<float> &xaxis, vector<float> &yaxis, int N, string plotType, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << xaxis[i] << ";" << yaxis[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());

	if (plotType == "dots")
		drawChartWithDots(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);
	else
		drawChart(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSignalChart(vector<float> &xaxis, vector<float> &yaxis, int N, string plotType, string title, string filename, string xlabel, string ylabel, int width, int height, float xOffset, float yOffset) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << xaxis[i] << ";" << yaxis[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());

	if (plotType == "dots")
		drawChartWithDots(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);
	else
		drawChart(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height, xOffset, yOffset);

	st.str("");
}

void drawSignalMultiChart(vector<float> &xaxis1, vector<float> &xaxis2, vector<float> &xaxis3, vector<float> &xaxis4, vector<float> &yaxis1, vector<float> &yaxis2, vector<float> &yaxis3, vector<float> &yaxis4, int N, string plotType, string title, string title1, string title2, string title3, string title4, string outputFile, string filename1, string filename2, string filename3, string filename4, string xlabel, string ylabel, string xlabel2, string ylabel2, int width, int height) {
	stringstream st;
	float ymin = *min_element(yaxis1.begin(), yaxis1.end());
	float ymax = *max_element(yaxis1.begin(), yaxis1.end());

	// yaxis1
	for (int i = 0; i < N; i++)
		st << xaxis1[i] << ";" << yaxis1[i] << ";\n";

	dataToCsv("csv/" + filename1 + ".csv", st.str());
	st.str("");

	// yaxis2
	for (int i = 0; i < N; i++)
		st << xaxis2[i] << ";" << yaxis2[i] << ";\n";

	dataToCsv("csv/" + filename2 + ".csv", st.str());
	st.str("");

	// x3
	for (int i = 0; i < N; i++)
		st << xaxis3[i] << ";" << yaxis3[i] << ";\n";

	dataToCsv("csv/" + filename3 + ".csv", st.str());
	st.str("");

	// x4
	for (int i = 0; i < N; i++)
		st << xaxis4[i] << ";" << yaxis4[i] << ";\n";

	dataToCsv("csv/" + filename4 + ".csv", st.str());
	st.str("");

	drawChartMultiPlot(title, title1, title2, title3, title4, xlabel, ylabel, xlabel2, ylabel2, "charts/" + outputFile + ".png", "csv/" + filename1 + ".csv", "csv/" + filename2 + ".csv", "csv/" + filename3 + ".csv", "csv/" + filename4 + ".csv", width, height, 0, ymax * 0.1);

	st.str("");
}

void drawSignalMultiChart(vector<float> &xaxis1, vector<float> &xaxis2, vector<float> &xaxis3, vector<float> &xaxis4, vector<float> &xaxis5, vector<float> &yaxis1, vector<float> &yaxis2, vector<float> &yaxis3, vector<float> &yaxis4, vector<float> &yaxis5, int N, string plotType, string title, string title1, string title2, string title3, string title4, string title5, string outputFile, string filename1, string filename2, string filename3, string filename4, string filename5, string xlabel, string ylabel, string xlabel2, string ylabel2, int width, int height) {
	stringstream st;
	float ymin = *min_element(yaxis1.begin(), yaxis1.end());
	float ymax = *max_element(yaxis1.begin(), yaxis1.end());

	// yaxis1
	for (int i = 0; i < N; i++)
		st << xaxis1[i] << ";" << yaxis1[i] << ";\n";

	dataToCsv("csv/" + filename1 + ".csv", st.str());
	st.str("");

	// yaxis2
	for (int i = 0; i < N; i++)
		st << xaxis2[i] << ";" << yaxis2[i] << ";\n";

	dataToCsv("csv/" + filename2 + ".csv", st.str());
	st.str("");

	// x3
	for (int i = 0; i < N; i++)
		st << xaxis3[i] << ";" << yaxis3[i] << ";\n";

	dataToCsv("csv/" + filename3 + ".csv", st.str());
	st.str("");

	// x4
	for (int i = 0; i < N; i++)
		st << xaxis4[i] << ";" << yaxis4[i] << ";\n";

	dataToCsv("csv/" + filename4 + ".csv", st.str());
	st.str("");

	// x5
	for (int i = 0; i < N; i++)
		st << xaxis5[i] << ";" << yaxis5[i] << ";\n";

	dataToCsv("csv/" + filename5 + ".csv", st.str());
	st.str("");

	drawChartMultiPlot(title, title1, title2, title3, title4, title5, xlabel, ylabel, xlabel2, ylabel2, "charts/" + outputFile + ".png", "csv/" + filename1 + ".csv", "csv/" + filename2 + ".csv", "csv/" + filename3 + ".csv", "csv/" + filename4 + ".csv", "csv/" + filename5 + ".csv", width, height, 0, ymax * 0.1);

	st.str("");
}

void drawSpectrumChart(vector<float> &xaxis, vector<float> &yaxis, int N, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;

	for (int i = 0; i < N; i++)
		st << xaxis[i] << ";" << yaxis[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSpectrumPartChart(vector<float> &xaxis, vector<float> &yaxis, int N, float threshold, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;

	for (int i = 0; i < N / 2; i++) {
		if (yaxis[i] < threshold) {
			st << xaxis[i] << ";" << 0 << ";\n";
		}
		else {
			st << xaxis[i] << ";" << yaxis[i] << ";\n";
		}
	}

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSpectrumPartChart(vector<float> &xaxis, vector<float> &yaxis, int N, float threshold, string title, string filename, string xlabel, string ylabel, int width, int height, int freqCap) {
	stringstream st;

	for (int i = 0; i < N / 2; i++) {
		if (yaxis[i] < threshold) {
			st << xaxis[i] << ";" << 0 << ";\n";
		}
		else {
			st << xaxis[i] << ";" << yaxis[i] << ";\n";
		}

		if (xaxis[i] > freqCap)
			break;
	}

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}

void drawSpectrumPartdBChart(vector<float> &xaxis, vector<float> &yaxis, int N, string title, string filename, string xlabel, string ylabel, int width, int height) {
	stringstream st;
	int counter = 0;

	for (int i = 0; i < N / 2; i++)
		st << xaxis[i] << ";" << yaxis[i] << ";\n";

	dataToCsv("csv/" + filename + ".csv", st.str());
	drawChartWithImpulses(title, xlabel, ylabel, "csv/" + filename + ".csv", "charts/" + filename + ".png", width, height);

	st.str("");
}
