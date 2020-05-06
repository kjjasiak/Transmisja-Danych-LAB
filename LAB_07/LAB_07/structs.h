#ifndef _STRUCTS
#define _STRUCTS

#include <vector>
#include <string>

using std::vector;
using std::string;

struct signal {
	vector<float> yaxis;
	vector<float> xaxis;
	string title;
	string filename;
	string xlabel;
	string ylabel;

	signal(vector<float> xaxis, vector<float> yaxis, string title, string filename, string xlabel, string ylabel) {
		this->xaxis = xaxis;
		this->yaxis = yaxis;
		this->title = title;
		this->filename = filename;
		this->xlabel = xlabel;
		this->ylabel = ylabel;
	}
};

#endif