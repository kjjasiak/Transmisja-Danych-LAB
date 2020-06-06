#include <iostream>
#include <algorithm>
#include <iomanip>
#include "helper_functions.h"
#include "coders_modulators.h"
#include "structs.h"


stringstream out;
string logFilename = "log.txt";

int main() {
	// Strumienie binarne umieszczone w generowanym przez program pliku log.txt

	// Etap 0.
	out << "Etap 0: " << endl;

	string source = "Kot.";
	out << "ASCII: " << source << endl;

	// Etap 1.
	out << endl << "Etap 1: " << endl;

	vector<bitset<8>> s1 = strToBinStream(source, "bigEndian");	
	out << "ASCII->BIN: " << bitsetVectorToStr(s1, true) << endl;

	// Etap 2.
	out << endl << "Etap 2: " << endl;
	vector<bitset<8>> scd = secded(s1);

	out << "BIN->SECDED: " << bitsetVectorToStr(scd) << endl;

	// Etap 3.
	out << endl << "Etap 3: " << endl;

	float Tb = 0.1;
	float f = 1.0 / Tb;
	float tN = Tb * (scd.size() * 8);
	float fi = 0;
	const int num_samples = 100 * f;

	string namePrefix = "LAB_09_", filename = "";
	stringstream title;

	vector<float> clk_xaxis, clk_yaxis,
		ttl_xaxis, ttl_yaxis, ttl_desampled_xaxis, ttl_desampled_yaxis,
		nrzi_xaxis, nrzi_yaxis, nrzi_desampled_xaxis, nrzi_desampled_yaxis,
		bami_xaxis, bami_yaxis, bami_desampled_xaxis, bami_desampled_yaxis,
		manch_xaxis, manch_yaxis, manch_desampled_xaxis, manch_desampled_yaxis;

	generateClockSignal(clk_xaxis, clk_yaxis, f, fi, tN);
	title << std::fixed << std::setprecision(2)
		<< "sygnal zegarowy (CLK)"
		<< ", f = " << f << " Hz";

	filename = namePrefix + "clk";
	signal clk(clk_xaxis, clk_yaxis, title.str(), filename, "t[s]", "A");
	title.str("");

	generateTTLSignal(ttl_xaxis, ttl_yaxis, scd, f);
	desampleSignal(ttl_xaxis, ttl_yaxis, ttl_desampled_xaxis, ttl_desampled_yaxis, num_samples);

	float ymax = *max_element(ttl_yaxis.begin(), ttl_yaxis.end());

	// wykres sygnalu TTL
	title << std::fixed << std::setprecision(2)
		<< "sygnal TTL"
		<< ", f = " << f << " Hz"
		<< ", T_b = " << Tb << " s";

	filename = namePrefix + "ttl";
	signal ttl(ttl_xaxis, ttl_yaxis, title.str(), filename, "t[s]", "A");
	title.str("");

	generateNRZISignal(nrzi_xaxis, nrzi_yaxis, clk_xaxis, clk_yaxis, ttl_xaxis, ttl_yaxis);
	desampleSignal(nrzi_xaxis, nrzi_yaxis, nrzi_desampled_xaxis, nrzi_desampled_yaxis, num_samples, 0.5 * num_samples);
	
	// wykres sygnalu NRZI
	title << std::fixed << std::setprecision(2)
		<< "sygnal NRZI"
		<< ", f = " << f << " Hz"
		<< ", T_b = " << Tb << " s";

	filename = namePrefix + "nrzi";
	signal nrzi(nrzi_xaxis, nrzi_yaxis, title.str(), filename, "t[s]", "A");
	title.str("");

	generateBAMISignal(bami_xaxis, bami_yaxis, clk_xaxis, clk_yaxis, ttl_xaxis, ttl_yaxis);
	desampleSignal(bami_xaxis, bami_yaxis, bami_desampled_xaxis, bami_desampled_yaxis, num_samples, 0.5 * num_samples);

	// wykres sygnalu BAMI
	title << std::fixed << std::setprecision(2)
		<< "sygnal BAMI"
		<< ", f = " << f << " Hz"
		<< ", T_b = " << Tb << " s";

	filename = namePrefix + "bami";
	signal bami(bami_xaxis, bami_yaxis, title.str(), filename, "t[s]", "A");
	title.str("");

	generateManchesterSignal(manch_xaxis, manch_yaxis, clk_xaxis, clk_yaxis, ttl_xaxis, ttl_yaxis);
	desampleSignal(manch_xaxis, manch_yaxis, manch_desampled_xaxis, manch_desampled_yaxis, num_samples, 0.5 * num_samples);

	// wykres sygnalu Manchester
	title << std::fixed << std::setprecision(2)
		<< "sygnal Manchester"
		<< ", f = " << f << " Hz"
		<< ", T_b = " << Tb << " s";

	filename = namePrefix + "manch";
	signal manch(manch_xaxis, manch_yaxis, title.str(), filename, "t[s]", "A");
	title.str("");

	// PDF z wykresami sygnalow
	vector<signal> signalsToPlot;

	signalsToPlot.push_back(clk);
	signalsToPlot.push_back(ttl);
	signalsToPlot.push_back(nrzi);
	signalsToPlot.push_back(bami);
	signalsToPlot.push_back(manch);

	drawSignalsMultiCharts(signalsToPlot, { 5 }, nrzi_yaxis.size(), "Modulacje - NRZI/BAMI/Manchester", namePrefix + "modulacje_wykresy", 0.1, 0, ymax * 0.1);

	out << "wykresy zapisane do pliku " << namePrefix + "modulacje_wykresy.pdf" << endl;

	// Etap 4.
	out << endl << "Etap 4: " << endl;
	vector<float> d_nrzi_xaxis, d_nrzi_yaxis, d_nrzi_desampled_xaxis, d_nrzi_desampled_yaxis;
	vector<float> d_bami_xaxis, d_bami_yaxis, d_bami_desampled_xaxis, d_bami_desampled_yaxis;
	vector<float> d_manch_xaxis, d_manch_yaxis, d_manch_desampled_xaxis, d_manch_desampled_yaxis;

	decodeNRZISignal(d_nrzi_xaxis, d_nrzi_yaxis, nrzi_xaxis, nrzi_yaxis, f);
	desampleSignal(d_nrzi_xaxis, d_nrzi_yaxis, d_nrzi_desampled_xaxis, d_nrzi_desampled_yaxis, num_samples, 0.5 * num_samples);

	out << "NRZI->demod.: ";
	for (int i = 0; i < d_nrzi_desampled_yaxis.size(); i++) {
		if ((i > 0) && (i % 8 == 0))
			out << " ";
		out << d_nrzi_desampled_yaxis[i];
	}
	out << endl;

	decodeBAMISignal(d_bami_xaxis, d_bami_yaxis, bami_xaxis, bami_yaxis);
	desampleSignal(d_bami_xaxis, d_bami_yaxis, d_bami_desampled_xaxis, d_bami_desampled_yaxis, num_samples, 0.5 * num_samples);

	out << "BAMI->demod.: ";
	for (int i = 0; i < d_bami_desampled_yaxis.size(); i++) {
		if ((i > 0) && (i % 8 == 0))
			out << " ";
		out << d_bami_desampled_yaxis[i];
	}
	out << endl;

	decodeManchesterSignal(d_manch_xaxis, d_manch_yaxis, manch_xaxis, manch_yaxis, f, tN);
	desampleSignal(d_manch_xaxis, d_manch_yaxis, d_manch_desampled_xaxis, d_manch_desampled_yaxis, num_samples, 0.75 * num_samples);

	out << "Man.->demod.: ";
	for (int i = 0; i < d_manch_desampled_yaxis.size(); i++) {
		if ((i > 0) && (i % 8 == 0))
			out << " ";
		out << d_manch_desampled_yaxis[i];
	}
	out << endl;
	
	// Etap 5.
	out << endl << "Etap 5: " << endl;

	vector<bitset<8>> nrzi_bitsets = signalToBitsets(d_nrzi_desampled_yaxis);
	vector<bitset<4>> nrzi_scd_decoded = decodeSECDED(nrzi_bitsets);

	out << "SECDED (NRZI)->BIN: ";
	for (int i = 0; i < nrzi_scd_decoded.size(); i++) {
		out << nrzi_scd_decoded[i] << " ";
	}
	out << endl;

	vector<bitset<8>> nrzi_scd_decoded_8bits;

	for (int i = 0; i < nrzi_scd_decoded.size(); i = i + 2) {
		stringstream st;

		st << nrzi_scd_decoded[i];
		st << nrzi_scd_decoded[i + 1];

		bitset<8> bits(st.str());
		nrzi_scd_decoded_8bits.push_back(bits);
		st.str("");
	}

	vector<bitset<8>> bami_bitsets = signalToBitsets(d_bami_desampled_yaxis);
	vector<bitset<4>> bami_scd_decoded = decodeSECDED(bami_bitsets);

	out << "SECDED (BAMI)->BIN: ";
	for (int i = 0; i < bami_scd_decoded.size(); i++) {
		out << bami_scd_decoded[i] << " ";
	}
	out << endl;

	vector<bitset<8>> bami_scd_decoded_8bits;

	for (int i = 0; i < bami_scd_decoded.size(); i = i + 2) {
		stringstream st;

		st << bami_scd_decoded[i];
		st << bami_scd_decoded[i + 1];

		bitset<8> bits(st.str());
		bami_scd_decoded_8bits.push_back(bits);
		st.str("");
	}

	vector<bitset<8>> manch_bitsets = signalToBitsets(d_manch_desampled_yaxis);
	vector<bitset<4>> manch_scd_decoded = decodeSECDED(manch_bitsets);

	out << "SECDED (Man.)->BIN: ";
	for (int i = 0; i < manch_scd_decoded.size(); i++) {
		out << manch_scd_decoded[i] << " ";
	}
	out << endl;

	vector<bitset<8>> manch_scd_decoded_8bits;

	for (int i = 0; i < manch_scd_decoded.size(); i = i + 2) {
		stringstream st;

		st << manch_scd_decoded[i];
		st << manch_scd_decoded[i + 1];

		bitset<8> bits(st.str());
		manch_scd_decoded_8bits.push_back(bits);
		st.str("");
	}

	// Etap 6.
	out << endl << "Etap 6: " << endl;

	string nrzi_s1r = binStreamToStr(nrzi_scd_decoded_8bits, "bigEndian");
	string bami_s1r = binStreamToStr(bami_scd_decoded_8bits, "bigEndian");
	string manch_s1r = binStreamToStr(manch_scd_decoded_8bits, "bigEndian");

	out << "BIN (NRZI)->ASCII: " << nrzi_s1r << endl;
	out << "BIN (BAMI)->ASCII: " << bami_s1r << endl;
	out << "BIN (Man.)->ASCII: " << manch_s1r << endl;

	dataToFile(logFilename, out.str());
	out.str("");

	cout << "Wygenerowano log z dzialania programu do pliku " << logFilename << endl;
	cout << "Nacisnij dowolny klawisz, aby zakonczyc..." << endl;

	getchar();
	return 0;
}