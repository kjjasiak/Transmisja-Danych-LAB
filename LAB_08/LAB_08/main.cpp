#include <iostream>
#include <string>
#include <vector>
#include <bitset>
#include <sstream>
#include <fstream>
#include <time.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::bitset;
using std::stringstream;
using std::to_string;
using std::ofstream;


stringstream out;
string filename = "log.txt";

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

	for (int i = 0; i < output.size(); i++)
		cout << output[i] << " ";
	cout << endl;

	return output;
}

bitset<7> negateBit(bitset<7> bits, bool randomPos, int bitPos = 1) {
	bits = reverseBitset(bits);

	if (randomPos)
		bitPos = round(rand() % 7 + 1);
	out << " " << bitPos;

	bits[bitPos - 1] = !bits[bitPos - 1];
	bits = reverseBitset(bits);

	return bits;
}

bitset<8> negateBits(bitset<8> bits, bool randomPos, vector<int> bitPos = { 1, 1 }) {
	bits = reverseBitset(bits);

	for (int i = 0; i < bitPos.size(); i++) {
		if (randomPos)
			bitPos[i] = round(rand() % 7 + 1);
		out << " " << bitPos[i];

		bits[bitPos[i] - 1] = !bits[bitPos[i] - 1];
	}
	out << ";";

	bits = reverseBitset(bits);

	return bits;
}

vector<bitset<7>> hamming(vector<bitset<8>> binStream) {
	vector<bitset<4>> dVec;

	for (int i = 0; i < binStream.size(); i++) {
		bitset<8> streamByte = binStream[i];
		stringstream st;

		streamByte = reverseBitset(streamByte);

		for (int j = 0; j < 4; j++)
			st << streamByte[j];

		bitset<4> bits(st.str());
		dVec.push_back(bits);
		st.str("");

		for (int j = 4; j < 8; j++)
			st << streamByte[j];

		bitset<4> bits2(st.str());
		dVec.push_back(bits2);
		st.str("");
	}

	vector<bitset<7>> hBitsVec;

	for (int i = 0; i < dVec.size(); i++) {
		bitset<4> d = dVec[i];

		d = reverseBitset(d);

		bool p1 = (d[1 - 1] + d[2 - 1] + d[4 - 1]) % 2,
			 p2 = (d[1 - 1] + d[3 - 1] + d[4 - 1]) % 2,
			 p3 = (d[2 - 1] + d[3 - 1] + d[4 - 1]) % 2;

		stringstream st;
		st << p1 << p2 << d[1 - 1] << p3 << d[2 - 1] << d[3 - 1] << d[4 - 1];

		bitset<7> h(st.str());
		hBitsVec.push_back(h);
		st.str("");
	}

	return hBitsVec;
}

vector<bitset<8>> secded(vector<bitset<8>> binStream) {
	vector<bitset<4>> dVec;

	for (int i = 0; i < binStream.size(); i++) {
		bitset<8> streamByte = binStream[i];
		stringstream st;

		streamByte = reverseBitset(streamByte);

		for (int j = 0; j < 4; j++)
			st << streamByte[j];

		bitset<4> bits(st.str());
		dVec.push_back(bits);
		st.str("");

		for (int j = 4; j < 8; j++)
			st << streamByte[j];

		bitset<4> bits2(st.str());
		dVec.push_back(bits2);
		st.str("");
	}

	vector<bitset<8>> hBitsVec;

	for (int i = 0; i < dVec.size(); i++) {
		bitset<4> d = dVec[i];

		d = reverseBitset(d);

		bool p1 = (d[1 - 1] + d[2 - 1] + d[4 - 1]) % 2,
			 p2 = (d[1 - 1] + d[3 - 1] + d[4 - 1]) % 2,
			 p3 = (d[2 - 1] + d[3 - 1] + d[4 - 1]) % 2,
			 p4 = (p1 + p2 + d[1 - 1] + p3 + d[2 - 1] + d[3 - 1] + d[4 - 1]) % 2;

		stringstream st;
		st << p1 << p2 << d[1 - 1] << p3 << d[2 - 1] << d[3 - 1] << d[4 - 1] << p4;

		bitset<8> h(st.str());
		hBitsVec.push_back(h);

		st.str("");
	}

	return hBitsVec;
}

vector<bitset<4>> decodeHamming74(vector<bitset<7>> hamming74) {
	vector<bitset<4>> dVec;
	stringstream st;

	for (int i = 0; i < hamming74.size(); i++) {
		bitset<7> h = hamming74[i];
		h = reverseBitset(h);

		bool p1 = (h[1 - 1] + h[3 - 1] + h[5 - 1] + h[7 - 1]) % 2,
			 p2 = (h[2 - 1] + h[3 - 1] + h[6 - 1] + h[7 - 1]) % 2,
			 p3 = (h[4 - 1] + h[5 - 1] + h[6 - 1] + h[7 - 1]) % 2;

		int n = p1 * 1 + p2 * 2 + p3 * 4;

		if (n > 0) {
			h[n - 1] = !h[n - 1];
			out << "#" << i + 1 << ":" << endl;
			out << "Wystapil blad w transmisji, dokonano korekty na bicie nr " << n << endl;
		}

		st << h[3 - 1] << h[5 - 1] << h[6 - 1] << h[7 - 1];

		bitset<4> d(st.str());

		dVec.push_back(d);
		st.str("");
	}

	return dVec;
}

vector<bitset<4>> decodeSECDED(vector<bitset<8>> scd) {
	vector<bitset<4>> dVec;
	stringstream st;

	for (int i = 0; i < scd.size(); i++) {
		bitset<8> h = scd[i];
		h = reverseBitset(h);

		bool p4 = (h[1 - 1] + h[2 - 1] + h[3 - 1] + h[4 - 1] + h[5 - 1] + h[6 - 1] + h[7 - 1]) % 2;

		if (p4 != h[8 - 1]) {
			out << "#" << i + 1 << ":" << endl;
			out << "Wystapil blad, nastapi proba naprawy..." << endl;
		}

		bool p1 = (h[1 - 1] + h[3 - 1] + h[5 - 1] + h[7 - 1]) % 2,
			 p2 = (h[2 - 1] + h[3 - 1] + h[6 - 1] + h[7 - 1]) % 2,
			 p3 = (h[4 - 1] + h[5 - 1] + h[6 - 1] + h[7 - 1]) % 2;

		int n = p1 * pow(2, 0) + p2 * pow(2, 1) + p3 * pow(2, 2);

		if (n > 0) {
			h[n - 1] = !h[n - 1];

			out << "#" << i + 1 << ":" << endl;
			out << "Korekta na bicie nr " << n;

			if (p4 == h[8 - 1])
				out << " (XOR numerow bitow z bledami)";

			out << endl;

			p4 = (h[1 - 1] + h[2 - 1] + h[3 - 1] + h[4 - 1] + h[5 - 1] + h[6 - 1] + h[7 - 1]) % 2;

			if (p4 == h[8 - 1]) {
				out << "Wystapil 1 blad w transmisji, ktory naprawiono" << endl << endl;
			}
			else {
				out << "Wykryto co najmniej 2 bledy w transmisji" << endl;
				out << "Pakiet odrzucono, zdekodowany strumien jest niepoprawny - konieczna ponowna transmisja" << endl << endl;
			}
		}

		st << h[3 - 1] << h[5 - 1] << h[6 - 1] << h[7 - 1];

		bitset<4> d(st.str());
		dVec.push_back(d);
		st.str("");
	}

	return dVec;
}


int main() {
	// Strumienie binarne umieszczone w generowanym przez program pliku log.txt

	srand(time(NULL));

	vector<bitset<8>> str = strToBinStream("KOT", "bigEndian");

	out << "Strumien zrodlowy = ";
	for (int i = 0; i < str.size(); i++)
		out << "[" << str[i] << "]";
	out << endl;

	out << "----------------------------------------------------------------" << endl;

	vector<bitset<7>> hamming74 = hamming(str);

	out << "Hamming74 = ";

	for (int i = 0; i < hamming74.size(); i++)
		out << "[" << hamming74[i] << "]";

	out << endl;
	out << "----------------------------------------------------------------" << endl;

	out << "Zanegowane bity w kolejnych pakietach:";

	for (int i = 0; i < hamming74.size(); i++) {
		hamming74[i] = negateBit(hamming74[i], true);
	}
	out << endl << endl;

	out << "Hamming74 po negacji bitow = ";

	for (int i = 0; i < hamming74.size(); i++)
		out << "[" << hamming74[i] << "]";
	out << endl;

	out << "----------------------------------------------------------------" << endl;

	vector<bitset<4>> hammingDecoded = decodeHamming74(hamming74);

	out << endl << "Zdekodowany Hamming74 = ";

	for (int i = 0; i < hammingDecoded.size(); i++)
		out << "[" << hammingDecoded[i] << "]";

	out << endl;
	out << "----------------------------------------------------------------" << endl;

	vector<bitset<8>> scd = secded(str);
	out << "SECDED = ";

	for (int i = 0; i < scd.size(); i++)
		out << "[" << scd[i] << "]";
	
	out << endl;
	out << "----------------------------------------------------------------" << endl;

	out << "Zanegowane bity w kolejnych pakietach:";
	for (int i = 0; i < scd.size(); i++)
		scd[i] = negateBits(scd[i], true);
	out << endl << endl;

	out << "SECDED po negacji bitow = ";

	for (int i = 0; i < scd.size(); i++)
		out << "[" << scd[i] << "]";
	out << endl;

	out << "----------------------------------------------------------------" << endl;

	vector<bitset<4>> scdDecoded = decodeSECDED(scd);

	out << "Zdekodowany SECDED = ";

	for (int i = 0; i < scdDecoded.size(); i++)
		out << "[" << scdDecoded[i] << "]";
		
	out << endl;

	dataToFile(filename, out.str());
	out.str("");

	cout << endl << "Wygenerowano log z dzialania programu do pliku " << filename << endl;
	cout << "Nacisnij dowolny klawisz, aby zakonczyc..." << endl;

	getchar();
	return 0;
}