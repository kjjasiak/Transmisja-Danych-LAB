#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <string>

using std::cout;
using std::endl;
using std::string;

const int a = 8, b = 4, c = 3;

void rysujWykres(string filename) {
	
}

float wyroznik(float a, float b, float c) {
	return b * b - 4 * a*c;
}

void mZerowe(float a, float b, float c) {
	float t1, t2, delta;

	delta = wyroznik(a, b, c);

	if (delta == 0) {
		t1 = -(b / (2 * a));
		cout << "Jedno miejsce zerowe: " << t1;
	}
	else if (delta > 0) {
		t1 = (-b + sqrt(delta)) / (2 * a);
		t2 = (-b - sqrt(delta)) / (2 * a);
		cout << "Miejsca zerowe: t1 = " << t1 << "; t2 = " << t2;
	}
	else
		cout << "Brak miejsc zerowych" << endl;
}

float x(float t) {
	return a * pow(t, 2) + b * t + c;
}

float y(float t) {
	return 2 * pow(x(t), 2) + 12 * cos(t);
}

float z(float t) {
	return sin(2 * M_PI * 7 * t)*x(t) - 0.2*log10(abs(y(t)) + M_PI);
}

float u(float t) {
	return sqrt(abs(y(t)*y(t)*z(t))) - 1.8*sin(0.4*t*z(t)*x(t));
}

float v(float t) {
	if ((t >= 0) && (t < 0.22))
		return (1-7*t)*sin((2*M_PI*t*10)/(t+0.04));

	if ((t < 0.7) && (t >= 0.22))
		return 0.63*t*sin(125*t);

	if ((t >= 0.7) && (t <= 1.0))
		return pow(t, -0.662) + 0.77*sin(8*t);
}

void zad1() {
	mZerowe(a, b, c);
	const float t0 = -10, tN = 10, dt = 1.0 / 100;
	float tn = t0;
	//float t0 = -10, tN = 10, tn = t0, dt = 0.01;
	int n = 0;

	while (tn <= tN) {
		cout << tn << ";" << x(tn) << ";" << endl;
		n++;
		tn = t0 + (n*dt);
	}
}

void zad2() {
	const float t0 = 0, tN = 1, dt = 1.0 / 22050;
	float tn = t0;
	//float t0 = 0, tN = 1, tn = t0, dt = 1.0 / 22050;
	int n = 0;

	// y(t)
	while (tn <= tN) {
		cout << tn << ";" << y(tn) << ";" << endl;
		n++;
		tn = t0 + (n*dt);
	}

	// z(t)
	tn = t0;

	while (tn <= tN) {
		cout << tn << ";" << z(tn) << ";" << endl;
		n++;
		tn = t0 + (n*dt);
	}

	// u(t)
	tn = t0;

	while (tn <= tN) {
		cout << tn << ";" << u(tn) << ";" << endl;
		n++;
		tn = t0 + (n*dt);
	}
}

int main() {
	//zad1();
	//zad2();
	rysujWykres("cos.csv");

	//getchar();

	return 0;
}