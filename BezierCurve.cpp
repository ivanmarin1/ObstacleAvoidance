#include "BezierCurve.h"
#include "Source.cpp"
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

BezierCurve::BezierCurve()
{
	c1 = 0; c2 = 0; c3 = 0;
	t1 = 0; t2 = 0;
	a = 0; 
	b = 0; 
	c = 0; 
	d = 0;
	e = 0;
	f = 0;
	g = 0;
	h = 0;
}

BezierCurve::~BezierCurve()
{
}

void BezierCurve::bezier2(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3)
{
	float c1 = distanceEarth(x0, y0, x1, y1);
	float c2 = distanceEarth(x1, y1, x2, y2);
	float c3 = distanceEarth(x3, y3, x2, y2);
	float t1 = c1 / (c1 + c2 + c3);
	float t2 = (c1 + c2) / (c1 + c2 + c3);
	float a = t1 * (1 - t1)*(1 - t1) * 3;
	float b = (1 - t1)*t1*t1 * 3;
	float d = t2 * (1 - t2)*(1 - t2) * 3;
	float e = (1 - t2)*t2*t2 * 3;
	float c = x1 - (x0*pow(1 - t1, 3.0)) - (x3*pow(t1, 3));
	float f = x2 - (x0*pow(1 - t2, 3.0)) - (x3*pow(t2, 3));
	float g = y1 - (y0*pow(1 - t1, 3.0)) - (y3*pow(t1, 3));
	float h = y2 - (y0*pow(1 - t2, 3.0)) - (y3*pow(t2, 3));
	x2 = (c - a / d * f) / (b - a * e / d);
	x1 = (c - (b * x2)) / a;
	y2 = (g - a / d * h) / (b - a * e / d);
	y1 = (g - (b * y2)) / a;
	bezierPath(x0, x1, x2, x3, y0, y1, y2, y3);
}

void BezierCurve::bezierPath(double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3) {
	double xu = 0.0, yu = 0.0, u = 0.0;
	int i = 0;
	for (u = 0.0; u <= 1.0; u += 0.1)
	{
		xu = pow(1 - u, 3)*x0 + 3 * u*pow(1 - u, 2)*x1 + 3 * pow(u, 2)*(1 - u)*x2
			+ pow(u, 3)*x3;
		yu = pow(1 - u, 3)*y0 + 3 * u*pow(1 - u, 2)*y1 + 3 * pow(u, 2)*(1 - u)*y2
			+ pow(u, 3)*y3;
		cout << "Izracunao " << endl;
		writeBezier(xu, yu, "bezier.waypoints");
		cout << "Zapisao " << endl;
	}
}

void BezierCurve::writeBezier(double x, double y, string fn) {
	ofstream myfile;
	myfile.open(fn, ios::app);
	myfile.precision(10);
	if (!myfile) {
		cout << "Can't open file for writting waypoints!" << endl;
		return;
	}
	myfile << "QGC WPL 500" << endl;
	myfile << "0" << "\t0\t0\t16\t0\t0\t0\t0\t" << x << "\t" << y << "\t100.000000\t1" << endl;
	myfile.close();
}