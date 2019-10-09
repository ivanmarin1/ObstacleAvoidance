#pragma once
#include <string>

class BezierCurve
{
	//additional parameters for bezier curve through 4 points
	double c1, c2, c3, t1, t2;
	double a, b, c, d, e, f, g, h;
public:
	BezierCurve();
	~BezierCurve();
	//finding parameters to make sure curve passes through our points
	void bezier2(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3);
	//making bezier curve
	void bezierPath(double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3);
	void writeBezier(double x, double y, std::string fn);
};s

