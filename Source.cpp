// A C++ Program to implement Obstacle Avoidance
//_____________________________________________________________________
//_____________________________________________________________________
#include <iostream>
#include <list>
#include <set>
#include <stack>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
//#include "BezierCurve.h"
using namespace std;

#define _USE_MATH_DEFINES
#define earthRadiusMiles 3958.8
#define earthRadiusKilometers 6371.07103
#define M_PI 3.141592653589793238462643383279502884

// Creating a shortcut for int, int pair type 
typedef pair<int, int> Pair;

// Creating a shortcut for pair<int, pair<int, int>> type 
typedef pair<double, pair<int, int>> pPair;

// A structure to hold the neccesary parameters 
struct cell {
	// Row and Column index of its parent 
	// Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1 
	int parent_i, parent_j;
	// f = g + h 
	double f, g, h;
	double lat, lon;
	int fly;
	int id, FirstCol, SecondCol, ThirdCol;
	int ForthCol, FiveCol, SixCol, SevenCol, LastCol;
	double altitude;
	cell() {
		parent_i = 0;
		parent_j = 0;
		f = 0;
		g = 0;
		h = 0;
		fly = 1;
		lat = 0;
		lon = 0;
		id = 0;
		FirstCol = 0;
		SecondCol = 0;
		ThirdCol = 16;
		ForthCol = 0;
		FiveCol = 0;
		SixCol = 0;
		SevenCol = 0;
		altitude = 100.00;
		LastCol = 1;
	}
};

struct waypoint :cell {
	int row, col;
	double lat, lon;
	waypoint() {
		row = 0;
		col = 0;
		lat = 0;
		lon = 0;
	}
};

struct obstacle :waypoint {
	double length;
	obstacle() {
		length = 0;
	}
};

struct path :waypoint {
	
};

//Opening blank files instead of appending data to the old ones
void resetFiles(string filename) {
	ofstream myfile;
	myfile.open(filename);
	if (filename == "files/path.waypoints")
		myfile << "QGC WPL 500" << endl;
	myfile.close();
}

void countElementsInFile(string fn, int& Num) {
	string word[12];
	string str;
	ifstream myfile;
	int i = 0;
	myfile.open(fn);
	if (!myfile) {
		cout << fn << " file is missing! " << endl;
		system("pause");
		exit(0);
	}
	else {
		myfile.exceptions(ifstream::badbit); // No need to check failbit
		try {
			while (getline(myfile, str)) {
				stringstream ss(str);
				ss >> word[0] >> word[1] >> word[2] >> word[3]
					>> word[4] >> word[5] >> word[6] >> word[7] >> word[8] >> word[9]
					>> word[10] >> word[11];
				if (word[0] == "QGC" || word[0] == "QGC " || word[0] == " " || word[1] == " ") {
				}
				else
					i++;
			}
			myfile.close();
			Num = i;
			// for line-oriented input use file.getline(s)
		}
		catch (const ifstream::failure& e) {
			cout << "Exception opening/reading file" << endl;
		}
	}
}

void readBoundaryFile(cell *boundaryPoint, string fn) {
	string line, word, word2;
	ifstream myfile;
	myfile.open(fn);
	if (!myfile) {
		cout << "Boundary file is missing! " << endl;
		system("pause");
		exit(0);
	}
	else {
		int i = 0;
		myfile.exceptions(ifstream::badbit); // No need to check failbit
		try {
			while (myfile >> word >> word2) {
				//cout << word << " " << word2 << endl;
				boundaryPoint[i].lat = atof(word.c_str());
				boundaryPoint[i].lon = atof(word2.c_str());
				i++;
			}
			myfile.close();
			// for line-oriented input use file.getline(s)
		}
		catch (const ifstream::failure& e) {
			cout << "Exception opening/reading file";
		}
	}
}

void readWaypointsFile(waypoint *waypoints, string fn) {
	ifstream myfile;
	string word[12];
	string str;
	myfile.open(fn);
	if (!myfile) {
		cout << fn << " file is missing! " << endl;
		system("pause");
		exit(0);
	}
	else {
		int i = 0;
		myfile.exceptions(ifstream::badbit); // No need to check failbit
		try {
			while (getline(myfile, str)) { // Read the entire line
				stringstream ss(str);
				ss >> word[0] >> word[1] >> word[2] >> word[3]
					>> word[4] >> word[5] >> word[6] >> word[7] >> word[8] >> word[9]
					>> word[10] >> word[11];
				if (word[0] == "QGC" || word[0] == "QGC ") {

				}
				else if (word[3] == "22" || (stoi(word[8]) == 0 && stoi(word[9]) == 0)) {
					waypoints[i].id = atof(word[0].c_str());
					waypoints[i].FirstCol = atof(word[1].c_str());
					waypoints[i].SecondCol = atof(word[2].c_str());
					waypoints[i].ThirdCol = atof(word[3].c_str());
					waypoints[i].ForthCol = atof(word[4].c_str());
					waypoints[i].FiveCol = atof(word[5].c_str());
					waypoints[i].SixCol = atof(word[6].c_str());
					waypoints[i].SevenCol = atof(word[7].c_str());
					waypoints[i].lat = waypoints[i - 1].lat;
					waypoints[i].lon = waypoints[i - 1].lon;
					waypoints[i].altitude = atof(word[10].c_str());
					waypoints[i].LastCol = atof(word[11].c_str());
					i++;
				}
				else {
					waypoints[i].id = atof(word[0].c_str());
					waypoints[i].FirstCol = atof(word[1].c_str());
					waypoints[i].SecondCol = atof(word[2].c_str());
					waypoints[i].ThirdCol = atof(word[3].c_str());
					waypoints[i].ForthCol = atof(word[4].c_str());
					waypoints[i].FiveCol = atof(word[5].c_str());
					waypoints[i].SixCol = atof(word[6].c_str());
					waypoints[i].SevenCol = atof(word[7].c_str());
					waypoints[i].lat = atof(word[8].c_str());
					waypoints[i].lon = atof(word[9].c_str());
					waypoints[i].altitude = atof(word[10].c_str());
					waypoints[i].LastCol = atof(word[11].c_str());
					i++;
				}
			}
			myfile.close();
			// for line-oriented input use file.getline(s)
		}
		catch (const ifstream::failure& e) {
			cout << "Exception opening/reading file";
		}
	}
}

void readPathFile(waypoint *waypoints, string fn, path *path1, int el) {
	ifstream myfile;
	string word[12];
	string str;
	myfile.open(fn);
	if (!myfile) {
		cout << fn << " file is missing! " << endl;
		system("pause");
		exit(0);
	}
	else {
		int i = 0;
		myfile.exceptions(ifstream::badbit); // No need to check failbit
		try {
			while (getline(myfile, str)) { // Read the entire line
				stringstream ss(str);
				ss >> word[0] >> word[1] >> word[2] >> word[3]
					>> word[4] >> word[5] >> word[6] >> word[7] >> word[8] >> word[9]
					>> word[10] >> word[11];
				if (word[0] == "QGC" || word[0] == "QGC ") {

				}
				else if (word[3] == "22" || (stoi(word[8]) == 0 && stoi(word[9]) == 0)) {
					int temp = i;
					path1[i].id = atof(word[0].c_str());
					path1[i].FirstCol = atof(word[1].c_str());
					path1[i].SecondCol = atof(word[2].c_str());
					path1[i].ThirdCol = atof(word[3].c_str());
					path1[i].ForthCol = atof(word[4].c_str());
					path1[i].FiveCol = atof(word[5].c_str());
					path1[i].SixCol = atof(word[6].c_str());
					path1[i].SevenCol = atof(word[7].c_str());
					path1[i].lat = waypoints[--temp].lat;
					path1[i].lon = waypoints[temp].lon;
					path1[i].altitude = atof(word[10].c_str());
					path1[i].LastCol = atof(word[11].c_str());
					i++;
				}
			}
			el = i;
			myfile.close();
			// for line-oriented input use file.getline(s)
		}
		catch (const ifstream::failure& e) {
			cout << "Exception opening/reading file";
		}
	}
}

void readObstaclesFile(obstacle *obstacles, string fn) {
	string line, word, word2, word3;
	ifstream myfile;
	myfile.open(fn);
	if (!myfile) {
		cout << "Obstacles file is missing! " << endl;
		system("pause");
		exit(0);
	}
	else {
		int i = 0;
		myfile.exceptions(ifstream::badbit); // No need to check failbit
		try {
			while (myfile >> word >> word2 >> word3) {
				//cout << word << " " << word2 << endl;
				obstacles[i].lat = atof(word.c_str());
				obstacles[i].lon = atof(word2.c_str());
				obstacles[i].length = atof(word3.c_str());
				i++;
			}
			myfile.close();
			// for line-oriented input use file.getline(s)
		}
		catch (const ifstream::failure& e) {
			cout << "Exception opening/reading file";
		}
	}
}

//writing grid, only lat and lon, not other parameters
void writeToFile(cell **grid, int nodeNum, string fn) {
	ofstream myfile;
	myfile.open(fn);
	myfile.precision(10);
	if (!myfile) {
		cout << "Can't open file for writting coordinates!" << endl;
		return;
	}
	myfile << "QGC WPL 500" << endl;
	int count = 0;
	for (int i = 0; i < nodeNum; i++) {
		for (int j = 0; j < nodeNum; j++) {
			myfile << count << "\t0\t0\t16\t0\t0\t0\t0\t" << grid[i][j].lat << "\t" << grid[i][j].lon << "\t100.000000\t1" << endl;
			count++;
		}
	}
	myfile.close();
}

void writeWPToFile(waypoint *waypoints, int wpNum, string fn) {
	ofstream myfile;
	myfile.open(fn);
	myfile.precision(10);
	if (!myfile) {
		cout << "Can't open file for writting waypoints!" << endl;
		return;
	}
	myfile << "QGC WPL 500" << endl;
	for (int i = 0; i < wpNum; i++) {
		myfile << waypoints[i].id << "\t" << waypoints[i].FirstCol << "\t" << waypoints[i].SecondCol << "\t" <<
			waypoints[i].ThirdCol << "\t" << waypoints[i].ForthCol << "\t" << waypoints[i].FiveCol << "\t" <<
			waypoints[i].SixCol << "\t" << waypoints[i].SevenCol << "\t" << waypoints[i].lat << "\t" <<
			waypoints[i].lon << "\t" << waypoints[i].altitude << "\t" << waypoints[i].LastCol << endl;
	}
	myfile.close();
}

void writeIndividualWPToFile(waypoint *waypoints, int i, string fn, int WPNum) {
	ofstream myfile;
	myfile.open(fn, ofstream::app);
	myfile.precision(10);
	if (!myfile) {
		cout << "Can't open file for writting waypoints!" << endl;
		return;
	}
	
	myfile << waypoints[i].id << "\t" << waypoints[i].FirstCol << "\t" << waypoints[i].SecondCol << "\t" <<
		waypoints[i].ThirdCol << "\t" << waypoints[i].ForthCol << "\t" << waypoints[i].FiveCol << "\t" <<
		waypoints[i].SixCol << "\t" << waypoints[i].SevenCol << "\t" << waypoints[i].lat << "\t" <<
		waypoints[i].lon << "\t" << waypoints[i].altitude << "\t" << waypoints[i].LastCol << endl;
	if (i == WPNum - 2)
		myfile << waypoints[i + 1].id << "\t" << waypoints[i + 1].FirstCol << "\t" << waypoints[i + 1].SecondCol << "\t" <<
		waypoints[i + 1].ThirdCol << "\t" << waypoints[i + 1].ForthCol << "\t" << waypoints[i + 1].FiveCol << "\t" <<
		waypoints[i + 1].SixCol << "\t" << waypoints[i + 1].SevenCol << "\t" << waypoints[i + 1].lat << "\t" <<
		waypoints[i + 1].lon << "\t" << waypoints[i + 1].altitude << "\t" << waypoints[i + 1].LastCol << endl;
	myfile.close();
}

void writeOBToFile(cell **grid, int nodeNum, string fn) {
	ofstream myfile;
	myfile.open(fn);
	myfile.precision(10);
	if (!myfile) {
		cout << "Can't open file for writting obstacles!" << endl;
		return;
	}
	myfile << "QGC WPL 500" << endl;
	int count = 0;
	for (int i = 0; i < nodeNum; i++) {
		for (int j = 0; j < nodeNum; j++) {
			if (grid[i][j].fly == 0) {
				myfile << count << "\t0\t0\t16\t0\t0\t0\t0\t" << grid[i][j].lat << "\t" << grid[i][j].lon << "\t100.000000\t1" << endl;
				count++;
			}
		}
	}
	myfile.close();
}

void writeBezier(double x, double y, string fn) {
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

void printWaypoints(waypoint *waypoints, int nodeNum) {
	cout.precision(10);
	for (int i = 0; i < nodeNum; i++) {
		cout << waypoints[i].lat << " " << waypoints[i].lon << "\t" << endl;
	}
}

void writeGPSGrid(cell **grid, int nodeNum, string fn) {
	ofstream myfile;
	myfile.open(fn);
	myfile.precision(10);
	if (!myfile) {
		cout << "Can't open file for writting obstacles!" << endl;
		return;
	}
	myfile << "QGC WPL 500" << endl;
	int count = 0;
	for (int i = 0; i < nodeNum; i++) {
		for (int j = 0; j < nodeNum; j++) {
			if(count % 100 == 0)
				myfile << count << "\t0\t0\t16\t0\t0\t0\t0\t" << grid[i][j].lat << "\t" << grid[i][j].lon << "\t100.000000\t1" << endl;
			count++;
		}
		cout << endl;
	}
	myfile.close();
}

// Method to compare which one is the more close. 
// We find the closest by taking the difference 
// between the target and both values. It assumes 
// that val2 is greater than val1 and target lies 
// between these two. 
double getClosest(double val1, double val2, double target) {
	if (abs(target - val1) <= abs(val2 - target))
		return val1;
	else
		return val2;
}

// Returns element closest to target in arr[] 
//orderNum ---- 1 ascending, 2 descending
double findClosest(double arr[], int n, double target, int &el, int orderNum) {
	double temp = 0;
	int order = orderNum;
	if (order == 1) {
		// Corner cases 
		if (target >= arr[0]) {
			el = 0;
			return arr[0];
		}
		if (target <= arr[n - 1]) {
			el = n - 1;
			return arr[n - 1];
		}
		// Doing binary search 
		int i = 0, j = n, mid = 0;
		while (i < j) {
			mid = (i + j) / 2;

			if (arr[mid] == target) {
				el = mid;
				return arr[mid];
			}
			/* If target is greater than array element,
				then search in left */
			if (target > arr[mid]) {

				// If target is greater than previous 
				// to mid, return closest of two 
				if (mid > 0 && target < arr[mid - 1]) {
					temp = getClosest(arr[mid - 1], arr[mid], target);
					if (temp == arr[mid - 1]) el = mid - 1;
					else el = mid;
					return temp;
				}
				/* Repeat for left half */
				j = mid;
			}

			// If target is less than mid 
			else {
				if (mid < n - 1 && target > arr[mid + 1]) {
					temp = getClosest(arr[mid], arr[mid + 1], target);
					if (temp == arr[mid]) el = mid;
					else el = mid + 1;
					return temp;
				}
				// update i 
				i = mid + 1;
			}
		}
		return arr[mid];
	}
	else if (order == 2) {
		// Corner cases 
		if (target <= arr[0]) {
			el = 0;
			return arr[0];
		}
		if (target >= arr[n - 1]) {
			el = n - 1;
			return arr[n - 1];
		}
		// Doing binary search 
		int i = 0, j = n, mid = 0;
		while (i < j) {
			mid = (i + j) / 2;

			if (arr[mid] == target) {
				el = mid;
				return arr[mid];
			}
			/* If target is less than array element,
				then search in left */
			if (target < arr[mid]) {

				// If target is greater than previous 
				// to mid, return closest of two 
				if (mid > 0 && target > arr[mid - 1]) {
					temp = getClosest(arr[mid - 1], arr[mid], target);
					if (temp == arr[mid - 1]) el = mid - 1;
					else el = mid;
					return temp;
				}
				/* Repeat for left half */
				j = mid;
			}
			// If target is greater than mid 
			else {
				if (mid < n - 1 && target < arr[mid + 1]) {
					temp = getClosest(arr[mid], arr[mid + 1], target);
					if (temp == arr[mid]) el = mid;
					else el = mid + 1;
					return temp;
				}
				// update i
				i = mid + 1;
			}
		}
		return arr[mid];
	}
	else {
		cout << "Wrong order parameters sent!" << endl;
		return 0;
	}
	// Only single element left after search 
}

double mToFeet(double meters) {
	return (meters * 1000 * 3.2808399);
}

// This function converts decimal degrees to radians
double deg2rad(double deg) {
	return (deg * M_PI / 180);
}

//  This function converts radians to decimal degrees
double rad2deg(double rad) {
	return (rad * 180 / M_PI);
}

double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
	double lat1r, lon1r, lat2r, lon2r, u, v;
	lat1r = deg2rad(lat1d);
	lon1r = deg2rad(lon1d);
	lat2r = deg2rad(lat2d);
	lon2r = deg2rad(lon2d);
	u = sin((lat2r - lat1r) / 2);
	v = sin((lon2r - lon1r) / 2);
	return 2.0 * earthRadiusKilometers * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

double calculateMaxDistance(cell *boundary, int boundNum, cell &boundaryPoint1, cell &boundaryPoint2) {
	double maxDist = 0;
	double distance = 0;
	for (int i = 0; i < boundNum; i++) {
		for (int j = 0; j < boundNum; j++) {
			distance = distanceEarth(boundary[i].lat, boundary[i].lon, boundary[j].lat, boundary[j].lon);
			if (distance > maxDist) {
				maxDist = distance;
				boundaryPoint1 = boundary[i];
				boundaryPoint2 = boundary[j];
			}
		}
	}
	return maxDist;
}

double bearing(double lat, double lon, double lat2, double lon2) {

	double teta1 = deg2rad(lat);
	double teta2 = deg2rad(lat2);
	double delta1 = deg2rad(lat2 - lat);
	double delta2 = deg2rad(lon2 - lon);

	//==================Heading Formula Calculation================//

	double y = sin(delta2) * cos(teta2);
	double x = cos(teta1)*sin(teta2) - sin(teta1)*cos(teta2)*cos(delta2);
	double brng = atan2(y, x);
	brng = rad2deg(brng);// radians to degrees
	brng = (((int)brng + 360) % 360);

	return brng;
}

cell fillGPSCoordinate(cell &home, double dist, double theta) {
	double theta2 = theta * M_PI / 180;							//theta facing north = 0, east = 90
	double homeLatRad, homeLonRad;
	cell nextNode;
	homeLatRad = deg2rad(home.lat);
	homeLonRad = deg2rad(home.lon);
	double lat2 = asin(sin(homeLatRad) * cos(dist / (earthRadiusMiles)) + cos(homeLatRad) * sin(dist / (earthRadiusMiles)) * cos(theta2));
	double lon2 = homeLonRad + atan2(sin(theta2) * sin(dist / (earthRadiusMiles)) * cos(homeLatRad), cos(dist / (earthRadiusMiles)) - sin(homeLatRad) * sin(lat2));
	nextNode.lat = rad2deg(lat2);
	nextNode.lon = rad2deg(lon2);
	nextNode.fly = 1;
	return nextNode;
}

// A Utility Function to check whether given cell (row, col) 
// is a valid cell or not. 
bool isValid(int row, int col, int nodeNum)
{
	// Returns true if row number and column number 
	// is in range 
	return (row >= 0) && (row < nodeNum) &&
		(col >= 0) && (col < nodeNum);
}

// A Utility Function to check whether the given cell is 
// blocked or not 
bool isUnBlocked(cell **grid, int row, int col)
{
	// Returns true if the cell is not blocked else false 
	if (grid[row][col].fly == 1)
		return (true);
	else
		return (false);
}

// A Utility Function to check whether destination cell has 
// been reached or not 
bool isDestination(int row, int col, Pair dest)
{
	if (row == dest.first && col == dest.second)
		return (true);
	else
		return (false);
}

bool isWithinAngle(double heading, double heading2, int delta) {
	double temp = heading - delta;
	double temp2 = heading + delta;
	if (temp < 0) {
		temp = temp + 360;
		if ((heading2 >= temp && heading2 < 360) || (heading2 >= 0 && heading2 <= temp2))
			return true;
		else
			return false;
	}	
	else if (temp2 > 360){
		temp2 = temp2 - 360;
		if ((heading2 >= temp && heading2 < 360) || (heading2 >= 0 && heading2 <= temp2))
			return true;
		else
			return false;
	}
	else {
		if (heading2 >= temp && heading2 <= temp2)
			return true;
		else
			return false;
	}
}

// A Utility Function to calculate the 'h' heuristics. 
double calculateHValue(int row, int col, Pair dest)
{
	// Return using the distance formula 
	return ((double)sqrt((row - dest.first)*(row - dest.first)
		+ (col - dest.second)*(col - dest.second)));
}

void addWaypoint(cell **grid, int nodeNum, waypoint *waypoint, int el) {
	double *lat = new double[nodeNum];
	double *latLargest = new double[nodeNum];
	double *lon = new double[nodeNum];
	double *lonLargest = new double[nodeNum];
	int row = 0, col = 0;
	for (int i = 0; i < nodeNum; i++) {
		for (int j = 0; j < nodeNum; j++) {
			lat[j] = grid[i][j].lat;
		}
		latLargest[i] = findClosest(lat, nodeNum, waypoint[el].lat, row, 1);	//1 - descending
	}
	latLargest[0] = findClosest(latLargest, nodeNum, waypoint[el].lat, row, 1);
	for (int i = 0; i < nodeNum; i++) {
		lon[i] = grid[row][i].lon;
	}
	lonLargest[0] = findClosest(lon, nodeNum, waypoint[el].lon, col, 2);		//2 - ascending
	waypoint[el].row = row;
	waypoint[el].col = col;
	grid[row][col].id = waypoint[el].id;
	grid[row][col].FirstCol = waypoint[el].FirstCol;
	grid[row][col].SecondCol = waypoint[el].SecondCol;
	grid[row][col].ThirdCol = waypoint[el].ThirdCol;
	grid[row][col].ForthCol = waypoint[el].ForthCol;
	grid[row][col].FiveCol = waypoint[el].FiveCol;
	grid[row][col].SixCol = waypoint[el].SixCol;
	grid[row][col].SevenCol = waypoint[el].SevenCol;
	grid[row][col].altitude = waypoint[el].altitude;
	grid[row][col].LastCol = waypoint[el].LastCol;

	delete[] lat;
	delete[] latLargest;
	delete[] lon;
	delete[] lonLargest;
}

void FindObstacleRowCol(cell **grid, int nodeNum, obstacle *obstacle, int el) {
	double *lat = new double[nodeNum];
	double *latLargest = new double[nodeNum];
	double *lon = new double[nodeNum];
	double *lonLargest = new double[nodeNum];
	int row = 0, col = 0;
	for (int i = 0; i < nodeNum; i++) {
		for (int j = 0; j < nodeNum; j++) {
			lat[j] = grid[i][j].lat;
		}
		latLargest[i] = findClosest(lat, nodeNum, obstacle[el].lat, row, 1);	//1 - descending
	}
	latLargest[0] = findClosest(latLargest, nodeNum, obstacle[el].lat, row, 1);
	for (int i = 0; i < nodeNum; i++) {
		lon[i] = grid[row][i].lon;
	}
	lonLargest[0] = findClosest(lon, nodeNum, obstacle[el].lon, col, 2);		//2 - ascending
	obstacle[el].row = row;
	obstacle[el].col = col;
	delete[] lat;
	delete[] latLargest;
	delete[] lon;
	delete[] lonLargest;
}

void addObstacle(waypoint *wp, int wpNum, int nodeNum, obstacle *obstacles, int el, cell **grid, int obNum) {
	int row = 0, col = 0, length = 0;
	//findClosestWP(wp, wpNum, obstacles, el, w1, w2, obNum);
	FindObstacleRowCol(grid, nodeNum, obstacles, el);
	length = obstacles[el].length / 10;
	length++;
	row = obstacles[el].row;
	col = obstacles[el].col;
	if (isValid(row, col, nodeNum) == false)
		grid[row][col].fly = 0;
	for (int i = 1; i < length; i++) {
		//left, right, up and down
		if (isValid(row - i, col, nodeNum) == true) grid[row - i][col].fly = 0;
		if (isValid(row, col - i, nodeNum) == true) grid[row][col - i].fly = 0;
		if (isValid(row + i, col, nodeNum) == true) grid[row + i][col].fly = 0;
		if (isValid(row, col + i, nodeNum) == true) grid[row][col + i].fly = 0;
		//diagonal
		if (isValid(row + i, col + i, nodeNum) == true) grid[row + i][col + i].fly = 0;
		if (isValid(row + i, col - i, nodeNum) == true) grid[row + i][col - i].fly = 0;
		if (isValid(row - i, col + i, nodeNum) == true) grid[row - i][col + i].fly = 0;
		if (isValid(row - i, col - i, nodeNum) == true) grid[row - i][col - i].fly = 0;
		//between diagonal and straight lines
		if (isValid(row + i, col + i - 1, nodeNum) == true) grid[row + i][col + i - 1].fly = 0;
		if (isValid(row + i, col - i - 1, nodeNum) == true) grid[row + i][col - i - 1].fly = 0;
		if (isValid(row - i + 1, col + i, nodeNum) == true) grid[row - i + 1][col + i].fly = 0;
		if (isValid(row - i - 1, col - i, nodeNum) == true) grid[row - i - 1][col - i].fly = 0;
	}
}

// A Utility Function to trace the path from the source 
// to destination 
/*void tracePath(cell **cellDetails, Pair dest, int nodeNum, cell **grid, string pathFN)
{
	printf("\nThe Path is ");
	int row = dest.first;
	int col = dest.second;
	ofstream myfile;
	myfile.open(pathFN, ofstream::app);
	myfile.precision(10);
	myfile << "QGC WPL 500" << endl;
	int count = 0;

	stack<Pair> Path;

	while (!(cellDetails[row][col].parent_i == row
		&& cellDetails[row][col].parent_j == col))
	{
		Path.push(make_pair(row, col));
		int temp_row = cellDetails[row][col].parent_i;
		int temp_col = cellDetails[row][col].parent_j;
		row = temp_row;
		col = temp_col;
	}

	Path.push(make_pair(row, col));
	while (!Path.empty())
	{
		pair<int, int> p = Path.top();
		Path.pop();
		printf("-> (%d,%d) ", p.first, p.second);
		myfile << count << "\t0\t0\t16\t0\t0\t0\t0\t" << grid[p.first][p.second].lat << "\t" << grid[p.first][p.second].lon << "\t100.000000\t1" << endl;
		count++;
	}
	myfile.close();
	return;
}
*/
//Function to print every corner waypoints on the way
/*void tracePath(cell **cellDetails, Pair dest, int nodeNum, cell **grid, string pathFN)
{
	printf("\nThe Path is ");
	int row = dest.first;
	int col = dest.second;
	int count = 0;
	waypoint prev;
	string direction = "", prevDir = "";
	string N = "N", E = "E", S = "S", W = "W", NW = "NW", NE = "NE", SE = "SE", SW = "SW";
	stack<Pair> Path;
	ofstream myfile;

	myfile.open(pathFN, ofstream::app);
	myfile.precision(10);
	myfile << "QGC WPL 500" << endl;
	prev.row = 0;
	prev.col = 0;
	while (!(cellDetails[row][col].parent_i == row
		&& cellDetails[row][col].parent_j == col))
	{
		Path.push(make_pair(row, col));
		int temp_row = cellDetails[row][col].parent_i;
		int temp_col = cellDetails[row][col].parent_j;
		row = temp_row;
		col = temp_col;
	}

	Path.push(make_pair(row, col));
	int i = 0;
	while (!Path.empty())
	{
		pair<int, int> p = Path.top();
		Path.pop();
		if (p.first == prev.row && p.second == (prev.col + 1))
			direction = "E";
		if (p.second == prev.col && p.first == (prev.row + 1))
			direction = "S";
		if (p.first == (prev.row - 1) && p.second == prev.col)
			direction = "N";
		if (p.first == prev.row && p.second == (prev.col - 1))
			direction = "W";
		if (p.first == (prev.row - 1) && p.second == (prev.col - 1))
			direction = "NW";
		if (p.first == (prev.row - 1) && p.second == (prev.col + 1))
			direction = "NE";
		if (p.first == (prev.row + 1) && p.second == (prev.col + 1))
			direction = "SE";
		if (p.first == (prev.row + 1) && p.second == (prev.col - 1))
			direction = "SW";

		if (prevDir != direction && prevDir != "") {
			myfile << count << "\t0\t0\t16\t0\t0\t0\t0\t" << grid[prev.row][prev.col].lat << "\t" << grid[prev.row][prev.col].lon << "\t100.000000\t1" << endl;
			count++;
			//myfile << count << "\t0\t0\t16\t0\t0\t0\t0\t" << grid[p.first][p.second].lat << "\t" << grid[p.first][p.second].lon << "\t100.000000\t1" << endl;
			//count++;
			printf("-> (%d,%d) ", prev.row, prev.col);
			//printf("-> (%d,%d) ", p.first, p.second);
		}
		else if (prevDir != direction && prevDir == "") {
			myfile << count << "\t0\t0\t16\t0\t0\t0\t0\t" << grid[p.first][p.second].lat << "\t" << grid[p.first][p.second].lon << "\t100.000000\t1" << endl;
			count++;
			printf("-> (%d,%d) ", p.first, p.second);
		}

		prev.row = p.first;
		prev.col = p.second;
		prevDir = direction;
		i++;
	}
	myfile.close();
	return;
}*/
//Function to print 15th corner waypoints on the way
void tracePath(cell **cellDetails, Pair dest, int nodeNum, cell **grid, string pathFN)
{
	printf("\nThe Path is ");
	int row = dest.first;
	int col = dest.second;
	int count = 0;
	waypoint prev;
	string direction = "", prevDir = "";
	string N = "N", E = "E", S = "S", W = "W", NW = "NW", NE = "NE", SE = "SE", SW = "SW";
	stack<Pair> Path;
	ofstream myfile;

	myfile.open(pathFN, ofstream::app);
	myfile.precision(10);
	
	while (!(cellDetails[row][col].parent_i == row
		&& cellDetails[row][col].parent_j == col))
	{
		Path.push(make_pair(row, col));
		int temp_row = cellDetails[row][col].parent_i;
		int temp_col = cellDetails[row][col].parent_j;
		row = temp_row;
		col = temp_col;
	}

	Path.push(make_pair(row, col));
	int i = 0;
	while (!Path.empty())
	{
		pair<int, int> p = Path.top();
		Path.pop();
		//pathWP[count].lat = grid[p.first][p.second].lat;
		//pathWP[count].lon = grid[p.first][p.second].lon;
		if (i % 15 == 0) {
			printf("-> (%d,%d) ", p.first, p.second);
			//myfile << count << "\t0\t0\t16\t0\t0\t0\t0\t" << grid[p.first][p.second].lat << "\t" << grid[p.first][p.second].lon << "\t100.000000\t1" << endl;
			myfile << grid[p.first][p.second].id << "\t" << grid[p.first][p.second].FirstCol << "\t" << grid[p.first][p.second].SecondCol << "\t" <<
				grid[p.first][p.second].ThirdCol << "\t" << grid[p.first][p.second].ForthCol << "\t" << grid[p.first][p.second].FiveCol << "\t" <<
				grid[p.first][p.second].SixCol << "\t" << grid[p.first][p.second].SevenCol << "\t" << grid[p.first][p.second].lat << "\t" <<
				grid[p.first][p.second].lon << "\t" << grid[p.first][p.second].altitude << "\t" << grid[p.first][p.second].LastCol << endl;
		}
		count++;
		i++;
	}
	myfile.close();
	return;
}

// A Function to find the shortest path between 
// a given source cell to a destination cell according 
// to A* Search Algorithm 
void aStarSearch(cell **grid, Pair src, Pair dest, int nodeNum, string pathFN, bool &print)
{
	// If the source is out of range 
	if (isValid(src.first, src.second, nodeNum) == false)
	{
		printf("Source is invalid\n");
		return;
	}

	// If the destination is out of range 
	if (isValid(dest.first, dest.second, nodeNum) == false)
	{
		printf("Destination is invalid\n");
		return;
	}

	// Either the source or the destination is blocked 
	if (isUnBlocked(grid, src.first, src.second) == false ||
		isUnBlocked(grid, dest.first, dest.second) == false)
	{
		printf("Source or the destination is blocked\n");
		return;
	}

	// If the destination cell is the same as source cell 
	if (isDestination(src.first, src.second, dest) == true)
	{
		printf("We are already at the destination\n");
		return;
	}

	// Create a closed list and initialise it to false which means 
	// that no cell has been included yet 
	// This closed list is implemented as a boolean 2D array 
	bool **closedList = new bool*[nodeNum];
	for (int i = 0; i < nodeNum; i++) {
		closedList[i] = new bool[nodeNum];
	}
	for (int i = 0; i < nodeNum; i++) {
		for (int j = 0; j < nodeNum; j++) {
			closedList[i][j] = false;
		}
	}
	cout << "ClosedList[0][0]: " << closedList[0][0] << endl;
	// Declare a 2D array of structure to hold the details 
	//of that cell 
	cell **cellDetails = new cell*[nodeNum];
	for (int i = 0; i < nodeNum; i++) {
		cellDetails[i] = new cell[nodeNum];
	}
	//cell cellDetails[ROW][COL];

	int i, j;

	for (i = 0; i < nodeNum; i++)
	{
		for (j = 0; j < nodeNum; j++)
		{
			cellDetails[i][j].f = FLT_MAX;
			cellDetails[i][j].g = FLT_MAX;
			cellDetails[i][j].h = FLT_MAX;
			cellDetails[i][j].parent_i = -1;
			cellDetails[i][j].parent_j = -1;
		}
	}

	// Initialising the parameters of the starting node
	i = src.first, j = src.second;
	cellDetails[i][j].f = 0.0;
	cellDetails[i][j].g = 0.0;
	cellDetails[i][j].h = 0.0;
	cellDetails[i][j].parent_i = i;
	cellDetails[i][j].parent_j = j;

	/*
	Create an open list having information as-
	<f, <i, j>>
	where f = g + h,
	and i, j are the row and column index of that cell
	Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
	This open list is implenented as a set of pair of pair.*/
	set<pPair> openList;

	// Put the starting cell on the open list and set its 
	// 'f' as 0 
	openList.insert(make_pair(0.0, make_pair(i, j)));

	// We set this boolean value as false as initially 
	// the destination is not reached. 
	bool foundDest = false;
	bool algorithNeeded = false;		//do we need AStar algorithm if there are no blocked cells
	/*double heading = bear;
	double delta = 46;
	double headDiff = 10000;
	double tempHead = bear;
	*/
	while (!openList.empty())
	{
		pPair p = *openList.begin();
		if (algorithNeeded == true)
			print = true;
		// Remove this vertex from the open list 
		openList.erase(openList.begin());

		// Add this vertex to the closed list 
		i = p.second.first;
		j = p.second.second;
		closedList[i][j] = true;
		if (isUnBlocked(grid, i, j) == false)
			algorithNeeded = true;
		/*
			Generating all the 8 successor of this cell

				N.W N N.E
				\ | /
				\ | /
				W----Cell----E
					/ | \
				/ | \
				S.W S S.E

			Cell-->Popped Cell (i, j)
			N --> North	 (i-1, j)
			S --> South	 (i+1, j)
			E --> East	 (i, j+1)
			W --> West		 (i, j-1)
			N.E--> North-East (i-1, j+1)
			N.W--> North-West (i-1, j-1)
			S.E--> South-East (i+1, j+1)
			S.W--> South-West (i+1, j-1)*/

			// To store the 'g', 'h' and 'f' of the 8 successors 
		double gNew, hNew, fNew;

		//----------- 1st Successor (North) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i - 1, j, nodeNum) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i - 1, j, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i - 1][j].parent_i = i;
				cellDetails[i - 1][j].parent_j = j;
				printf("The destination cell is found\n");
				if (algorithNeeded == true)
					tracePath(cellDetails, dest, nodeNum, grid, pathFN);
				foundDest = true;
				return;
			}
			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i - 1][j] == false &&
				isUnBlocked(grid, i - 1, j) == true)
			{
				gNew = cellDetails[i][j].g + 1.0;
				hNew = calculateHValue(i - 1, j, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i - 1][j].f == FLT_MAX ||
					cellDetails[i - 1][j].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i - 1, j)));

					// Update the details of this cell 
					cellDetails[i - 1][j].f = fNew;
					cellDetails[i - 1][j].g = gNew;
					cellDetails[i - 1][j].h = hNew;
					cellDetails[i - 1][j].parent_i = i;
					cellDetails[i - 1][j].parent_j = j;
				}
			}
			if (isUnBlocked(grid, i - 1, j) == false)
				algorithNeeded = true;
		}

		//----------- 2nd Successor (South) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i + 1, j, nodeNum) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i + 1, j, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i + 1][j].parent_i = i;
				cellDetails[i + 1][j].parent_j = j;
				printf("The destination cell is found\n");
				if (algorithNeeded == true)
					tracePath(cellDetails, dest, nodeNum, grid, pathFN);
				foundDest = true;
				return;
			}
			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i + 1][j] == false &&
				isUnBlocked(grid, i + 1, j) == true)
			{
				gNew = cellDetails[i][j].g + 1.0;
				hNew = calculateHValue(i + 1, j, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i + 1][j].f == FLT_MAX ||
					cellDetails[i + 1][j].f > fNew)
				{
					openList.insert(make_pair(fNew, make_pair(i + 1, j)));
					// Update the details of this cell 
					cellDetails[i + 1][j].f = fNew;
					cellDetails[i + 1][j].g = gNew;
					cellDetails[i + 1][j].h = hNew;
					cellDetails[i + 1][j].parent_i = i;
					cellDetails[i + 1][j].parent_j = j;
				}
			}
			if (isUnBlocked(grid, i + 1, j) == false)
				algorithNeeded = true;
		}

		//----------- 3rd Successor (East) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i, j + 1, nodeNum) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i, j + 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i][j + 1].parent_i = i;
				cellDetails[i][j + 1].parent_j = j;
				printf("The destination cell is found\n");
				if (algorithNeeded == true)
					tracePath(cellDetails, dest, nodeNum, grid, pathFN);
				foundDest = true;
				return;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i][j + 1] == false &&
				isUnBlocked(grid, i, j + 1) == true)
			{
				gNew = cellDetails[i][j].g + 1.0;
				hNew = calculateHValue(i, j + 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i][j + 1].f == FLT_MAX ||
					cellDetails[i][j + 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i, j + 1)));

					// Update the details of this cell 
					cellDetails[i][j + 1].f = fNew;
					cellDetails[i][j + 1].g = gNew;
					cellDetails[i][j + 1].h = hNew;
					cellDetails[i][j + 1].parent_i = i;
					cellDetails[i][j + 1].parent_j = j;
				}
			}
			if (isUnBlocked(grid, i, j + 1) == false)
				algorithNeeded = true;
		}

		//----------- 4th Successor (West) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i, j - 1, nodeNum) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i, j - 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i][j - 1].parent_i = i;
				cellDetails[i][j - 1].parent_j = j;
				printf("The destination cell is found\n");
				if (algorithNeeded == true)
					tracePath(cellDetails, dest, nodeNum, grid, pathFN);
				foundDest = true;
				return;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i][j - 1] == false &&
				isUnBlocked(grid, i, j - 1) == true)
			{
				gNew = cellDetails[i][j].g + 1.0;
				hNew = calculateHValue(i, j - 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i][j - 1].f == FLT_MAX ||
					cellDetails[i][j - 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i, j - 1)));

					// Update the details of this cell 
					cellDetails[i][j - 1].f = fNew;
					cellDetails[i][j - 1].g = gNew;
					cellDetails[i][j - 1].h = hNew;
					cellDetails[i][j - 1].parent_i = i;
					cellDetails[i][j - 1].parent_j = j;
				}
			}
			if (isUnBlocked(grid, i, j - 1) == false)
				algorithNeeded = true;
		}

		//----------- 5th Successor (North-East) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i - 1, j + 1, nodeNum) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i - 1, j + 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i - 1][j + 1].parent_i = i;
				cellDetails[i - 1][j + 1].parent_j = j;
				printf("The destination cell is found\n");
				if (algorithNeeded == true)
					tracePath(cellDetails, dest, nodeNum, grid, pathFN);
				foundDest = true;
				return;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i - 1][j + 1] == false &&
				isUnBlocked(grid, i - 1, j + 1) == true)
			{
				gNew = cellDetails[i][j].g + 1.414;
				hNew = calculateHValue(i - 1, j + 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i - 1][j + 1].f == FLT_MAX ||
					cellDetails[i - 1][j + 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i - 1, j + 1)));

					// Update the details of this cell 
					cellDetails[i - 1][j + 1].f = fNew;
					cellDetails[i - 1][j + 1].g = gNew;
					cellDetails[i - 1][j + 1].h = hNew;
					cellDetails[i - 1][j + 1].parent_i = i;
					cellDetails[i - 1][j + 1].parent_j = j;
				}
			}
			if (isUnBlocked(grid, i - 1, j + 1) == false)
				algorithNeeded = true;
		}

		//----------- 6th Successor (North-West) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i - 1, j - 1, nodeNum) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i - 1, j - 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i - 1][j - 1].parent_i = i;
				cellDetails[i - 1][j - 1].parent_j = j;
				printf("The destination cell is found\n");
				if (algorithNeeded == true)
					tracePath(cellDetails, dest, nodeNum, grid, pathFN);
				foundDest = true;
				return;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i - 1][j - 1] == false &&
				isUnBlocked(grid, i - 1, j - 1) == true)
			{
				gNew = cellDetails[i][j].g + 1.414;
				hNew = calculateHValue(i - 1, j - 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i - 1][j - 1].f == FLT_MAX ||
					cellDetails[i - 1][j - 1].f > fNew)
				{
					openList.insert(make_pair(fNew, make_pair(i - 1, j - 1)));
					// Update the details of this cell 
					cellDetails[i - 1][j - 1].f = fNew;
					cellDetails[i - 1][j - 1].g = gNew;
					cellDetails[i - 1][j - 1].h = hNew;
					cellDetails[i - 1][j - 1].parent_i = i;
					cellDetails[i - 1][j - 1].parent_j = j;
				}
			}
			if (isUnBlocked(grid, i - 1, j - 1) == false)
				algorithNeeded = true;
		}

		//----------- 7th Successor (South-East) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i + 1, j + 1, nodeNum) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i + 1, j + 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i + 1][j + 1].parent_i = i;
				cellDetails[i + 1][j + 1].parent_j = j;
				printf("The destination cell is found\n");
				if (algorithNeeded == true)
					tracePath(cellDetails, dest, nodeNum, grid, pathFN);
				foundDest = true;
				return;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i + 1][j + 1] == false &&
				isUnBlocked(grid, i + 1, j + 1) == true)
			{
				gNew = cellDetails[i][j].g + 1.414;
				hNew = calculateHValue(i + 1, j + 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i + 1][j + 1].f == FLT_MAX ||
					cellDetails[i + 1][j + 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i + 1, j + 1)));

					// Update the details of this cell 
					cellDetails[i + 1][j + 1].f = fNew;
					cellDetails[i + 1][j + 1].g = gNew;
					cellDetails[i + 1][j + 1].h = hNew;
					cellDetails[i + 1][j + 1].parent_i = i;
					cellDetails[i + 1][j + 1].parent_j = j;
				}
			}
			if (isUnBlocked(grid, i + 1, j + 1) == false)
				algorithNeeded = true;
		}

		//----------- 8th Successor (South-West) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i + 1, j - 1, nodeNum) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i + 1, j - 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i + 1][j - 1].parent_i = i;
				cellDetails[i + 1][j - 1].parent_j = j;
				printf("The destination cell is found\n");
				if (algorithNeeded == true)
					tracePath(cellDetails, dest, nodeNum, grid, pathFN);
				foundDest = true;
				return;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i + 1][j - 1] == false &&
				isUnBlocked(grid, i + 1, j - 1) == true)
			{
				gNew = cellDetails[i][j].g + 1.414;
				hNew = calculateHValue(i + 1, j - 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//			 OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i + 1][j - 1].f == FLT_MAX ||
					cellDetails[i + 1][j - 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i + 1, j - 1)));

					// Update the details of this cell 
					cellDetails[i + 1][j - 1].f = fNew;
					cellDetails[i + 1][j - 1].g = gNew;
					cellDetails[i + 1][j - 1].h = hNew;
					cellDetails[i + 1][j - 1].parent_i = i;
					cellDetails[i + 1][j - 1].parent_j = j;
				}
			}
			if (isUnBlocked(grid, i + 1, j - 1) == false)
				algorithNeeded = true;
		}
	}
	// When the destination cell is not found and the open 
	// list is empty, then we conclude that we failed to 
	// reach the destiantion cell. This may happen when the 
	// there is no way to destination cell (due to blockages) 
	if (foundDest == false)
		printf("Failed to find the Destination Cell\n");
	delete[] * cellDetails;
	delete[] * closedList;
	return;
}

void findQuickestWay(cell **grid, waypoint *waypoints, int WPNum, int nodeNum, string pathFN) { 
	bool print = false;
	double heading = 0;
	for (int i = 0; i < WPNum - 1; i++) {
		int temp = i;
		print = false;
		heading = 0;
		if (i == 0) {
			Pair src = make_pair(waypoints[temp].row, waypoints[temp].col);
			Pair wpFirst = make_pair(waypoints[temp + 1].row, waypoints[temp + 1].col);
			heading = bearing(waypoints[temp].lat, waypoints[temp].lon, waypoints[temp + 1].lat, waypoints[temp + 1].lon);
			cout << "Targets are: " << src.first << " " << src.second << " -> " << wpFirst.first << " " << wpFirst.second << endl;
			//aStarSearch(grid, src, wpFirst, nodeNum, pathFN, print, heading);
			aStarSearch(grid, src, wpFirst, nodeNum, pathFN, print);
			if (print == false) {
				writeIndividualWPToFile(waypoints, i, pathFN, WPNum);
			}
		}
		else if (i == WPNum - 2) {
			Pair wpOneBeforeLast = make_pair(waypoints[temp].row, waypoints[temp].col);
			Pair wpLast = make_pair(waypoints[temp + 1].row, waypoints[temp + 1].col);
			heading = bearing(waypoints[temp - 1].lat, waypoints[temp - 1].lon, waypoints[temp].lat, waypoints[temp].lon);
			cout << "Targets are: " << wpOneBeforeLast.first << " " << wpOneBeforeLast.second << " -> " << wpLast.first << " " << wpLast.second << endl;
			//aStarSearch(grid, src, wpFirst, nodeNum, pathFN, print, heading);
			aStarSearch(grid, wpOneBeforeLast, wpLast, nodeNum, pathFN, print);
			if (print == false)
				writeIndividualWPToFile(waypoints, i, pathFN, WPNum);
		}
		else {
			Pair wpPrev = make_pair(waypoints[temp].row, waypoints[temp].col);
			Pair wpNext = make_pair(waypoints[temp + 1].row, waypoints[temp + 1].col);
			heading = bearing(waypoints[temp - 1].lat, waypoints[temp - 1].lon, waypoints[temp].lat, waypoints[temp].lon);
			cout << "Targets are: " << wpPrev.first << " " << wpPrev.second << " -> " << wpNext.first << " " << wpNext.second << endl;
			//aStarSearch(grid, wpPrev, wpNext, nodeNum, pathFN, print, heading);
			aStarSearch(grid, wpPrev, wpNext, nodeNum, pathFN, print);
			if (print == false)
				writeIndividualWPToFile(waypoints, i, pathFN, WPNum);
		}
	}
}

void bezierPath(double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3)
{
	double xu = 0.0, yu = 0.0, u = 0.0;
	int i = 0;
	for (u = 0.0; u <= 1.0; u += 0.1)
	{
		xu = pow(1 - u, 3)*x0 + 3 * u*pow(1 - u, 2)*x1 + 3 * pow(u, 2)*(1 - u)*x2
			+ pow(u, 3)*x3;
		yu = pow(1 - u, 3)*y0 + 3 * u*pow(1 - u, 2)*y1 + 3 * pow(u, 2)*(1 - u)*y2
			+ pow(u, 3)*y3;
		//cout << "Izracunao " << endl;
		writeBezier(xu, yu, "bezier.waypoints");
		//cout << "Zapisao " << endl;
	}
}

void bezier2(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
	double c1 = distanceEarth(x0, y0, x1, y1);
	double c2 = distanceEarth(x1, y1, x2, y2);
	double c3 = distanceEarth(x3, y3, x2, y2);
	double t1 = c1 / (c1 + c2 + c3);
	double t2 = (c1 + c2) / (c1 + c2 + c3);
	double a = t1 * (1 - t1)*(1 - t1) * 3;
	double b = (1 - t1)*t1*t1 * 3;
	double d = t2 * (1 - t2)*(1 - t2) * 3;
	double e = (1 - t2)*t2*t2 * 3;
	double c = x1 - (x0*pow(1 - t1, 3.0)) - (x3*pow(t1, 3));
	double f = x2 - (x0*pow(1 - t2, 3.0)) - (x3*pow(t2, 3));
	double g = y1 - (y0*pow(1 - t1, 3.0)) - (y3*pow(t1, 3));
	double h = y2 - (y0*pow(1 - t2, 3.0)) - (y3*pow(t2, 3));
	x2 = (c - a / d * f) / (b - a * e / d);
	x1 = (c - (b * x2)) / a;
	y2 = (g - a / d * h) / (b - a * e / d);
	y1 = (g - (b * y2)) / a;
	bezierPath(x0, x1, x2, x3, y0, y1, y2, y3);
}

int main() {
	/* Description of the Grid-
	1--> The cell is not blocked
	0--> The cell is blocked */
	int nodeNum = 0, WPNum = 0, boundNum = 0, obNum = 0, obPrev = 0, pathNum = 0;
	double feet = 0, miles = 0, maxDist = 0, distance = 0;
	bool newObstacle = false;
	string wpInFN = "files/waypointsInput.waypoints", boundInFN = "files/boundary.waypoints";
	string obInFN = "files/obstacles.waypoints", pathOutFN = "files/path.waypoints";
	string wpOutFN = "files/waypointOutput.waypoints", obOutFN = "files/obstaclesOutput.waypoints";
	string bezierOutFN = "bezier.waypoints";
	cell middleGround, startNode, middleGround2;
	cell *boundaryPoint, boundaryPoint1, boundaryPoint2;
	waypoint *waypoints, *wp1, *wp2;
	obstacle *obstacles;
	waypoint *pathWP;
	//path *path1;

	resetFiles(pathOutFN);
	resetFiles(wpOutFN);
	resetFiles(obOutFN);
	//resetFiles(bezierOutFN);

	countElementsInFile(wpInFN, WPNum);
	countElementsInFile(boundInFN, boundNum);
	countElementsInFile(obInFN, obNum);
	//countElementsInFile(pathOutFN, pathNum);

	//path1 = new path[nodeNum];
	obstacles = new obstacle[obNum];
	boundaryPoint = new cell[boundNum];
	waypoints = new waypoint[WPNum];
	//pathWP = new waypoint[pathNum];

	readBoundaryFile(boundaryPoint, boundInFN);
	readWaypointsFile(waypoints, wpInFN);
	readObstaclesFile(obstacles, obInFN);
	//readPathFile(pathWP, pathOutFN, path1, pathNum);

	maxDist = calculateMaxDistance(boundaryPoint, boundNum, boundaryPoint1, boundaryPoint2);
	feet = mToFeet(maxDist);
	miles = feet * 0.00018939393939393939393939;
	nodeNum = (int)round((feet * 2) / 10) + 1;

	middleGround = fillGPSCoordinate(boundaryPoint1, miles, 0);
	middleGround2 = fillGPSCoordinate(middleGround, miles, 270);
	startNode = middleGround2;

	//initialising grid with coordinates
	cell **grid = new cell*[nodeNum];
	for (int i = 0; i < nodeNum; i++) {
		grid[i] = new cell[nodeNum];
	}
	grid[0][0] = startNode;
	for (int i = 0; i < nodeNum; i++) {
		if (i != 0)
			grid[i][0] = fillGPSCoordinate(startNode, i * 0.000189393939 * 10, 180);
		for (int j = 0; j < nodeNum; j++) {
			grid[i][j] = fillGPSCoordinate(grid[i][0], j * 0.000189393939 * 10, 90);
		}
	}
	writeGPSGrid(grid, nodeNum, "grid.waypoints");
	//initialising waypoints
	for (int i = 0; i < WPNum; i++) {
		addWaypoint(grid, nodeNum, waypoints, i);
	}
	
	writeWPToFile(waypoints, WPNum, wpOutFN);
	//printWaypoints(waypoints, WPNum);

	//reading obstacle's file and performing if necessary obstacle avoidance
	cout << "\n\nWaiting for new entry in obstacles file...\n" << endl;
	while (1) {
		countElementsInFile(obInFN, obNum);
		newObstacle = false;
		if (obNum > 0 && obNum != obPrev) {
			obstacles = new obstacle[obNum];
			readObstaclesFile(obstacles, obInFN);
			resetFiles(pathOutFN);
			newObstacle = true;
			//initialising obstacle
			for (int i = 0; i < obNum; i++) {
				addObstacle(waypoints, WPNum, nodeNum, obstacles, i, grid, obNum);
			}

			//Calculating fastest way
			findQuickestWay(grid, waypoints, WPNum, nodeNum, pathOutFN);

			writeOBToFile(grid, nodeNum, obOutFN);
			obPrev = obNum;
			
			delete[] obstacles;
		}
		if (newObstacle)
			cout << "\n\nWaiting for new entry in obstacles file...\n" << endl;
	}
	//writeToFile(grid, nodeNum, "grid.waypoints");

	delete[] boundaryPoint, waypoints, wp1, wp2, pathWP;
	delete[] * grid;
	system("pause");
	return(0);
}
