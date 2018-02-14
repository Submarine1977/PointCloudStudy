#pragma once
#include <vector>
using namespace std;

class CTraceFile
{
public:
	CTraceFile();
	~CTraceFile();
public:
	struct trace_point
	{
		double gps_time;
		double x, y, z;
		double heading;
	};
protected:
	vector<struct trace_point> m_vecTrace;
	bool DPReduce(int s, int e, double err, vector<int> &vecIndexsToKeep);
public:
	bool Load(char *strTraceFile, int start_line, int line_num);
	bool DPReduce(double err);
public:
	bool TransferToTraceCoords(double *x, double *y, double *z);
};

