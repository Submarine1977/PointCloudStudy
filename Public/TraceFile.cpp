#include "TraceFile.h"
#include <algorithm>

void project_to_segment(double x1, double y1, double x2, double y2, double x, double y, double *rx, double *ry)
{
	double xx, yy;
	double f = (x2 - x1) * (x - x1) + (y2 - y1) * (y - y1);       //dot product
	double l2 = (x2 - x1) *(x2 - x1) + (y2 - y1)*(y2 - y1);      //square of length
	double s = (x1 - x)*(y2 - y) - (y1 - y)*(x2 - x);            //area
	if (f < 0)
	{
		*rx = 0;
		*ry = sqrt((x - x1)* (x - x1) + (y - y1)*(y - y1));
	}
	else if (f > l2)
	{
		*rx = l2;
		*ry = sqrt((x - x2)* (x - x2) + (y - y2)*(y - y2));
	}
	else
	{
		xx = x1 + f / l2 * (x2 - x1);
		yy = y1 + f / l2 * (y2 - y1);
		*rx = sqrt((x1 - xx)* (x1 - xx) + (y1 - yy)*(y1 - yy));
		*ry = sqrt((x - xx)* (x - xx) + (y - yy)*(y - yy));
	}
	if (s < 0)
	{
		*ry = -*ry;
	}
}

int split_string(char* buffer, char** field)
{
	int m = 0;
	int i = 0;

	//skip space;
	while (buffer[i] == '\t' || buffer[i] == ' ')
	{
		i++;
	}
	if (buffer[i] == '\0')
	{
		return 0;
	}

	field[m] = buffer + i;
	m++;
	while (buffer[i] != '\0')
	{
		if (buffer[i] == ' ')
		{
			buffer[i] = '\0';
			field[m] = buffer + i + 1;
			m++;
			i++;
			//skip space
			while (buffer[i] == '\t' || buffer[i] == ' ')
			{
				i++;
			}
		}
		else
		{
			i++;
		}
	}
	return m;
}

CTraceFile::CTraceFile()
{
}


CTraceFile::~CTraceFile()
{
}

bool CTraceFile::Load(char *strTraceFile, int start_line, int line_num)
{
	FILE *f = NULL;
	char buffer[256];
	char *field[32];
	int n;
	errno_t err = fopen_s(&f, strTraceFile, "rt");
	if (f == NULL)
	{
		printf("Failed to open file %s, err = %d\n", strTraceFile, err);
		return false;
	}
	while (start_line > 0)
	{
		if (fgets(buffer, sizeof(buffer), f) == NULL)
		{
			printf("Failed to read file %s\n", strTraceFile);
			fclose(f);
			return false;
		}
		else
		{
			start_line--;
		}
	}
	while (line_num > 0)
	{
		if (fgets(buffer, sizeof(buffer), f) == NULL)
		{
			printf("Failed to read file %s\n", strTraceFile);
			fclose(f);
			return false;
		}
		else
		{
			line_num--;
			n = split_string(buffer, field);
			if (n != 14)
			{
				printf("Wrong Line:\n%s\n", buffer);
			}
			else
			{
				struct trace_point pt;
				pt.gps_time = atof(field[1]);
				pt.x = atof(field[3]);
				pt.y = atof(field[2]);
				pt.z = atof(field[4]);
				pt.heading = atof(field[10]);
				m_vecTrace.push_back(pt);
			}
		}
	}
	fclose(f);
	return true;
}

bool CTraceFile::DPReduce(double err)
{
	unsigned int i;
	int s, e;
	vector<int> vecIndexsToKeep;
	vector<struct trace_point> vecTrace;

	if (m_vecTrace.size() < 3)
	{
		return true;	
	}
	s = 0;
	e = m_vecTrace.size() - 1;
	vecIndexsToKeep.push_back(s);
	vecIndexsToKeep.push_back(e);
	DPReduce(s, e, err, vecIndexsToKeep);
	sort(vecIndexsToKeep.begin(), vecIndexsToKeep.end());
	for (i = 0; i < vecIndexsToKeep.size(); i++)
	{
		vecTrace.push_back(m_vecTrace[vecIndexsToKeep[i]]);
	}
	m_vecTrace = vecTrace;
	return true;
}

bool  CTraceFile::DPReduce(int s, int e, double err, vector<int> &vecIndexsToKeep)
{
	double max_dist = 0, dist, xx, yy;
	int    index = -1;
	int    i;

	if (e - s < 1)
	{
		return true;
	}
	for (i = s + 1; i < e; i++)
	{
		project_to_segment(m_vecTrace[s].x, m_vecTrace[s].y, m_vecTrace[e].x, m_vecTrace[e].y, m_vecTrace[i].x, m_vecTrace[i].y, &xx, &yy);
		dist = abs(yy);
		if (dist > max_dist)
		{
			max_dist = dist;
			index = i;
		}
	}

	if (max_dist > err && index >= 0)
	{
		vecIndexsToKeep.push_back(index);
		DPReduce(s, index, err, vecIndexsToKeep);
		DPReduce(index, e, err, vecIndexsToKeep);
	}
	return true;
}

bool CTraceFile::TransferToTraceCoords(double *x, double *y, double *z)
{
	double tyy = DBL_MAX, ty, tyy1, ty1;
	double d = 0, d1, tx, tx1, tz;
	unsigned int i;
	for (i = 0; i < m_vecTrace.size() - 1; i++)
	{
		tyy1 = (*x - (m_vecTrace[i].x + m_vecTrace[i + 1].x) / 2) * (*x - (m_vecTrace[i].x + m_vecTrace[i + 1].x) / 2)
			+ (*y - (m_vecTrace[i].y + m_vecTrace[i + 1].y) / 2) * (*y - (m_vecTrace[i].y + m_vecTrace[i + 1].y) / 2);
		d1 = sqrt((m_vecTrace[i].x - m_vecTrace[i + 1].x) * (m_vecTrace[i].x - m_vecTrace[i + 1].x) +
			(m_vecTrace[i].y - m_vecTrace[i + 1].y) * (m_vecTrace[i].y - m_vecTrace[i + 1].y));
		if (tyy1 < tyy)
		{
			project_to_segment(m_vecTrace[i].x, m_vecTrace[i].y, m_vecTrace[i+1].x, m_vecTrace[i + 1].y, *x, *y, &tx1, &ty1);
			if (ty1* ty1 < tyy)
			{
				tx = d + tx1;
				tz = *z - (m_vecTrace[i].z + (m_vecTrace[i + 1].z - m_vecTrace[i].z) * tx1 / d1);
				ty = ty1;
				tyy = ty1 * ty1;
			}
		}
		d += d1;
	}
	*x = tx;
	*y = ty;
	*z = tz;
	return true;
}