// CreateATPC.cpp: 定义控制台应用程序的入口点。
// create the along trace point cloud

#include "stdafx.h"
#include <algorithm>
#include "..\Public\LASFile.h"
#include "..\Public\TraceFile.h"

using namespace std;

int main(int argc, char* argv[])
{
	if (argc != 6)
	{
		printf("Usage %s Input_LASFile TraceFile StartLine NumLines Output_LASFile\n", argv[0]);
		return -1;
	}

	CLASFile lf;
	lf.Load(argv[1]);

	CTraceFile tf;
	tf.Load(argv[2], atoi(argv[3]), atoi(argv[4]));
	tf.DPReduce(1);

	unsigned char buffer[64];
	unsigned int n;

	CLASFile::las_header lh,newlh;
	lf.GetHeader(&lh);
	lf.GetHeader(&newlh);

	newlh.x_offset = newlh.y_offset = newlh.z_offset = 0;

	double x, y, z;
	for (n = 0; n < lh.number_point_records; n++)
	{
		lf.GetPointInfo(n, buffer);
		x = *((int*)(buffer)) * lh.x_scale + lh.x_offset;
		y = *((int*)(buffer + 4)) * lh.y_scale + lh.y_offset;
		z = *((int*)(buffer + 8)) * lh.z_scale + lh.z_offset;
		tf.TransferToTraceCoords(&x, &y, &z);
		*((int*)(buffer)) = (int)((x - newlh.x_offset) / lh.x_scale + 0.5);
		*((int*)(buffer + 4)) = (int)((y - newlh.y_offset) / lh.y_scale + 0.5);
		*((int*)(buffer + 8)) = (int)((z - newlh.z_offset) / lh.z_scale + 0.5);
		lf.SetPointInfo(n, buffer);
		if (n == 0)
		{
			newlh.minx = newlh.maxx = x;
			newlh.miny = newlh.maxy = y;
			newlh.minz = newlh.maxz = z;
		}
		else
		{
			newlh.minx = min(newlh.minx, x);
			newlh.maxx = max(newlh.maxx, x);
			newlh.miny = min(newlh.miny, y);
			newlh.maxy = max(newlh.maxy, y);
			newlh.minz = min(newlh.minz, z);
			newlh.maxz = max(newlh.maxz, z);
		}
	}
	lf.SetHeader(newlh);
	lf.SaveTo(argv[5]);
    return 0;
}

