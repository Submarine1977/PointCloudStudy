// ResamplePC.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "..\Public\LASFile.h"
#include <vector>
#include <map>
#include <algorithm>
#include <windows.h> 

using namespace std;

struct point_grids
{
	int grid_x, grid_y, grid_z;
	int seq;
};

int cmp(struct point_grids a, struct point_grids b) 
{
	if (a.grid_x != b.grid_x) 
		return a.grid_x < b.grid_x;
	if (a.grid_y != b.grid_y)
		return a.grid_y < b.grid_y;
	if (a.grid_z != b.grid_z)
		return a.grid_z < b.grid_z;
	return 0;
}

int main(int argc, char* argv[])
{
	long start;
	long end;


	unsigned long i;
	unsigned long   nClasses;
	double dAccurateX, dAccurateY, dAccurateZ;
	CLASFile lf;
	CLASFile::las_header lh;

	if (argc != 7)
	{
		printf("Usage: %s [Src LAS File] [Dest LAS File] [Classes] [Resample Accurate X (cm)] [Resample Accurate Y (cm)] [Resample Accurate Z (cm)]\n", argv[0]);
		return 1;
	}
	if (_strnicmp(argv[3], "0x", 2) == 0)
		nClasses = strtoul(argv[3], NULL, 16);
	else
		nClasses = atoi(argv[3]);
	dAccurateX = atof(argv[4]) / 100;
	dAccurateY = atof(argv[5]) / 100;
	dAccurateZ = atof(argv[6]) / 100;

	start = GetTickCount();
	printf("Reading input ... \n");
	if (!lf.Load(argv[1]))
	{
		printf("Failed to load las file %s\n", argv[1]);
		return 1;
	}
	end = GetTickCount();
	printf("Reading input finished,  %d seconds cost \n", (end - start + 500) / 1000);
	lf.GetHeader(&lh);

	unsigned char buffer[64];
	double x, y, z;
	int d = 0;
	int key;

	map<int, vector<struct point_grids>> points;

	start = GetTickCount();

	for (i = 0; i < lh.number_point_records; i++)
	{
		lf.GetPointInfo(i, buffer);
		if (nClasses != -1)
		{
			if (((1 << buffer[15]) & nClasses) == 0)
				continue;
		}

		x = *((int*)(buffer)) * lh.x_scale + lh.x_offset;
		y = *((int*)(buffer + 4)) * lh.y_scale + lh.y_offset;
		z = *((int*)(buffer + 8)) * lh.z_scale + lh.z_offset;

		struct point_grids pg;
		pg.grid_x = (int)((x - lh.minx) / dAccurateX);
		pg.grid_y = (int)((y - lh.miny) / dAccurateY);
		pg.grid_z = (int)((z - lh.minz) / dAccurateZ);
		pg.seq = i;

		key = (pg.grid_x / 1000) + ((pg.grid_y / 1000) << 10) + ((pg.grid_z / 1000) << 20);
		points[key].push_back(pg);

		if ((i + 1) % 1000000 == 0)
		{
			end = GetTickCount();
			printf("%d points loaded, %d seconds passed\n", i + 1, (end - start + 500) / 1000);
		}
	}
	end = GetTickCount();
	printf("%d points loaded, %d seconds passed\n", i, (end - start + 500) / 1000);


	start = GetTickCount();
	printf("Resampleing ... \n");
	map<int, vector<struct point_grids>>::iterator it;
	vector<struct point_grids>::iterator cur, prev;
	for (it = points.begin(); it != points.end(); ++it)
	{
		sort(it->second.begin(), it->second.end(), cmp);
		cur = it->second.begin();
		prev = cur;
		cur++;

		memset(buffer, 0, sizeof(buffer));
		for (; cur != it->second.end(); cur++)
		{
			if (cur->grid_x == prev->grid_x && cur->grid_y == prev->grid_y && cur->grid_z == prev->grid_z)
			{
				lf.SetPointInfo(cur->seq, buffer);
				d++;
			}
			prev = cur;
		}
    }
	lh.number_point_records -= d;
	lf.SetHeader(lh);
	end = GetTickCount();
	printf("Resampleing finished %d points removed, %d seconds cost!\n", d, (end - start + 500) / 1000);

	start = GetTickCount();
	printf("Writing result ... \n");
	lf.SaveTo(argv[2]);
	end = GetTickCount();
	printf("Writing result  %d seconds cost \n", (end - start + 500) / 1000);

    return 0;
}

