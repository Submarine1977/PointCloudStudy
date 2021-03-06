// CreatePCImage.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "..\Public\LASFile.h"
#include <algorithm>
#include <Windows.h>

#define GET_ALIGN(x) (((x+3)/4)*4)

struct statistic_data
{
	BYTE  d;
	int   count;
	float intensity;
	float z_mean,z_variance;
};
bool write_to_bmp(char* file, struct statistic_data **data, int width, int height)
{
	BITMAPFILEHEADER bfh;
	BITMAPINFOHEADER bmh;
	RGBQUAD bmiColors[256];

	FILE *fp;
	int i, j;
	fopen_s(&fp, file, "wb");
	if (fp == NULL)
	{
		return false;
	}
	bfh.bfType = 0x4d42;
	bfh.bfSize = GET_ALIGN(width) *height + 1078;
	bfh.bfOffBits = 1078;
	bfh.bfReserved1 = 0;
	bfh.bfReserved2 = 0;

	bmh.biSize = 40;
	bmh.biWidth = width;
	bmh.biHeight = height;
	bmh.biPlanes = 1;
	bmh.biBitCount = 8;
	bmh.biCompression = 0;
	bmh.biSizeImage = GET_ALIGN(width)*height;
	bmh.biXPelsPerMeter = 0;
	bmh.biYPelsPerMeter = 0;
	bmh.biClrUsed = 0;
	bmh.biClrImportant = 0;

	bmiColors[0].rgbRed   = 128;
	bmiColors[0].rgbGreen = 128;
	bmiColors[0].rgbBlue = 255;
	for (i = 1; i<256; i++)
	{
		bmiColors[i].rgbRed = (BYTE)i;
		bmiColors[i].rgbGreen = (BYTE)i;
		bmiColors[i].rgbBlue = (BYTE)i;
	}
	fwrite(&bfh, sizeof(BITMAPFILEHEADER), 1, fp); 
	fwrite(&bmh, sizeof(BITMAPINFOHEADER), 1, fp); 
	fwrite(bmiColors, sizeof(RGBQUAD), 256, fp);   

	BYTE d;
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < GET_ALIGN(width); j++)
		{
			if (j < width)
			{
				d = (BYTE)data[j][i].d;
			}
			else
			{
				d = 0;
			}
			fwrite(&d, 1, 1, fp);
		}
	}

	fclose(fp);
	return true;
}

int main(int argc, char *argv[])
{
	CLASFile lf;
	CLASFile::las_header lh;

	lf.Load(argv[1]);
	lf.GetHeader(&lh);

	int i, j;
	double x, y, z;
	unsigned short intensity;
	unsigned char  buffer[128];
	double minx = 0, miny = -20, maxx = 400, maxy = 20;
	double resolution = 0.1;
	
	int width = (int)((maxx - minx) / resolution + 0.5);
	int height = (int)((maxy - miny) / resolution + 0.5);

	int row, column;
	float dmax, dmin;
	struct statistic_data **data = new struct statistic_data*[width];
	for (i = 0; i < width; i++)
	{
		data[i] = new struct statistic_data[height];
		memset(data[i], 0, sizeof(struct statistic_data)*height);
	}
	for(i = 0; i < lh.number_point_records; i++)
	{
		lf.GetPointInfo(i, buffer);
		x = *((int*)(buffer)) * lh.x_scale + lh.x_offset;
		y = *((int*)(buffer + 4)) * lh.y_scale + lh.y_offset;
		z = *((int*)(buffer + 8)) * lh.z_scale + lh.z_offset;
		intensity = *(short*)(buffer + 12);

		column = (x - minx) / (maxx - minx) * width;
		row = (y - miny) / (maxy - miny) * height;
		if (column >= 0 && column < width && row >= 0 && row < height)
		{
			if (z < 0 && z >-2)
			{
				data[column][row].intensity += intensity;
				data[column][row].count++;
				data[column][row].z_mean += z;
			}
		}
	}
	for (i = 0; i < width; i++)
	{
		for (j = 0; j < height; j++)
		{
			if (data[i][j].count > 0)
			{
				data[i][j].z_mean /= data[i][j].count;
			}
		}
	}
	for (i = 0; i < lh.number_point_records; i++)
	{
		lf.GetPointInfo(i, buffer);
		x = *((int*)(buffer)) * lh.x_scale + lh.x_offset;
		y = *((int*)(buffer + 4)) * lh.y_scale + lh.y_offset;
		z = *((int*)(buffer + 8)) * lh.z_scale + lh.z_offset;

		column = (x - minx) / (maxx - minx) * width;
		row = (y - miny) / (maxy - miny) * height;
		if (column >= 0 && column < width && row >= 0 && row < height)
		{
			if (z < 0 && z >-2)
			{
				if (data[column][row].count > 0)
				{
					data[column][row].z_variance += (z - data[column][row].z_mean) * (z - data[column][row].z_mean) / data[column][row].count;
				}
			}
		}
	}
	for (i = 0; i < width; i++)
	{
		for (j = 0; j < height; j++)
		{
			if (data[i][j].count > 0)
			{
				data[i][j].z_variance = sqrt(data[i][j].count);
			}
		}
	}

	dmax = FLT_MIN;
	dmin = FLT_MAX;
	for (i = 0; i < width; i++)
	{
		for (j = 0; j < height; j++)
		{
			if (data[i][j].count > 0)
			{
				dmax = max(dmax, data[i][j].count);
				dmin = min(dmin, data[i][j].count);
			}
		}

	}
	for (i = 0; i < width; i++)
	{
		for (j = 0; j < height; j++)
		{
			if (data[i][j].count > 0)
			{
				data[i][j].d = (data[i][j].z_variance - dmin) / (dmax - dmin) * 254 + 1;
			}
			else
			{
				data[i][j].d = 0;
			}
		}
	}

	write_to_bmp(argv[2], data, width, height);

	for (i = 0; i < height; i++)
	{
		delete []data[i];
	}
	delete []data;

    return 0;
}

