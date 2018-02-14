// ColorPC.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "..\Public\LASFile.h"


#define PI 3.14159265

void transfer_coords(double center_x, double center_y, double center_z, double roll, double pitch, double heading, double *x, double *y, double *z)
{
	double xx, yy, zz;

	*x -= center_x;
	*y -= center_y;
	*z -= center_z;

	pitch = pitch * PI / 180;
	heading = heading * PI / 180;
	roll = roll * PI / 180;
	
	//pitch   -- rotate by x axis;
    //x′ = x
	//y′ = ycosθ−zsinθ
	//z′ = ysinθ + zcosθ
	yy = *y * cos(pitch) - *z * sin(pitch);
	zz = *y * sin(pitch) + *z * cos(pitch);
	*y = yy;
	*z = zz;

	//roll    -- rotate by y axis;
	//x′ = zsinθ + xcosθ
	//y′ = y
	//z′ = zcosθ−xsinθ
	xx = *z * sin(roll) + *x * cos(roll);
	zz = *z * cos(roll) - *x * sin(roll);
	*x = xx;
	*z = zz;

	//heading -- rotate by z axis;
	//x′ = xcosθ−ysinθ
	//y′ = xsinθ + ycosθ
	//z′ = z
	xx = *x*cos(heading) - *y *sin(heading);
	yy = *x*sin(heading) + *y *cos(heading);
	*x = xx;
	*y = yy;
};
void transfer_rowcol(double x, double y, double z, short *row, short *col)
{
	double lon, lat;
	lon = atan2(y, x) / PI * 180 - 90; // -270. 90
	while (lon < -180)
		lon += 360;
	lat = atan2(z, sqrt(x*x + y * y)) / PI * 180;
	*row = (int)((lat + 90) / 180 * 4096);
	*col = (int)((lon + 180) / 360 * 8192);
};

void transfer_rowcol2(double x, double y, double z, short *row, short *col)
{

};

void get_rgb(unsigned char* bmp_buffer, short row, short col, unsigned short *r, unsigned short *g, unsigned short *b)
{
	row = row;
	*r = bmp_buffer[row * 8192 * 3 + col * 3];
	*g = bmp_buffer[row * 8192 * 3 + col * 3 + 1];
	*b = bmp_buffer[row * 8192 * 3 + col * 3 + 2];
};


int color_point_cloud(CLASFile *lf, const char* bmpFile,double center_x, double center_y, double center_z, double roll, double pitch, double heading)
{
	unsigned long i;
	double x, y, z;
	short row, col;
	unsigned short r, g, b;

	CLASFile::las_header lh;
	lf->GetHeader(&lh);
	unsigned char buffer[64];

	FILE *bf;
	size_t s;
	fopen_s(&bf, bmpFile, "rb");
	unsigned  char *bmp_buffer = new unsigned char[4096 * 8192 * 3];
	fseek(bf, 54, SEEK_SET);
	s = fread(bmp_buffer, 1, 4096 * 8192 * 3, bf);
	fclose(bf);

	for (i = 0; i < lh.number_point_records; i++)
	{
		lf->GetPointInfo(i, buffer);
		x = *((int*)(buffer)) * lh.x_scale + lh.x_offset;
		y = *((int*)(buffer + 4)) * lh.y_scale + lh.y_offset;
		z = *((int*)(buffer + 8)) * lh.z_scale + lh.z_offset;
		
		transfer_coords(center_x, center_y, center_z, roll, pitch, heading, &x, &y, &z);
		
		if (   y > 4 && y < 24 
			//|| y < -10 && y > -30
			)
		{
			transfer_rowcol(x, y, z, &row, &col);
			get_rgb(bmp_buffer, row, col, &r, &g, &b);

			memcpy(buffer + 28, &r, 2);
			memcpy(buffer + 30, &g, 2);
			memcpy(buffer + 32, &b, 2);

			lf->SetPointInfo(i, buffer);
		}
	}
	delete[]bmp_buffer;
    return 0;
}

/*
4592      38 - 10141002170728123205030  20170728123205030  16343.03100  123205.030  518341.665  3377657.819  19.639    30.518947616   114.191098021  1.2175441     0.8574699     125.0879528
4593      38 - 10141002170728123205843  20170728123205843  16343.84300  123205.843  518353.917  3377649.134  19.697    30.518869091   114.191225516  1.1541861     0.8623973     125.0730876
4594      38 - 10141002170728123206660  20170728123206660  16344.66100  123206.660  518366.168  3377640.433  19.768    30.518790416   114.191353008  0.9805725     0.8422101     125.2558904
4595      38 - 10141002170728123207485  20170728123207485  16345.48500  123207.485  518378.370  3377631.694  19.837    30.518711404   114.191479986  1.0851449     0.8854849     125.3950756
4596      38 - 10141002170728123208309  20170728123208309  16346.31000  123208.309  518390.578  3377622.920  19.893    30.518632072   114.191607022  1.2450461     0.9543102     125.4136008
*/

int main()
{
	const char * lasFile = "D:\\PointCloud\\data\\yuanshi\\38-10141002170728-123148-00209.las";
	const char * lasOutputFile = "D:\\PointCloud\\data\\yuanshi\\38-10141002170728-123148-00209_Color.las";

	const char * bmpFile[] = {
		"D:\\PointCloud\\data\\yuanshi\\pano\\38-10141002170728123205030.BMP",
		"D:\\PointCloud\\data\\yuanshi\\pano\\38-10141002170728123205843.BMP",
		"D:\\PointCloud\\data\\yuanshi\\pano\\38-10141002170728123206660.BMP"
	};

	double center_x[] = { 518341.665,518353.917,518366.168 };
	double center_y[] = { 3377657.819,3377649.134,3377640.433 };
	double center_z[] = { 19.639,19.697,19.768 };
	double roll[] = { 1.2175441,1.1541861, 0.9805725 };
	double pitch[] = { 0.8623973,0.8623973,0.8422101 };
	double heading[] = { 125.0879528,125.0730876, 125.2558904 };

	int i;
	CLASFile lf;
	lf.Load(lasFile);
	for (i = 0; i < sizeof(bmpFile) / sizeof(const char *); i++)
	{
		color_point_cloud(&lf, bmpFile[i], center_x[i], center_y[i], center_z[i], roll[i], pitch[i], heading[i]);
	}
	lf.SaveTo(lasOutputFile);
}