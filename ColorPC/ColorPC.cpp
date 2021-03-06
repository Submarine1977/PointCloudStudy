// ColorPC.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "..\Public\LASFile.h"
#include "Matrix.h"
#include "windows.h"
#include "Z:\TechnicalResearch\cximage702_full\CxImage\ximage.h"
//#define PI 3.14159265

#include <ladybug.h>
#include <ladybugrenderer.h>
#include <ladybuggeom.h>
#include <ladybugstream.h>
#include <ladybuggps.h>
#include <ladybugvideo.h>

#define _CHECK_ERROR \
    if( error != LADYBUG_OK ) \
{ \
    return error; \
} \

LadybugContext context;
LadybugStreamHeadInfo streamHeaderInfo;

//测试
//0.031200   1.457700   0.499600 - 0.174971 - 0.281440 - 0.540180
double imu2cam_cx = 0.031200, imu2cam_cy = 1.457700, imu2cam_cz = 0.499600;
double imu2cam_roll = -0.174971, imu2cam_pitch = -0.281440, imu2cam_heading = -0.540180;

//char pszInputStream[_MAX_PATH] = "\\\\192.168.3.110\\f\\DATA\\colour_cloud\\1304-1-002-170509\\13041002170509100827-000000.pgr";
//const char * imageDir          = "\\\\192.168.3.110\\f\\DATA\\colour_cloud\\1304-1-002-170509\\data\\img\\";
//const char * imagepostFile     = "\\\\192.168.3.110\\f\\DATA\\color_pointCloud_data_for_test\\1304-1-002-170509\\data\\imgpost\\40-13041002170509-imgpost.txt";
//const char * lasFile           = "\\\\192.168.3.110\\f\\DATA\\colour_cloud\\1304-1-002-170509\\data\\las\\40-13041002170509-155449-00349.las";
//const char * lasOutputFile     = "Z:\\40-13041002170509-155449-00349.las";



////机场高速
//double imu2cam_cx = 0.0142, imu2cam_cy = 1.4542, imu2cam_cz = 0.43385;
//double imu2cam_roll = 0.1, imu2cam_pitch = -0.5, imu2cam_heading = -0.05;
//
//char pszInputStream[_MAX_PATH] = "Z:\\TechnicalResearch\\PointCloud\\1001-1-004-180129-Pre\\1001-1-004-180129\\Photo\\10011004180129112955-000000.pgr";
//const char * imageDir          = "Z:\\TechnicalResearch\\PointCloud\\1001-1-004-180129\\data\\img\\";
//const char * imagepostFile     = "Z:\\TechnicalResearch\\PointCloud\\1001-1-004-180129\\data\\imgpost\\39-10011004180129-imgpost.txt";
//const char * lasFile           = "Z:\\TechnicalResearch\\PointCloud\\1001-1-004-180129\\data\\las\\39-10011004180129-122353-00151.las";
//const char * lasOutputFile     = "Z:\\TechnicalResearch\\PointCloud\\1001-1-004-180129\\data\\las\\39-10011004180129-122353-00151_Color.las";



//北清路
//-0.239852 0.075247 - 0.264791 
//- 0.031200 1.457700 0.499600
//double imu2cam_cx = -0.031200, imu2cam_cy = 1.457700, imu2cam_cz = 0.499600;
//double imu2cam_roll = -0.239852, imu2cam_pitch = 0.075247, imu2cam_heading = -0.264791;
//
//char pszInputStream[_MAX_PATH] = "Z:\\TechnicalResearch\\PointCloud\\1001-1-001-180515\\Photo\\10011001180515145748-000000.pgr";
//const char * imageDir          = "Z:\\TechnicalResearch\\PointCloud\\1001-1-001-180515\\1001-1-001-180515\\data\\img\\";
//const char * imagepostFile     = "Z:\\TechnicalResearch\\PointCloud\\1001-1-001-180515\\1001-1-001-180515\\data\\imgpost\\39-10011001180515-imgpost.txt";
//const char * lasFile           = "Z:\\TechnicalResearch\\PointCloud\\1001-1-001-180515\\1001-1-001-180515\\data\\las\\39-10011001180515-145829-00002.las";
//const char * lasOutputFile     = "Z:\\TechnicalResearch\\PointCloud\\1001-1-001-180515\\1001-1-001-180515\\data\\las\\39-10011001180515-145829-00002_Color.las";


//小路
//-0.239852 0.075247 - 0.264791 
//- 0.031200 1.457700 0.499600
//double imu2cam_cx = -0.031200, imu2cam_cy = 1.457700, imu2cam_cz = 0.499600;
//double imu2cam_roll = -0.239852, imu2cam_pitch = 0.075247, imu2cam_heading = -0.264791;
//
//char pszInputStream[_MAX_PATH] = "Z:\\TechnicalResearch\\PointCloud\\1001-1-001-180516\\Photo\\10011001180516104017-000000.pgr";
//const char * imageDir          = "Z:\\TechnicalResearch\\PointCloud\\1001-1-001-180516\\1001-1-001-180516\\data\\img\\";
//const char * imagepostFile     = "Z:\\TechnicalResearch\\PointCloud\\1001-1-001-180516\\1001-1-001-180516\\data\\imgpost\\39-10011001180516-imgpost.txt";
//const char * lasFile           = "Z:\\TechnicalResearch\\PointCloud\\1001-1-001-180516\\1001-1-001-180516\\data\\las\\39-10011001180516-104038-00001.las";
//const char * lasOutputFile     = "Z:\\TechnicalResearch\\PointCloud\\1001-1-001-180516\\1001-1-001-180516\\data\\las\\39-10011001180516-104038-00001_Color.las";



CMatrix RIMU2CAM(3, 3), RPC2IMU(3, 3), RIMU2LEN0(3,3);

void build_imu_len0_matrix()
{
	double rx = 1.566671, ry = 0.003142, rz = -1.566895;
	CMatrix R(3, 3), Rx(3, 3), Ry(3, 3), Rz(3, 3);

	Rx(0, 0) = 1.0; Rx(0, 1) = 0.0;   Rx(0, 2) = 0.0;
	Rx(1, 0) = 0.0; Rx(1, 1) = cos(rx);   Rx(1, 2) = -sin(rx);
	Rx(2, 0) = 0.0; Rx(2, 1) = sin(rx);   Rx(2, 2) = cos(rx);

	Ry(0, 0) = cos(ry); Ry(0, 1) = 0.0;   Ry(0, 2) = sin(ry);
	Ry(1, 0) = 0.0; Ry(1, 1) = 1.0;   Ry(1, 2) = 0.0;
	Ry(2, 0) = -sin(ry); Ry(2, 1) = 0.0;   Ry(2, 2) = cos(ry);

	Rz(0, 0) = cos(rz); Rz(0, 1) = -sin(rz);   Rz(0, 2) = 0.0;
	Rz(1, 0) = sin(rz); Rz(1, 1) = cos(rz);   Rz(1, 2) = 0.0;
	Rz(2, 0) = 0.0; Rz(2, 1) = 0.0;   Rz(2, 2) = 1.0;
	RIMU2LEN0 = Rz * Ry * Rx;
}

void build_imu_cam_maxtrix()
{
	imu2cam_pitch = imu2cam_pitch * PI / 180;
	imu2cam_heading = imu2cam_heading * PI / 180;
	imu2cam_roll = imu2cam_roll * PI / 180;

	CMatrix R(3, 3), Rx(3, 3), Ry(3, 3), Rz(3, 3);

	Rx(0, 0) = 1.0; Rx(0, 1) = 0.0;   Rx(0, 2) = 0.0;
	Rx(1, 0) = 0.0; Rx(1, 1) = cos(imu2cam_pitch);   Rx(1, 2) = -sin(imu2cam_pitch);
	Rx(2, 0) = 0.0; Rx(2, 1) = sin(imu2cam_pitch);   Rx(2, 2) = cos(imu2cam_pitch);

	Ry(0, 0) = cos(imu2cam_roll); Ry(0, 1) = 0.0;   Ry(0, 2) = sin(imu2cam_roll);
	Ry(1, 0) = 0.0; Ry(1, 1) = 1.0;   Ry(1, 2) = 0.0;
	Ry(2, 0) = -sin(imu2cam_roll); Ry(2, 1) = 0.0;   Ry(2, 2) = cos(imu2cam_roll);

	Rz(0, 0) = cos(imu2cam_heading); Rz(0, 1) = -sin(imu2cam_heading);   Rz(0, 2) = 0.0;
	Rz(1, 0) = sin(imu2cam_heading); Rz(1, 1) = cos(imu2cam_heading);   Rz(1, 2) = 0.0;
	Rz(2, 0) = 0.0; Rz(2, 1) = 0.0;   Rz(2, 2) = 1.0;

	RIMU2CAM = Ry * Rx * Rz;
}

void build_pc_imu_maxtrix(double roll, double pitch, double heading)
{
	double yaw;
	pitch = pitch *  PI / 180;
	yaw = -heading * PI / 180;
	roll = roll *    PI / 180;

	CMatrix R(3, 3);

	R(0, 0) = cos(yaw)*cos(roll) - sin(yaw)*sin(pitch)*sin(roll);
	R(0, 1) = -sin(yaw)*cos(pitch);
	R(0, 2) = cos(yaw)*sin(roll) + sin(yaw)*sin(pitch)*cos(roll);

	R(1, 0) = sin(yaw)*cos(roll) + cos(yaw)*sin(pitch)*sin(roll);
	R(1, 1) = cos(yaw)*cos(pitch);
	R(1, 2) = sin(yaw)*sin(roll) - cos(yaw)*sin(pitch)*cos(roll);

	R(2, 0) = -cos(pitch)*sin(roll);
	R(2, 1) = sin(pitch);
	R(2, 2) = cos(pitch)*cos(roll);

	RPC2IMU = R;
}

void transfer_coords_pc2imu(double center_x, double center_y, double center_z, double *x, double *y, double *z)
{
	double xx, yy, zz;
	//double yaw;

	*x -= center_x;
	*y -= center_y;
	*z -= center_z;

	//pitch = pitch * PI / 180;
	//yaw   = -heading * PI / 180;
	//roll  = roll * PI / 180;

	//CMatrix R(3, 3);

	//R(0, 0) = cos(yaw)*cos(roll) - sin(yaw)*sin(pitch)*sin(roll);
	//R(0, 1) = -sin(yaw)*cos(pitch);
	//R(0, 2) = cos(yaw)*sin(roll) + sin(yaw)*sin(pitch)*cos(roll);

	//R(1, 0) = sin(yaw)*cos(roll) + cos(yaw)*sin(pitch)*sin(roll);
	//R(1, 1) = cos(yaw)*cos(pitch);
	//R(1, 2) = sin(yaw)*sin(roll) - cos(yaw)*sin(pitch)*cos(roll);

	//R(2, 0) = -cos(pitch)*sin(roll);
	//R(2, 1) = sin(pitch);
	//R(2, 2) = cos(pitch)*cos(roll);

	xx = *x * RPC2IMU(0, 0) + *y *RPC2IMU(1, 0) + *z * RPC2IMU(2, 0);
	yy = *x * RPC2IMU(0, 1) + *y *RPC2IMU(1, 1) + *z * RPC2IMU(2, 1);
	zz = *x * RPC2IMU(0, 2) + *y *RPC2IMU(1, 2) + *z * RPC2IMU(2, 2);

	*x = xx;
	*y = yy;
	*z = zz;
}

void transfer_coords_imu2cam(double *x, double *y, double *z)
{
	double xx, yy, zz;

	*x -= imu2cam_cx;
	*y -= imu2cam_cy;
	*z -= imu2cam_cz;

	//pitch = pitch * PI / 180;
	//heading = heading * PI / 180;
	//roll  = roll * PI / 180;

	//CMatrix R(3, 3), Rx(3, 3), Ry(3, 3), Rz(3, 3);

	//Rx(0, 0) = 1.0; Rx(0, 1) = 0.0;   Rx(0, 2) = 0.0;
	//Rx(1, 0) = 0.0; Rx(1, 1) = cos(pitch);   Rx(1, 2) = -sin(pitch);
	//Rx(2, 0) = 0.0; Rx(2, 1) = sin(pitch);   Rx(2, 2) = cos(pitch);

	//Ry(0, 0) = cos(roll); Ry(0, 1) = 0.0;   Ry(0, 2) = sin(roll);
	//Ry(1, 0) = 0.0; Ry(1, 1) = 1.0;   Ry(1, 2) = 0.0;
	//Ry(2, 0) = -sin(roll); Ry(2, 1) = 0.0;   Ry(2, 2) = cos(roll);

	//Rz(0, 0) = cos(heading); Rz(0, 1) = -sin(heading);   Rz(0, 2) = 0.0;
	//Rz(1, 0) = sin(heading); Rz(1, 1) = cos(heading);   Rz(1, 2) = 0.0;
	//Rz(2, 0) = 0.0; Rz(2, 1) = 0.0;   Rz(2, 2) = 1.0;

	//R = Ry * Rx * Rz;


	xx = *x * RIMU2CAM(0, 0) + *y *RIMU2CAM(1, 0) + *z * RIMU2CAM(2, 0);
	yy = *x * RIMU2CAM(0, 1) + *y *RIMU2CAM(1, 1) + *z * RIMU2CAM(2, 1);
	zz = *x * RIMU2CAM(0, 2) + *y *RIMU2CAM(1, 2) + *z * RIMU2CAM(2, 2);

	*x = xx;
	*y = yy;
	*z = zz;
}

bool transfer_coords_imu2len0(double x, double y, double z, double *r, double *c, double *d2)
{
	double xx, yy, zz;
	double cx = -0.003752, cy = 1.539427, cz = 0.529968;
	double fl = 466.865424;
	double cx0 = 1225.930864;
	double cy0 = 1040.092896;
	
	x -= cx;
	y -= cy;
	z -= cz;

	xx = x * RIMU2LEN0(0, 0) + y * RIMU2LEN0(0, 1) + z * RIMU2LEN0(0, 2);
	yy = x * RIMU2LEN0(1, 0) + y * RIMU2LEN0(1, 1) + z * RIMU2LEN0(1, 2);
	zz = x * RIMU2LEN0(2, 0) + y * RIMU2LEN0(2, 1) + z * RIMU2LEN0(2, 2);

	*r = xx / zz * fl + cx0;
	*c = yy / zz * fl + cy0;
	*d2 = xx * xx + yy * yy + zz * zz;
	return zz > 0;
}

/*
void transfer_coords(double center_x, double center_y, double center_z, double roll, double pitch, double heading, double *x, double *y, double *z)
{
	double xx, yy, zz;

	*x -= center_x;
	*y -= center_y;
	*z -= center_z;

	pitch = pitch * PI / 180;
	heading = heading * PI / 180;
	roll = roll * PI / 180;

	//heading -- rotate by z axis;
	//x′ = xcosθ−ysinθ
	//y′ = xsinθ + ycosθ
	//z′ = z
	xx = *x*cos(heading) - *y *sin(heading);
	yy = *x*sin(heading) + *y *cos(heading);
	*x = xx;
	*y = yy;

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

};
*/
/*
void makeTransformation(const double rotX, const double rotY, const double rotZ,
	double matrix[3][3])
{
	double cosX, sinX, cosY, sinY, cosZ, sinZ;

	cosX = cos(rotX);		sinX = sin(rotX);
	cosY = cos(rotY);		sinY = sin(rotY);
	cosZ = cos(rotZ);		sinZ = sin(rotZ);

	// cz*cy;
	matrix[0][0] = cosZ * cosY;
	// cz*sy*sx - sz*cx;
	matrix[0][1] = cosZ * sinY * sinX - sinZ * cosX;
	// cz*sy*cx + sz*sx;
	matrix[0][2] = cosZ * sinY * cosX + sinZ * sinX;

	// sz*cy;
	matrix[1][0] = sinZ * cosY;
	// sz*sy*sx + cz*cx;
	matrix[1][1] = sinZ * sinY * sinX + cosZ * cosX;
	// sz*sy*cx - cz*sx;
	matrix[1][2] = sinZ * sinY * cosX - cosZ * sinX;

	//-sy;
	matrix[2][0] = -sinY;
	//cy*sx;
	matrix[2][1] = cosY * sinX;
	//cy*cx;
	matrix[2][2] = cosY * cosX;
}
*/
/*
void applyTransformation(const double matrix[3][3], double & x, double & y, double & z)
{
	double outPoint[3];
	for (int r = 0; r < 3; r++)
	{
		outPoint[r] = 
			+ matrix[0][r] * x
			+ matrix[1][r] * y
			+ matrix[2][r] * z;
	}
	x = outPoint[0];
	y = outPoint[1];
	z = outPoint[2];
}
*/
/*
void transfer_xyz2rowcol(double x, double y, double z, short *row, short *col)
{
	// Map the rectified coordinate to camera-local 3D coordinate (dLocalX, dLocalY, dLocalZ)
	//
	// Here, we solve the following equations:
	//
	//  (Pin-hole camera model)
	//  kRectifiedX = ( dLocalX / dLocalZ) * dFocalLen + dCameraCenterX
	//  kRectifiedY = ( dLocalY / dLocalZ) * dFocalLen + dCameraCenterY
	//
	//  (Set the constraint that the point is on a sphere with the radius of kSphereSize
	//  dLocalX^2 + dLocalY^2 + dLocalZ^2 = kSphereSize^2

	//#CamToLadybugEulerZYX Rx Ry Rz Tx Ty Tz
	//	CamToLadybugEulerZYX 2.388081 1.570075 2.388081 0.061129 - 0.000777 - 0.000016
    double dExtrinsics[] = { 2.388081, 1.570075, 2.388081, 0.061129, - 0.000777, - 0.000016 };
	
	//focalLength 0.187709 0.224371
	double dFocalLen = 0.187709;

	//Center  19.999997 - 0.000777 - 0.010511   0.508992   0.507099	;
	double dCameraCenterX = 19.999997;
	double dCameraCenterY = -0.000777;

	const double rotX = dExtrinsics[0];
	const double rotY = dExtrinsics[1];
	const double rotZ = dExtrinsics[2];
	const double transX = dExtrinsics[3];
	const double transY = dExtrinsics[4];
	const double transZ = dExtrinsics[5];

	double toGlobalCoords[3][3]; // Craig's Matrix
	x -= transX;
	y -= transY;
	z -= transZ;
	makeTransformation(rotX, rotY, rotZ, toGlobalCoords);
	applyTransformation(toGlobalCoords, x, y, z);

	double dSphereSize = sqrt(x * x + y * y + z * z);
	*col = (x / z) * dFocalLen + dCameraCenterX;
	*row = (y / z) * dFocalLen + dCameraCenterY;
}
*/
/*
void transfer_rowcol(double x, double y, double z, short *row, short *col)
{
	double lon, lat;
	lon = atan2(y, x) / PI * 180 - 90;
	while (lon < -180)
		lon += 360;
	lat = atan2(z, sqrt(x*x + y * y)) / PI * 180;
	*row = (int)((lat + 90) / 180 * 4096);
	*col = 8192 - (int)((lon + 180) / 360 * 8192);
};

void transfer_rowcol2(double x, double y, double z, short *row, short *col)
{

};
*/
void get_rgb(CxImage* img, unsigned char* mask_buffer, short width, short hight, short row, short col, unsigned short *r, unsigned short *g, unsigned short *b)
{
	row = row;
	if ((mask_buffer[row * 308 + col / 8] & (1 << (7 - (col % 8)))) == 0)
	{
		*b = *g = *r = 0;
		return;
	}
	RGBQUAD color = img->GetPixelColor(col, row);
	*b = color.rgbBlue;
	*g = color.rgbGreen;
	*r = color.rgbRed;
};


int color_point_cloud(CLASFile *lf, float dist2[], const char* bmpFile, unsigned char** mask_buffer, double center_x, double center_y, double center_z)
{
	int ncamera;
	unsigned long i;
	double x, y, z;
	short row, col;
	unsigned short r, g, b;

	CLASFile::las_header lh;
	lf->GetHeader(&lh);
	unsigned char buffer[64];

	CxImage img[5];
	char fileName[512];
	TCHAR name[512];
	for (i = 0; i < 5; i++)
	{
		sprintf_s(fileName, sizeof(fileName)/sizeof(char), ("%s-%d.jpg"), bmpFile, i);
		name[::MultiByteToWideChar(CP_ACP, 0, fileName, strlen(fileName), (LPWCH)name, sizeof(name)/sizeof(TCHAR))] = '\0';
		if (!img[i].Load(name))
		{
			printf("Failed to load image %s\n", fileName);
			return 0;
		}
	}

	double dist;
	double dr, dc;
	double angle;

	for (i = 0; i < lh.number_point_records; i++)
	{
		lf->GetPointInfo(i, buffer);
		x = *((int*)(buffer)) * lh.x_scale + lh.x_offset;
		y = *((int*)(buffer + 4)) * lh.y_scale + lh.y_offset;
		z = *((int*)(buffer + 8)) * lh.z_scale + lh.z_offset;

		if (fabs(x - center_x) > 50 || fabs(y - center_y) > 50 || fabs(z - center_z) > 50)
		{
			continue;
		}

		transfer_coords_pc2imu(center_x, center_y, center_z, &x, &y, &z); //image trace

		//x = 435212.741367;
		//y = 4437906.925457;
		//z = 34.563;
		//transfer_coords_pc2imu(center_x, center_y, center_z, &x, &y, &z); //image trace


		////transfer_coords_imu2len0(-1.5798,6.5432,-1.7731, &dr, &dc, &dist);
		//transfer_coords_imu2len0(x,y,z, &dr, &dc, &dist);
		////transfer_coords_imu2len0(6.6924,48.3193,-0.9704, &dr, &dc, &dist);


		ncamera = 0;
		if (!transfer_coords_imu2len0(x, y, z, &dr, &dc, &dist))
			continue;

		//transfer_coords_imu2cam(&x, &y, &z); //pp

		//angle = atan2(y, x)  * 180 / PI;
		//if (angle < 0)
		//	angle += 360;

		//if (angle > 54 && angle <= 126)
		//	ncamera = 0;
		//if ( angle > 0 && angle <= 54 || 
		//	 angle > 342 && angle <= 0)
		//	ncamera = 1;
		//if (angle > 270 && angle <= 342)
		//	ncamera = 2;
		//if (angle > 198 && angle <= 270)
		//	ncamera = 3;
		//if (angle > 126 && angle <= 199)
		//	ncamera = 4;
		//
		//if (ncamera != 0)
		//	continue;

		//if (y > 4 && y < 34
			//	|| y < -10 && y > -30
		//	)
			//if( y < 4 && y > -4)
		{
			//dist = x * x + y * y + z * z;
			if (dist < dist2[i])
			{
				//dist2[i] = dist;
				//transfer_xyz2rowcol(x, y, z, &row, &col);

				//LadybugError  error = ::ladybugXYZtoRC(context, y, -x, z, ncamera, &dr, &dc, NULL);

				if (dr > 0 && dc > 0)
				{
					row = (short)(dc + 0.5);
					col = (short)(dr + 0.5);
					row = 2048 - row - 1;
					//col = 2448 - col - 1;
					if (row >= 0 && row < 2048 && col >= 0 && col < 2448)
					{
						get_rgb(&img[ncamera], mask_buffer[ncamera], 2448, 2048, row, col, &r, &g, &b);
						if (r != 0 || g != 0 || b != 0)
						{
							dist2[i] = dist;
							memcpy(buffer + 28, &r, 2);
							memcpy(buffer + 30, &g, 2);
							memcpy(buffer + 32, &b, 2);

							lf->SetPointInfo(i, buffer);
						}
					}
				}
			}
		}
	}
    
	return 0;
}


bool get_las_start_end_time(CLASFile *lf, double *start, double *end)
{
	int i;
	CLASFile::las_header lh;
	lf->GetHeader(&lh);
	unsigned char buffer[64];
	double t = 0;
	for (i = 0; i < lh.number_point_records; i++)
	{
		lf->GetPointInfo(i, buffer);
		memcpy(&t, buffer + 20, 8);
		if (i == 0)
		{
			*start = *end = t;
		}
		else
		{
			*start = min(*start, t);
			*end = max(*end, t);
		}
	}
	return i > 0;
}

struct image_info
{
	char name[32];
	double center_x, center_y, center_z;
	double roll, pitch, heading;
};

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
		if (buffer[i] == '\t' || buffer[i] == ' ')
		{
			buffer[i] = '\0';
			i++;
			//skip space
			while (buffer[i] == '\t' || buffer[i] == ' ')
			{
				i++;
			}
			field[m] = buffer + i;
			m++;
		}
		else
		{
			i++;
		}
	}
	return m;
}

bool load_pp(const char *ppfile)
{
	FILE *f;
	char buffer[1024];
	char *field[64];
	double t;
	int i, m;
	fopen_s(&f, ppfile, "r");
	if (f == NULL)
	{
		return false;
	}
	bool first_line = true;
	
	while (fgets(buffer, 1024, f))
	{
		i = strlen(buffer);
		while (buffer[i - 1] == ' ' || buffer[i - 1] == '\t' || buffer[i - 1] == '\n' || buffer[i - 1] == '\r')
		{
			i--;
		}
		buffer[i] = '\0';

		m = split_string(buffer, field);
		if (m == 6)
		{
			imu2cam_cx = atof(field[0]);
			imu2cam_cy = atof(field[1]);
			imu2cam_cz = atof(field[2]);
			imu2cam_roll = atof(field[3]);
			imu2cam_pitch = atof(field[4]);
			imu2cam_heading = atof(field[5]);
			break;
		}
		if (m == 3)
		{
			if (first_line)
			{
				imu2cam_roll = atof(field[0]);
				imu2cam_pitch = atof(field[1]);
				imu2cam_heading = atof(field[2]);
				first_line = false;
			}
			else
			{
				imu2cam_cx = atof(field[0]);
				imu2cam_cy = atof(field[1]);
				imu2cam_cz = atof(field[2]);
			}
		}

	}
	fclose(f);
	return true;
}

int get_image_info(struct image_info *image, int size, double start_time, double end_time, const char *imgpostfile)
{
	FILE *f;
	char buffer[1024];
	char *field[64];
	double t, prev_t;
	int i, j = 0;
	bool first;
	struct image_info prev_image;
	fopen_s(&f, imgpostfile, "r");

	if (f == NULL)
	{
		printf("Failed to open file %s\n", imgpostfile);
		return 0;
	}

	first = true;
	while (fgets(buffer, 1024, f))
	{
		i = strlen(buffer);
		while (buffer[i - 1] == ' ' || buffer[i - 1] == '\t' || buffer[i - 1] == '\n' || buffer[i - 1] == '\r')
		{
			i--;
        }
		buffer[i] = '\0';

		if (split_string(buffer, field) != 13)
			continue;

		t = atof(field[3]);
		if (t < start_time)
		{
			strcpy_s(prev_image.name, sizeof(prev_image.name), field[1]);
			prev_image.center_x = atof(field[5]);
			prev_image.center_y = atof(field[6]);
			prev_image.center_z = atof(field[7]);
			prev_image.roll = atof(field[10]);
			prev_image.pitch = atof(field[11]);
			prev_image.heading = atof(field[12]);
		}
		else if (t >= start_time && t <= end_time)
		{
			if (!first && prev_t < start_time)
			{
				memcpy(&image[j], &prev_image, sizeof(prev_image));
				j++;
			}
			strcpy_s(image[j].name, sizeof(image[j].name), field[1]);
			image[j].center_x = atof(field[5]);
			image[j].center_y = atof(field[6]);
			image[j].center_z = atof(field[7]);
			image[j].roll = atof(field[10]);
			image[j].pitch = atof(field[11]);
			image[j].heading = atof(field[12]);
			j++;
		}
		else
		{
			break;
		}
		prev_t = t;
		first = false;
	}
	fclose(f);

	return j;
}

void Test1()
{
	CxImage img[5],img_o;
	img[0].Load(_T("Z:\\TechnicalResearch\\PointCloud\\1001-1-004-180129\\data\\img\\39-10011004180129113106397-0.jpg"));
	img[1].Load(_T("Z:\\TechnicalResearch\\PointCloud\\1001-1-004-180129\\data\\img\\39-10011004180129113106397-1.jpg"));
	img[2].Load(_T("Z:\\TechnicalResearch\\PointCloud\\1001-1-004-180129\\data\\img\\39-10011004180129113106397-2.jpg"));
	img[3].Load(_T("Z:\\TechnicalResearch\\PointCloud\\1001-1-004-180129\\data\\img\\39-10011004180129113106397-3.jpg"));
	img[4].Load(_T("Z:\\TechnicalResearch\\PointCloud\\1001-1-004-180129\\data\\img\\39-10011004180129113106397-4.jpg"));

	img_o.Create(8000, 4000, 24);

	double  x , y, z = -1.580;
	double  angle;
	int     ncamera;
	int i, j;
	double dr, dc;
	int row, col;
	for (i = 0; i < 8000; i++)
	{
		x = -40 + i / 100.0;

		for (j = 0; j < 4000; j++)
		{
			y = -20 + j / 100.0;

			angle = atan2(y, x) * 180 / PI;
			if (angle < 0)
				angle += 360;
			if (angle >= 360)
				angle -= 360;

			if (angle >= 0 && angle <= 36 ||
				angle > 324 && angle <= 360)
				ncamera = 0;
			if (angle > 252 && angle <= 324)
				ncamera = 1;
			if (angle > 180 && angle <= 252)
				ncamera = 2;
			if (angle > 108 && angle <= 180)
				ncamera = 3;
			if (angle > 36 && angle <= 108)
				ncamera = 4;
			ladybugXYZtoRC(context, x, y, z, ncamera, &dr, &dc, NULL);
			row = (int)(dr + 0.5);
			col = (int)(dc + 0.5);

			RGBQUAD cr = img[ncamera].GetPixelColor(col, 2048 - row - 1);
			int g = (int)(cr.rgbBlue * 0.114 + cr.rgbGreen * 0.587 + cr.rgbRed * 0.299);
			//g = g / 16 * 16;
			cr.rgbBlue = g;
			cr.rgbGreen = g;
			cr.rgbRed = g;
			img_o.SetPixelColor(i, j, cr);
		}
	}
	img_o.Save(_T("Z:\\test.jpg"), CXIMAGE_FORMAT_JPG);
}

void Test()
{
	CxImage imag0,img1,img2, img3, img4, img;
	img2.Load(_T("Z:\\TechnicalResearch\\PointCloud\\1001-1-004-180129\\data\\img\\39-10011004180129112958078-2.jpg"));
	img3.Load(_T("Z:\\TechnicalResearch\\PointCloud\\1001-1-004-180129\\data\\img\\39-10011004180129112958078-3.jpg"));
	//img2.Load(_T("Z:\\TechnicalResearch\\PointCloud\\1001-1-001-180515\\1001-1-001-180515\\data\\img\\39-10011001180515145752604-2.jpg"));
	//img3.Load(_T("Z:\\TechnicalResearch\\PointCloud\\1001-1-001-180515\\1001-1-001-180515\\data\\img\\39-10011001180515145752604-3.jpg"));

	img.Create(8000, 8000, 24);
	double  x = -1.5, y, z;
	int i, j;
	double dr, dc;
	int row, col;
	for (i = 0; i < 8000; i++)
	{
		y = -0.5 + i / 8000.0;

		for (j = 0; j < 8000; j++)
		{
			z = -0.5 + j / 8000.0;
			if (y < 0)
			{
				ladybugXYZtoRC(context, x, y, z, 2, &dr, &dc, NULL);
				row = (int)(dr + 0.5);
				col = (int)(dc + 0.5);
				img.SetPixelColor(i, j, img2.GetPixelColor(col, 2048 - row - 1));
			}
			else
			{
				ladybugXYZtoRC(context, x, y, z, 3, &dr, &dc, NULL);
				row = (int)(dr + 0.5);
				col = (int)(dc + 0.5);
				img.SetPixelColor(i, j, img3.GetPixelColor(col, 2048 - row - 1));
			}
		}
	}
	img.Save(_T("Z:\\test.jpg"), CXIMAGE_FORMAT_JPG);
}

int main(int argc, char *argv[])
{
	char       imageFile[512];
	double     start_time, end_time;
	struct     image_info image[50];
	int        image_count;

	if (argc != 8)
	{
		printf("Usage: %s [InputLasFile] [ImgPath] [ImgPostFile] [PPFile] [PGRFile] [MASKFile] [OutputDir]\n", argv[0]);
		return 0;
	}


	CLASFile lf;
	if (!lf.Load(argv[1]))
	{
		printf("Failed to load las File %s\n", argv[1]);
		return 0;
	}
	printf("Las file loaded. \n");
	if (!load_pp(argv[4]))
	{
		printf("Failed to load PP File %s\n", argv[4]);
		return 0;
	}
	printf("PP file loaded. \n");
	printf("imu2cam_cx = %f, imu2cam_cy = %f, imu2cam_cz = %f, imu2cam_roll = %f, imu2cam_pitch = %f, imu2cam_heading = %f\n", 
		imu2cam_cx, imu2cam_cy, imu2cam_cz, imu2cam_roll, imu2cam_pitch, imu2cam_heading);

	get_las_start_end_time(&lf, &start_time, &end_time);
	image_count = get_image_info(image, 50, start_time, end_time, argv[3]);
	printf("start_time = %f, end_time = %f, image_count = %d\n", start_time, end_time, image_count);

	LadybugError error;
	error = ladybugCreateContext(&context);
	_CHECK_ERROR;
	error = ladybugInitializeStreamForReading(context, argv[5], false);
	_CHECK_ERROR;
	error = ladybugGetStreamConfigFile(context, "temp.cal");
	_CHECK_ERROR;
	error = ladybugLoadConfig(context, "temp.cal");
	_CHECK_ERROR;
	double dr, dc;
	error = ladybugConfigureOutputImages(context, LADYBUG_PANORAMIC);
	_CHECK_ERROR;
	error = ::ladybugXYZtoRC(context, -10, 10, -1.5, 3, &dr, &dc, NULL);
	_CHECK_ERROR;
	double fl, cx, cy;
	// Offscreen image size needs to be set before using ladybugGetCameraUnitFocalLength and ladybugGetCameraUnitImageCenter
	const int kRectImageWidth = 400; // this doesn't affect the result
	const int kRectImageHeight = 300; // this doesn't affect the result
	error = ladybugSetOffScreenImageSize(context, LADYBUG_PANORAMIC, kRectImageWidth, kRectImageHeight);
	::ladybugGetCameraUnitFocalLength(context, 0, &fl);
	::ladybugGetCameraUnitImageCenter(context, 0, &cx, &cy);
	double x, y, z;
	x = -1.5798;
	y = 6.5432;
	z = -1.7731;
	x = -1.573;
	y = 15.334;
	z = -1.7825;
	x = 6.6924;
	y = 48.3193;
	z = -0.9704;
	build_imu_cam_maxtrix();
	transfer_coords_imu2cam(&x, &y, &z);
	error = ::ladybugXYZtoRC(context, y, -x, z, 0, &dr, &dc, NULL);

	//Test();
	//Test1();
	//return 1;

	FILE *bf;
	size_t s;

	unsigned  char *mask_buffer[5];
	int i;
	char fileName[512];
	double angle;
	for (i = 0; i < 5; i++)
	{
		sprintf_s(fileName, sizeof(fileName), "%sMASK-%d.BMP", argv[6], i);
		fopen_s(&bf, fileName, "rb");
		mask_buffer[i] = new unsigned char[308 * 2048];
		fseek(bf, 62, SEEK_SET);
		s = fread(mask_buffer[i], 1, 308 * 2048, bf);
		fclose(bf);
	}

	printf("Loading mask file finished.\n");

	CLASFile::las_header lh;

	lf.GetHeader(&lh);
	float *dist2 = new float[lh.number_point_records];
	for (i = 0; i < lh.number_point_records; i++)
	{
		dist2[i] = 1000000000;
	}

	wchar_t name[512];
	
	printf("start coloring point cloud ....\n");
	build_imu_cam_maxtrix();
	build_imu_len0_matrix();
	for (i = 0; i < image_count; i++)
	{	
		//if (i == 0)
		{
			printf("%s\n", image[i].name);
			sprintf_s(imageFile, sizeof(imageFile), "%s%s", argv[2], image[i].name);
			build_pc_imu_maxtrix(image[i].roll, image[i].pitch, image[i].heading);
			color_point_cloud(&lf, dist2, imageFile, mask_buffer, image[i].center_x, image[i].center_y, image[i].center_z);
		}
		printf("%d/%d image processed\n", i + 1, image_count);
	}

	char stroutput[512], *p, *lasname = argv[1];
	p = lasname + strlen(lasname);
	while (*p != '\\' && p > lasname)
		p--;
	sprintf_s(stroutput, sizeof(stroutput), "%s%s", argv[7], p + 1);
	lf.SaveTo(stroutput);
	delete []dist2;
	for (i = 0; i < 5; i++)
	{
		delete[]mask_buffer[i];
	}
}