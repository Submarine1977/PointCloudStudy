// DumpLASFile.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include <string.h>

struct las_header
{
	char signature[4];
	unsigned short source_id, global_encoding;
	unsigned char  guid[16];
	unsigned char  major_version;
	unsigned char  minor_version;
	char system_identifies[32];
	char generating_software[32];
	unsigned short creation_day;
	unsigned short creation_year;
	unsigned short header_size;
	unsigned long offset_to_point_data;
	unsigned long number_variable_length_records;
	unsigned char point_data_format_id;
	unsigned short point_data_record_length;
	unsigned long number_point_records;
	unsigned long number_point_records_return[5];
	double x_scale, y_scale, z_scale;
	double x_offset, y_offset, z_offset;
	double maxx, minx, maxy, miny, maxz, minz;

};

int read_las_header(FILE *f, struct las_header *lh)
{
	if (fread(&lh->signature, sizeof(lh->signature), 1, f) < 1 ||
		fread(&lh->source_id, sizeof(lh->source_id), 1, f) < 1 ||
		fread(&lh->global_encoding, sizeof(lh->global_encoding), 1, f) < 1 ||
		fread(&lh->guid, sizeof(lh->guid), 1, f) < 1 ||
		fread(&lh->major_version, sizeof(lh->major_version), 1, f) < 1 ||
		fread(&lh->minor_version, sizeof(lh->minor_version), 1, f) < 1 ||
		fread(&lh->system_identifies, sizeof(lh->system_identifies), 1, f) < 1 ||
		fread(&lh->generating_software, sizeof(lh->generating_software), 1, f) < 1 ||
		fread(&lh->creation_day, sizeof(lh->creation_day), 1, f) < 1 ||
		fread(&lh->creation_year, sizeof(lh->creation_year), 1, f) < 1 ||
		fread(&lh->header_size, sizeof(lh->header_size), 1, f) < 1 ||
		fread(&lh->offset_to_point_data, sizeof(lh->offset_to_point_data), 1, f) < 1 ||
		fread(&lh->number_variable_length_records, sizeof(lh->number_variable_length_records), 1, f) < 1 ||
		fread(&lh->point_data_format_id, sizeof(lh->point_data_format_id), 1, f) < 1 ||
		fread(&lh->point_data_record_length, sizeof(lh->point_data_record_length), 1, f) < 1 ||
		fread(&lh->number_point_records, sizeof(lh->number_point_records), 1, f) < 1 ||
		fread(&lh->number_point_records_return, sizeof(lh->number_point_records_return), 1, f) < 1 ||
		fread(&lh->x_scale, sizeof(lh->x_scale), 1, f) < 1 ||
		fread(&lh->y_scale, sizeof(lh->y_scale), 1, f) < 1 ||
		fread(&lh->z_scale, sizeof(lh->z_scale), 1, f) < 1 ||
		fread(&lh->x_offset, sizeof(lh->x_offset), 1, f) < 1 ||
		fread(&lh->y_offset, sizeof(lh->y_offset), 1, f) < 1 ||
		fread(&lh->z_offset, sizeof(lh->z_offset), 1, f) < 1 ||
		fread(&lh->maxx, sizeof(lh->maxx), 1, f) < 1 ||
		fread(&lh->minx, sizeof(lh->minx), 1, f) < 1 ||
		fread(&lh->maxy, sizeof(lh->maxy), 1, f) < 1 ||
		fread(&lh->miny, sizeof(lh->miny), 1, f) < 1 ||
		fread(&lh->maxz, sizeof(lh->maxz), 1, f) < 1 ||
		fread(&lh->minz, sizeof(lh->minz), 1, f) < 1)
	{
		printf("Failed to las header, file size is less than las hearder length\n");
		return -1;
	}
	return 0;
}

void dump_las_header(char * las_file)
{
	struct las_header lh;
	FILE *f;
	fopen_s(&f, las_file, "rb");
	if (f == NULL)
	{
		printf("Failed to open file %s\n", las_file);
		return;
	}
	//cannot read the whole structure together because bit alignment issue
	if(read_las_header(f, &lh) < 0)
	{
		printf("Failed to las header %s\n", las_file);
		return;
	}
	printf("\
signature         = %s \n\
source_id         = %d \n\
global_encoding   = %d \n\
guid              = (%x%x%x%x%x%x%x%x%x%x%x%x%x%x%x%x) \n\
version           = %d.%d \n\
system_identifies = %s \n\
char generating_software = %s \n\
creation_day = %d \n\
creation_year = %d \n\
header_size = %d \n\
offset_to_point_data = %d \n\
number_variable_length_records = %d \n\
point_data_format_id = %d \n\
point_data_record_length = %d \n\
number_point_records = %d \n\
number_point_records_return = (%d,%d,%d,%d,%d) \n\
x_scale =  %f\n\
y_scale =  %f\n\
z_scale =  %f\n\
x_offset =  %f\n\
y_offset =  %f\n\
z_offset =  %f\n\
maxx =  %f\n\
minx =  %f\n\
maxy =  %f\n\
miny =  %f\n\
maxz =  %f\n\
minz =  %f\n\
",
lh.signature,
lh.source_id,
lh.global_encoding,
lh.guid[0], lh.guid[1], lh.guid[2], lh.guid[3], lh.guid[4], lh.guid[5], lh.guid[6], lh.guid[7], lh.guid[8], lh.guid[9], lh.guid[10], lh.guid[11], lh.guid[12], lh.guid[13], lh.guid[14], lh.guid[15],
lh.major_version, lh.minor_version,
lh.system_identifies,
lh.generating_software,
lh.creation_day,
lh.creation_year,
lh.header_size,
lh.offset_to_point_data,
lh.number_variable_length_records,
lh.point_data_format_id,
lh.point_data_record_length,
lh.number_point_records,
lh.number_point_records_return[0], lh.number_point_records_return[1], lh.number_point_records_return[2], lh.number_point_records_return[3], lh.number_point_records_return[4],
lh.x_scale, lh.y_scale, lh.z_scale,
lh.x_offset, lh.y_offset, lh.z_offset,
lh.maxx, lh.minx, lh.maxy, lh.miny, lh.maxz, lh.minz);
}

void dump_las_point(char *las_file, int start_pos, int count)
{
	struct las_header lh;
	FILE *f;
	int i, t;
	unsigned char buffer[64];
	fopen_s(&f, las_file, "rb");
	if (f == NULL)
	{
		printf("Failed to open file %s\n", las_file);
		return;
	}
	//cannot read the whole structure together because bit alignment issue
	if (read_las_header(f, &lh) < 0)
	{
		printf("Failed to read las header %s\n", las_file);
		return;
	}
	
	fseek(f, lh.offset_to_point_data + lh.point_data_record_length * start_pos, SEEK_SET);
	for (i = 0; i < count; i++)
	{
		if (fread(buffer, 1, lh.point_data_record_length, f) < lh.point_data_record_length)
		{
			printf("Failed to read las point %s, %d\n", las_file, start_pos + i);
			return;
		}
		long x, y, z;
		
		short intensity;
		unsigned char rnse; //return number, numbr of returns, scan direction flag, edge of flight line;
		unsigned char classification, scan_angle_rank, user_data;
		unsigned short point_source_id;
		double gps_time;
		unsigned short r, g, b;
		
		t = 0;
		memcpy(&x, buffer + t, sizeof(x));
		t += sizeof(x);
		memcpy(&y, buffer + t, sizeof(y));
		t += sizeof(y);
		memcpy(&z, buffer + t, sizeof(z));
		t += sizeof(z);
		memcpy(&intensity, buffer + t, sizeof(intensity));
		t += sizeof(intensity);
		memcpy(&rnse, buffer + t, sizeof(rnse));
		t += sizeof(rnse);
		memcpy(&classification, buffer + t, sizeof(classification));
		t += sizeof(classification);
		memcpy(&scan_angle_rank, buffer + t, sizeof(scan_angle_rank));
		t += sizeof(scan_angle_rank);
		memcpy(&user_data, buffer + t, sizeof(user_data));
		t += sizeof(user_data);
		memcpy(&point_source_id, buffer + t, sizeof(point_source_id));
		t += sizeof(point_source_id);
		memcpy(&gps_time, buffer + t, sizeof(gps_time));
		t += sizeof(gps_time);
		memcpy(&r, buffer + t, sizeof(r));
		t += sizeof(r);
		memcpy(&g, buffer + t, sizeof(g));
		t += sizeof(g);
		memcpy(&b, buffer + t, sizeof(b));
		t += sizeof(b);

		printf("%f,%f,%f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%d,%d,%d\n", 
			x * lh.x_scale + lh.x_offset, y * lh.y_scale + lh.y_offset, z * lh.z_scale + lh.z_offset, intensity,
			(rnse & 0x7), (rnse >> 3) & 0x7, (rnse & 64)?1:0, (rnse & 128)?1:0, 
			classification, scan_angle_rank, user_data, point_source_id, gps_time, r, g, b);

	}

}

int main(int argc, char * argv[])
{
	if (argc == 3 && strcmp(argv[2], "-h") == 0)
	{
		dump_las_header(argv[1]);
	}
	else if (argc == 5 && strcmp(argv[2], "-p") == 0)
	{
		dump_las_point(argv[1], atoi(argv[3]), atoi(argv[4]));
	}
    return 0;
}

