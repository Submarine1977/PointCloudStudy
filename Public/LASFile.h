#pragma once

#include <vector>
using namespace std;
class CLASFile
{
public:
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
public:
	CLASFile();
	~CLASFile();
private:
	char              m_strFileName[512];
	struct las_header m_lhHeader;
	vector<unsigned char*> m_pPointBuffer;
	const int              BUFFER_SIZE = 65536;
private:
	void Clear();
public:
	bool LoadHeader(const char *strFile);
	bool Load(const char *strFile);
	bool GetPointInfo(unsigned int i, unsigned char* pBuffer);
	bool SetPointInfo(unsigned int i, unsigned char* pBuffer);
	bool GetHeader(struct las_header *pHeader);
	bool SetHeader(struct las_header lasHeader);
	bool Save();
	bool SaveTo(const char *strFileName);
};

