#include "LASFile.h"
#include "stdio.h"
#include "stdlib.h"

CLASFile::CLASFile()
{
}

CLASFile::~CLASFile()
{
	Clear();
}
void CLASFile::Clear()
{
	unsigned int i;
	for (i = 0; i < m_pPointBuffer.size(); i++)
	{
		delete[]m_pPointBuffer[i];
	}
	m_pPointBuffer.clear();
}

bool CLASFile::LoadHeader(const char *strFile)
{
	FILE *f = NULL;
	fopen_s(&f, strFile, "rb");
	if (f == NULL)
	{
		printf("Failed to open file %s\n", strFile);
		return false;
	}
	if (fread(&m_lhHeader.signature, sizeof(m_lhHeader.signature), 1, f) < 1 ||
		fread(&m_lhHeader.source_id, sizeof(m_lhHeader.source_id), 1, f) < 1 ||
		fread(&m_lhHeader.global_encoding, sizeof(m_lhHeader.global_encoding), 1, f) < 1 ||
		fread(&m_lhHeader.guid, sizeof(m_lhHeader.guid), 1, f) < 1 ||
		fread(&m_lhHeader.major_version, sizeof(m_lhHeader.major_version), 1, f) < 1 ||
		fread(&m_lhHeader.minor_version, sizeof(m_lhHeader.minor_version), 1, f) < 1 ||
		fread(&m_lhHeader.system_identifies, sizeof(m_lhHeader.system_identifies), 1, f) < 1 ||
		fread(&m_lhHeader.generating_software, sizeof(m_lhHeader.generating_software), 1, f) < 1 ||
		fread(&m_lhHeader.creation_day, sizeof(m_lhHeader.creation_day), 1, f) < 1 ||
		fread(&m_lhHeader.creation_year, sizeof(m_lhHeader.creation_year), 1, f) < 1 ||
		fread(&m_lhHeader.header_size, sizeof(m_lhHeader.header_size), 1, f) < 1 ||
		fread(&m_lhHeader.offset_to_point_data, sizeof(m_lhHeader.offset_to_point_data), 1, f) < 1 ||
		fread(&m_lhHeader.number_variable_length_records, sizeof(m_lhHeader.number_variable_length_records), 1, f) < 1 ||
		fread(&m_lhHeader.point_data_format_id, sizeof(m_lhHeader.point_data_format_id), 1, f) < 1 ||
		fread(&m_lhHeader.point_data_record_length, sizeof(m_lhHeader.point_data_record_length), 1, f) < 1 ||
		fread(&m_lhHeader.number_point_records, sizeof(m_lhHeader.number_point_records), 1, f) < 1 ||
		fread(&m_lhHeader.number_point_records_return, sizeof(m_lhHeader.number_point_records_return), 1, f) < 1 ||
		fread(&m_lhHeader.x_scale, sizeof(m_lhHeader.x_scale), 1, f) < 1 ||
		fread(&m_lhHeader.y_scale, sizeof(m_lhHeader.y_scale), 1, f) < 1 ||
		fread(&m_lhHeader.z_scale, sizeof(m_lhHeader.z_scale), 1, f) < 1 ||
		fread(&m_lhHeader.x_offset, sizeof(m_lhHeader.x_offset), 1, f) < 1 ||
		fread(&m_lhHeader.y_offset, sizeof(m_lhHeader.y_offset), 1, f) < 1 ||
		fread(&m_lhHeader.z_offset, sizeof(m_lhHeader.z_offset), 1, f) < 1 ||
		fread(&m_lhHeader.maxx, sizeof(m_lhHeader.maxx), 1, f) < 1 ||
		fread(&m_lhHeader.minx, sizeof(m_lhHeader.minx), 1, f) < 1 ||
		fread(&m_lhHeader.maxy, sizeof(m_lhHeader.maxy), 1, f) < 1 ||
		fread(&m_lhHeader.miny, sizeof(m_lhHeader.miny), 1, f) < 1 ||
		fread(&m_lhHeader.maxz, sizeof(m_lhHeader.maxz), 1, f) < 1 ||
		fread(&m_lhHeader.minz, sizeof(m_lhHeader.minz), 1, f) < 1)
	{
		printf("Failed to las header, file size is less than las hearder length\n");
		fclose(f);
		return false;
	}
	fclose(f);
	return true;
}
bool CLASFile::Load(const char *strFile)
{
	unsigned int i;
	Clear();
	if (!LoadHeader(strFile))
	{
		return false;
	}
	FILE *f = NULL;
	fopen_s(&f, strFile, "rb");
	if (f == NULL)
	{
		printf("Failed to open file %s\n", strFile);
		return false;
	}

	fseek(f, m_lhHeader.offset_to_point_data, SEEK_SET);
	unsigned char *p;
	size_t m;
	for (i = 0; i < m_lhHeader.number_point_records; i += BUFFER_SIZE)
	{
		p = new unsigned char[BUFFER_SIZE * m_lhHeader.point_data_record_length];
		m = fread(p, m_lhHeader.point_data_record_length, BUFFER_SIZE, f);
		m_pPointBuffer.push_back(p);
	}
	fclose(f);
	strcpy_s(m_strFileName, strFile);
	return true;
}

bool CLASFile::GetPointInfo(unsigned int i, unsigned char* pBuffer)
{
	if (i < 0 || i >= m_lhHeader.number_point_records)
	{
		printf("Index out of range\n");
		return false;
	}
	memcpy(pBuffer, m_pPointBuffer[i / BUFFER_SIZE] + (i % BUFFER_SIZE) * m_lhHeader.point_data_record_length, m_lhHeader.point_data_record_length);
	return true;
}

bool CLASFile::SetPointInfo(unsigned int i, unsigned char* pBuffer)
{
	if (i < 0 || i >= m_lhHeader.number_point_records)
	{
		printf("Index out of range\n");
		return false;
	}
	memcpy(m_pPointBuffer[i / BUFFER_SIZE] + (i % BUFFER_SIZE) * m_lhHeader.point_data_record_length, pBuffer, m_lhHeader.point_data_record_length);
	return true;
}

bool CLASFile::GetHeader(struct las_header *pHeader)
{
	memcpy(pHeader, &m_lhHeader, sizeof(struct las_header));
	return true;
}
bool CLASFile::SetHeader(struct las_header lasHeader)
{
	memcpy(&m_lhHeader, &lasHeader, sizeof(struct las_header));
	return true;
}
bool CLASFile::Save()
{
	return SaveTo(m_strFileName);
}
bool CLASFile::SaveTo(const char *strFileName)
{
	FILE *f;
	unsigned int i;
	fopen_s(&f, strFileName, "wb");
	if (f == NULL)
	{
		printf("Failed to open file %s\n", strFileName);
		return false;
	}
	if (fwrite(&m_lhHeader.signature, sizeof(m_lhHeader.signature), 1, f) < 1 ||
		fwrite(&m_lhHeader.source_id, sizeof(m_lhHeader.source_id), 1, f) < 1 ||
		fwrite(&m_lhHeader.global_encoding, sizeof(m_lhHeader.global_encoding), 1, f) < 1 ||
		fwrite(&m_lhHeader.guid, sizeof(m_lhHeader.guid), 1, f) < 1 ||
		fwrite(&m_lhHeader.major_version, sizeof(m_lhHeader.major_version), 1, f) < 1 ||
		fwrite(&m_lhHeader.minor_version, sizeof(m_lhHeader.minor_version), 1, f) < 1 ||
		fwrite(&m_lhHeader.system_identifies, sizeof(m_lhHeader.system_identifies), 1, f) < 1 ||
		fwrite(&m_lhHeader.generating_software, sizeof(m_lhHeader.generating_software), 1, f) < 1 ||
		fwrite(&m_lhHeader.creation_day, sizeof(m_lhHeader.creation_day), 1, f) < 1 ||
		fwrite(&m_lhHeader.creation_year, sizeof(m_lhHeader.creation_year), 1, f) < 1 ||
		fwrite(&m_lhHeader.header_size, sizeof(m_lhHeader.header_size), 1, f) < 1 ||
		fwrite(&m_lhHeader.offset_to_point_data, sizeof(m_lhHeader.offset_to_point_data), 1, f) < 1 ||
		fwrite(&m_lhHeader.number_variable_length_records, sizeof(m_lhHeader.number_variable_length_records), 1, f) < 1 ||
		fwrite(&m_lhHeader.point_data_format_id, sizeof(m_lhHeader.point_data_format_id), 1, f) < 1 ||
		fwrite(&m_lhHeader.point_data_record_length, sizeof(m_lhHeader.point_data_record_length), 1, f) < 1 ||
		fwrite(&m_lhHeader.number_point_records, sizeof(m_lhHeader.number_point_records), 1, f) < 1 ||
		fwrite(&m_lhHeader.number_point_records_return, sizeof(m_lhHeader.number_point_records_return), 1, f) < 1 ||
		fwrite(&m_lhHeader.x_scale, sizeof(m_lhHeader.x_scale), 1, f) < 1 ||
		fwrite(&m_lhHeader.y_scale, sizeof(m_lhHeader.y_scale), 1, f) < 1 ||
		fwrite(&m_lhHeader.z_scale, sizeof(m_lhHeader.z_scale), 1, f) < 1 ||
		fwrite(&m_lhHeader.x_offset, sizeof(m_lhHeader.x_offset), 1, f) < 1 ||
		fwrite(&m_lhHeader.y_offset, sizeof(m_lhHeader.y_offset), 1, f) < 1 ||
		fwrite(&m_lhHeader.z_offset, sizeof(m_lhHeader.z_offset), 1, f) < 1 ||
		fwrite(&m_lhHeader.maxx, sizeof(m_lhHeader.maxx), 1, f) < 1 ||
		fwrite(&m_lhHeader.minx, sizeof(m_lhHeader.minx), 1, f) < 1 ||
		fwrite(&m_lhHeader.maxy, sizeof(m_lhHeader.maxy), 1, f) < 1 ||
		fwrite(&m_lhHeader.miny, sizeof(m_lhHeader.miny), 1, f) < 1 ||
		fwrite(&m_lhHeader.maxz, sizeof(m_lhHeader.maxz), 1, f) < 1 ||
		fwrite(&m_lhHeader.minz, sizeof(m_lhHeader.minz), 1, f) < 1)
	{
		printf("Failed to las header, file size is less than las hearder length\n");
		fclose(f);
		return false;
	}
	for (i = 0; i < m_pPointBuffer.size() - 1; i++)
	{
		fwrite(m_pPointBuffer[i], m_lhHeader.point_data_record_length, BUFFER_SIZE, f);
	}
	fwrite(m_pPointBuffer[i], m_lhHeader.point_data_record_length, (m_lhHeader.number_point_records - i * BUFFER_SIZE), f);
	fclose(f);
	return true;
}