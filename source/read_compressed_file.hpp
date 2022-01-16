#ifndef H_READ_COMPRESSED_FILE_H
#define H_READ_COMPRESSED_FILE_H 1

#include <fstream>
#include <string>
#include <sstream>

using namespace std;

class autodecompress_file_t {
	public:
		autodecompress_file_t(const string& file_path);
		bool getline(string& line);
	private:
		bool compressed;
		stringstream decompressed_file_content; // buffer containing file content if file is compressed
		ifstream uncompressed_file; // file handle if file is not compressed
		string file_path;
};

class tsv_stream_t {
	public:
		tsv_stream_t(const string& s, const char d='\t'): data(&s), delimiter(d), position(0), failbit(false) {};
		tsv_stream_t& operator>>(string& out);
		tsv_stream_t& operator>>(int& out);
		bool fail() const { return failbit; };
	private:
		const string* data;
		char delimiter;
		size_t position;
		bool failbit;
};

#endif /* H_READ_COMPRESSED_FILE_H */
