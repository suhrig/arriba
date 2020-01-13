#ifndef _H_READ_COMPRESSED_FILE_H
#define _H_READ_COMPRESSED_FILE_H 1

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

#endif /* _H_READ_COMPRESSED_FILE_H */
