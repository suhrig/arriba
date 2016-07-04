#ifndef _H_READ_COMPRESSED_FILE_H
#define _H_READ_COMPRESSED_FILE_H 1

#include <string>
#include <sstream>

using namespace std;

void autodecompress_file(const string& file_path, stringstream& file_content);

#endif /* _H_READ_COMPRESSED_FILE_H */
