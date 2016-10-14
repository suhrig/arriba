#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include "sam.h"
#include "read_compressed_file.hpp"

using namespace std;

void autodecompress_file(const string& file_path, stringstream& file_content) {

	if (file_path.length() >= 3 && file_path.substr(file_path.length() - 3) == ".gz") {

		// allocate some memory for buffered reading
		const unsigned int buffer_size = 64*1024;
		char* buffer = (char*) malloc(buffer_size+1);
		if (buffer == NULL) {
			cerr << "ERROR: failed to allocate memory for decompression." << endl;
			exit(1);
		}

		// open compressed file
		BGZF* compressed_file;
		compressed_file = bgzf_open(file_path.c_str(), "rb");
		if(compressed_file == NULL) {
			cerr << "ERROR: failed to open/decompress file '" << file_path << "'." << endl;
			exit(1);
		}

		// read data from file and put it into stringstream object
		int bytes_read;
		do {
			bytes_read = bgzf_read(compressed_file, buffer, buffer_size);
			if (bytes_read < 0) {
				cerr << "ERROR: failed to decompress file '" << file_path << "'." << endl;
				exit(1);
			}
			buffer[bytes_read] = '\0';
			file_content << buffer;
			if (!file_content.good()) {
				cerr << "ERROR: failed to load file '" << file_path << "' into memory." << endl;
				exit(1);
			}
		} while (bytes_read == buffer_size);

		// free resources
		bgzf_close(compressed_file);
		free(buffer);

	} else { // file is not compressed
		
		// copy file content to stringstream object
		ifstream uncompressed_file(file_path);
		if (!uncompressed_file.good()) {
			cerr << "ERROR: failed to open file '" << file_path << "'." << endl;
			exit(1);
		}
		file_content << uncompressed_file.rdbuf();
		if (!file_content.good()) {
			cerr << "ERROR: failed to load file '" << file_path << "' into memory." << endl;
			exit(1);
		}
		uncompressed_file.close();
	}
}

