#include <string>
#include <iostream>
#include "bgzf.h"
#include "sam.h"
#include "common.hpp"
#include "read_compressed_file.hpp"

using namespace std;

autodecompress_file_t::autodecompress_file_t(const string& file_path) {

	compressed = file_path.length() >= 3 && file_path.substr(file_path.length() - 3) == ".gz";
	if (compressed) {

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
		if (compressed_file == NULL) {
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
			decompressed_file_content << buffer;
			if (!decompressed_file_content.good()) {
				cerr << "ERROR: failed to load file '" << file_path << "' into memory." << endl;
				exit(1);
			}
		} while (bytes_read == (unsigned int) buffer_size);

		// free resources
		bgzf_close(compressed_file);
		free(buffer);

	} else {
		// the file is not compressed => only open it now and read it on-demand
		this->file_path = file_path;
		uncompressed_file.open(file_path);
		if (!uncompressed_file.is_open()) {
			cerr << "ERROR: failed to open file '" << file_path << "'." << endl;
			exit(1);
		}
	}
}

bool autodecompress_file_t::getline(string& line) {
	if (compressed) {
		if (std::getline(decompressed_file_content, line))
			return true;
		else
			return false;
	} else {
		if (std::getline(uncompressed_file, line)) {
			return true;
		} else {
			if (uncompressed_file.bad()) {
				cerr << "ERROR: failed to load file '" << file_path << "' into memory." << endl;
				exit(1);
			} else if (uncompressed_file.eof()) {
				uncompressed_file.close();
			}
			return false;
		}
	}
}

tsv_stream_t& tsv_stream_t::operator>>(string& out) {
	if (position >= data->size()) {
		failbit = true;
	} else {
		size_t start_position = position;
		position = data->find(delimiter, start_position);
		out = data->substr(start_position, position - start_position);
		position++; // skip delimiter
	}
	return *this;
}

tsv_stream_t& tsv_stream_t::operator>>(int& out) {
	if (position >= data->size()) {
		failbit = true;
	} else {
		size_t start_position = position;
		position = data->find(delimiter, start_position);
		if (!str_to_int(data->substr(start_position, position - start_position).c_str(), out))
			failbit = true;
		position++; // skip delimiter
	}
	return *this;
}
