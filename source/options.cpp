#include <cstring>
#include <libgen.h>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include "options.hpp"

using namespace std;

string wrap_help(const string& option, const string& text, const unsigned short int max_line_width) {
	istringstream iss(text);
	string result = " " + option + "  ";
	string indent(result.length(), ' ');
	unsigned short int line_width = result.length();
	while (iss && iss.tellg() != -1) {
		if (iss.str().substr(iss.tellg(), 1) == "\n") {
			result += "\n" + indent;
			line_width = indent.length();
		}
		string word;
		iss >> word;
		if (word.length() + line_width > max_line_width) {
			result += "\n" + indent;
			line_width = indent.length();
		}
		result += word + " ";
		line_width += word.length();
	}
	return result + "\n\n";
}

bool output_directory_exists(const string& output_file) {
	char* output_file_c_str = strdup(output_file.c_str()); // we need to make a copy, because dirname does not take const argument
	char* output_directory = dirname(output_file_c_str);
	struct stat file_info;
	return (stat(output_directory, &file_info) == 0 && (file_info.st_mode & S_IFDIR));
}

bool validate_int(const char* optarg, int& value, const int min_value, const int max_value) {
	value = atoi(optarg);
	if (optarg != "0" && value == 0)
		return false;
	return value >= min_value && value <= max_value;
}

bool validate_int(const char* optarg, unsigned int& value, const unsigned int min_value, const unsigned int max_value) {
	int signed_int;
	return validate_int(optarg, signed_int, min_value, max_value);
}

bool validate_float(const char* optarg, float& value, const float min_value, const float max_value) {
	value = atof(optarg);
	if (optarg != "0" && value == 0)
		return false;
	return value >= min_value && value <= max_value;
}

