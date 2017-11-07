#ifndef _OPTIONS_H
#define _OPTIONS_H 1

#include <cfloat>
#include <climits>
#include <string>

using namespace std;

const string HELP_CONTACT = "s.uhrig@dkfz.de";
const string ARRIBA_VERSION = "0.11.0";

string wrap_help(const string& option, const string& text, const unsigned short int max_line_width = 80);

bool output_directory_exists(const string& output_file);

bool validate_int(const char* optarg, int& value, const int min_value = LONG_MIN, const int max_value = LONG_MAX);
bool validate_int(const char* optarg, unsigned int& value, const unsigned int min_value = 0, const unsigned int max_value = LONG_MAX);

bool validate_float(const char* optarg, float& value, const float min_value = FLT_MIN, const float max_value = FLT_MAX);

#endif /* _OPTIONS_H */
