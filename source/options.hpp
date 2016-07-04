#ifndef _OPTIONS_H
#define _OPTIONS_H 1

#include <string>

using namespace std;

const string HELP_CONTACT = "s.uhrig@dkfz.de";
const string ARIBA_VERSION = "0.4";

string wrap_help(const string option, const string text, const unsigned short int max_line_width = 80);

#endif /* _OPTIONS_H */
