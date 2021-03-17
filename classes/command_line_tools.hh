#ifndef _command_line_tools_hh_
#define _command_line_tools_hh_

#include <string>

extern int cl_error;

extern bool TestBoolParameter(int argc, const char **argv, int &argi, const char *tag, bool &param1, bool &param2);

extern bool TestBoolParameter(int argc, const char **argv, int &argi, const char *tag, bool &param);

extern bool TestUIntParameter(int argc, const char **argv, int &argi, const char *tag, unsigned int &param);

extern bool TestDoubleParameter(int argc, const char **argv, int &argi, const char *tag, double &param1, double &param2);

extern bool TestDoubleParameter(int argc, const char **argv, int &argi, const char *tag, double &param);

extern bool TestStringParameter(int argc, const char **argv, int &argi, const char *tag, std::string &param1, std::string &param2);

extern bool TestStringParameter(int argc, const char **argv, int &argi, const char *tag, std::string &param);

#endif
