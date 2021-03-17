#include "classes/command_line_tools.hh"

#include <cstring>

using namespace std;

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: program <option> <option>\n");
	printf("OPTIONS:\n");
	printf("    -seed-min <integer>        first seed to process\n");
	printf("    -seed-max <integer>        last seed to process\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	unsigned int seed_min = 0;
	unsigned int seed_max = 0;

	// parse command line
	for (int argi = 1; (argi < argc) && (cl_error == 0); ++argi)
	{
		if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0)
		{
			cl_error = 1;
			continue;
		}

		if (TestUIntParameter(argc, argv, argi, "-seed-min", seed_min)) continue;
		if (TestUIntParameter(argc, argv, argi, "-seed-max", seed_max)) continue;

		printf("ERROR: unknown option '%s'.\n", argv[argi]);
		cl_error = 1;
	}

	if (cl_error)
	{
		PrintUsage();
		return 1;
	}

	// TODO

	return 0;
}
