/**
 * @file
 *
 * This is the main file. It contains functions to parse the commandline
 * arguments, read files etc.
 *
 * @brief The main file
 * @author Fabian Klötzl
 *
 * @section License
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <err.h>
#include <errno.h>
#include <getopt.h>
#include <string.h>
#include <unistd.h>

#include "config.h"
#include "global.h"
#include "io.h"
#include "process.h"
#include "sequence.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int FLAGS = flags::none;
int THREADS = 1;
double RANDOM_ANCHOR_PROP = 0.05;

void usage(int);
void version(void);

int main(int argc, char *argv[])
{
	int version_flag = 0;
	auto ref_name = std::string("longest");

	static struct option long_options[] = {
		{"version", no_argument, &version_flag, 1},
		{"help", no_argument, NULL, 'h'},
		{"verbose", no_argument, NULL, 'v'},
		{"threads", required_argument, NULL, 't'},
		{0, 0, 0, 0}};

#ifdef _OPENMP
	// Use all available processors by default.
	THREADS = omp_get_num_procs();
#endif

	// parse arguments
	while (1) {
		int c = getopt_long(argc, argv, "hr:t:v", long_options, NULL);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: break;
			case 'h': usage(EXIT_SUCCESS); break;
			case 'r': {
				ref_name = std::string(optarg);
				break;
			}
			case 't': {
#ifdef _OPENMP
				errno = 0;
				char *end;
				long unsigned int threads = strtoul(optarg, &end, 10);

				if (errno || end == optarg || *end != '\0') {
					warnx("Expected a number for -t argument, but '%s' was "
						  "given. Ignoring -t argument.",
						  optarg);
					break;
				}

				if (threads > (long unsigned int)omp_get_num_procs()) {
					warnx(
						"The number of threads to be used, is greater then the "
						"number of available processors; Ignoring -t %lu "
						"argument.",
						threads);
					break;
				}

				THREADS = threads;
#else
				warnx(
					"This version of phylonium was built without OpenMP "
					"and thus "
					"does not support multi threading. Ignoring -t argument.");
#endif
				break;
			}
			case 'v':
				FLAGS |= FLAGS & flags::verbose ? flags::extra_verbose
												: flags::verbose;
				break;
			case '?': /* intentional fall-through */
			default: usage(EXIT_FAILURE); break;
		}
	}

	if (version_flag) {
		version();
	}

	argc -= optind;
	argv += optind;

	// at least one file name must be given
	if (argc == 0) {
		errx(1, "At least one filename needs to be supplied.");
	}

	auto file_names = std::vector<std::string>(argv, argv + argc);

	if (file_names.size() < 2) {
		file_names.push_back("-"); // if no files are supplied, read from stdin
	}

	if (ref_name != "longest" &&
		std::find(file_names.begin(), file_names.end(), ref_name) ==
			file_names.end()) {
		file_names.push_back(ref_name);
	}

	// at max `argc` many files have to be read.
	auto genomes = std::vector<genome>();
	genomes.reserve(argc);

	// read all genomes
	std::transform(file_names.begin(), file_names.end(),
				   std::back_inserter(genomes), read_genome);

	auto it = max_element(begin(genomes), end(genomes),
						  [](const genome &a, const genome &b) {
							  return a.get_length() < b.get_length();
						  });
	auto derp = std::find(file_names.begin(), file_names.end(), ref_name);
	const auto &ref = (ref_name == "longest")
						  ? *it
						  : *(derp - file_names.begin() + genomes.begin());

	process(ref, genomes);

	return 0;
}

/**@brief
 * Prints the usage to stdout and then exits successfully.
 */
void usage(int status)
{
	const char str[] = {
		"Usage: phylonium [-lv] [-t INT] FILES...\n"
		"\tFILES... can be any sequence of FASTA files. If no files are "
		"supplied, stdin is used instead.\n"
		"Options:\n"
		"  -r longest|FILENAME   Use the sequence from FILENAME as reference;"
		"default: longest\n"
		"  -v, --verbose     Prints additional information\n"
#ifdef _OPENMP
		"  -t, --threads <INT> \n"
		"                    The number of threads to be used; by default, all "
		"available processors are used\n"
#endif
		"  -h, --help        Display this help and exit\n"
		"      --version     Output version information and acknowledgments\n"};

	fprintf(status == EXIT_SUCCESS ? stdout : stderr, "%s", str);
	exit(status);
}

/**
 * This function just prints the version string and then aborts
 * the program.
 */
void version(void)
{
	const char str[] = {
		"phylonium " VERSION "\n"
		"Copyright (C) 2017 Fabian Klötzl\n"
		"License GPLv3+: GNU GPL version 3 or later "
		"<http://gnu.org/licenses/gpl.html>\n"
		"This is free software: you are free to change and redistribute it.\n"
		"There is NO WARRANTY, to the extent permitted by law.\n\n"};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}
