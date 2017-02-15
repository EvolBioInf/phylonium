/**
 * @file
 * @brief This file contains the definitions for various io methods.
 */

#include <err.h>
#include <fcntl.h>
#include <unistd.h>

#include <pfasta.h>

#include <string>
#include <vector>

#include "io.h"
#include "sequence.h"

/** @brief extracts the genome name from a file path
 *
 * We try to be clever about the genome name. Given the file
 * path we extract just the file name. ie. path/file.ext -> file
 * This obviously fails on Windows.
 *
 * @param s_file_name - The file path
 * @returns just the genome name
 */
std::string extract_genome(const std::string &s_file_name)
{
	// find the last path separator
	auto left = s_file_name.rfind('/');
	left = (left == std::string::npos) ? 0 : left + 1;
	// left is the position one of to the right of the path separator

	// find the extension
	auto right = s_file_name.find('.', left);
	right = (right == std::string::npos) ? s_file_name.size() : right;

	// copy only the file name, not its path or extension
	return s_file_name.substr(left, right - left);
}

/**
 * @brief This function reads sequences from a file.
 * @param file_name - The file to read.
 * @param dsa - (output parameter) An array that holds found sequences.
 */
std::vector<sequence> read_fasta(std::string s_file_name, std::string prefix)
{
	std::vector<sequence> sequences{};
	const char *file_name = s_file_name.c_str();

	int file_descriptor =
		s_file_name != "-" ? open(file_name, O_RDONLY) : STDIN_FILENO;

	if (file_descriptor < 0) {
		warn("%s", file_name);
		throw "derp.";
	}

	int l;
	pfasta_file pf;

	if ((l = pfasta_parse(&pf, file_descriptor)) != 0) {
		warnx("%s: %s", file_name, pfasta_strerror(&pf));
		goto fail;
	}

	pfasta_seq ps;
	while ((l = pfasta_read(&pf, &ps)) == 0) {
		sequences.emplace_back(prefix + ps.name, filter_nucl(ps.seq));
		pfasta_seq_free(&ps);
	}

	if (l < 0) {
		warnx("%s: %s", file_name, pfasta_strerror(&pf));
		pfasta_seq_free(&ps);
	}

fail:
	pfasta_free(&pf);
	close(file_descriptor);

	return sequences;
}

genome read_genome(std::string file_name)
{
	std::string species{extract_genome(file_name.c_str())};

	return genome{species, read_fasta(file_name, species + ".")};
}
