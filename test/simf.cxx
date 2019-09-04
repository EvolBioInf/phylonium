/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2019 © Fabian Klötzl
 */
/**
 * This program can create genome sequences with a specific distance.
 */

#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <random>
#include <string>

using namespace std;

void usage(int);

template <typename OutputIt>
void print_seq(OutputIt &, unsigned, unsigned, size_t, size_t, double);

int main(int argc, char *argv[])
{
	random_device rd{};
	auto seed = rd();
	size_t length = 1000;
	size_t line_length = 70;
	bool raw = 0;
	auto prefix = std::string{};

	auto seqs = vector<double>{0};

	int check;
	while ((check = getopt(argc, argv, "d:hl:L:p:rs:")) != -1) {
		if (check == 'd') {
			seqs.push_back(stod(optarg));
		} else if (check == 'h') {
			usage(0);
		} else if (check == 'l') {
			length = stoi(optarg);
		} else if (check == 'L') {
			line_length = stoi(optarg);
		} else if (check == 'p') {
			prefix = optarg;
		} else if (check == 'r') {
			raw = true;
		} else if (check == 's') {
			seed = static_cast<unsigned int>(stol(optarg));
			if (seed == 0) {
				seed = rd();
			}
		} else {
			usage(1);
		}
	}

	if (seqs.size() < 2) {
		seqs.push_back(0.1);
	}

	if (!raw) {
		for (auto &dist : seqs) {
			auto d = dist;
			auto p = 0.75 - 0.75 * exp(-(4.0 / 3.0) * d);
			dist = p;
		}
	}

	auto base_seed = seed;

	for (size_t i = 0; i < seqs.size(); i++) {
		if (prefix != "") {
			auto file_name = prefix + to_string(i) + ".fasta";
			auto out = std::ofstream(file_name);
			out << ">S" << i << " (base_seed: " << base_seed << ")" << endl;
			print_seq(out, base_seed, seed++, length, line_length, seqs[i]);
		} else {
			cout << ">S" << i << " (base_seed: " << base_seed << ")" << endl;
			print_seq(cout, base_seed, seed++, length, line_length, seqs[i]);
		}
	}

	return 0;
}

static auto ACGT = "ACGT";
static auto NO_A = "CGT";
static auto NO_C = "AGT";
static auto NO_G = "ACT";
static auto NO_T = "ACG";

template <typename OutputIt>
void print_seq(OutputIt &out, unsigned base_seed, unsigned mut_seed,
			   size_t length, size_t line_length, double divergence)
{
	char line[line_length + 1];
	line[line_length] = '\0';

	auto base_rand = default_random_engine{base_seed};
	auto base_dist = uniform_int_distribution<int>{0, 3};
	auto base_acgt = [&] { return ACGT[base_dist(base_rand)]; };

	auto mut_rand = default_random_engine{mut_seed};
	auto mut_dist = uniform_real_distribution<double>{0, 1};
	auto mut = bind(mut_dist, mut_rand);
	auto mut_acgt = uniform_int_distribution<int>{0, 2};
	auto mutate = [&](char c) {
		int idx = mut_acgt(mut_rand);
		switch (c) {
			case 'A': return NO_A[idx];
			case 'C': return NO_C[idx];
			case 'G': return NO_G[idx];
			case 'T': return NO_T[idx];
			default: return 'X';
		}
	};

	double nucleotides = (double)length;
	double mutations = nucleotides * divergence;

	for (size_t i = length, j; i > 0; i -= j) {
		j = min(line_length, i);

		for (size_t k = 0; k < j; k++) {
			char c = base_acgt();

			if (mut() < mutations / nucleotides) {
				c = mutate(c);
				mutations--;
			}

			line[k] = c;
			nucleotides--;
		}

		line[j] = '\0';
		out << line << endl;
	}
}

void usage(int exit_code)
{
	const static char *str = {"usage: test_fasta [-d dist...] [-l length] [-L "
							  "line length] [-p prefix] [-r raw] [-s seed]\n"};

	if (exit_code == 0) {
		cout << str;
	} else {
		cerr << str;
	}
	exit(exit_code);
}
