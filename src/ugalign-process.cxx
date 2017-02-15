/**
 * @file
 * @brief This file contains various distance methods.
 */

#include <algorithm>
#include <iostream>

#include <math.h>
#include <stdio.h>
#include "esa.h"
#include "global.h"
#include "homology.h"
#include "io.h"
#include "process.h"
#include "sequence.h"

#include <time.h>

double shuprop(size_t x, double g, size_t l);
void print_pile(const sequence &, size_t, std::vector<homology> &, size_t,
				size_t);

/**
 * @brief Calculates the minimum anchor length.
 *
 * Given some parameters calculate the minimum length for anchors according
 * to the distribution from Haubold et al. (2009).
 *
 * @param p - The probability with which an anchor is allowed to be random.
 * @param g - The the relative amount of GC in the subject.
 * @param l - The length of the subject.
 * @returns The minimum length of an anchor.
 */
size_t minAnchorLength(double p, double g, size_t l)
{
	size_t x = 1;

	double prop = 0.0;
	while (prop < 1 - p) {
		prop = shuprop(x, g / 2, l);
		x++;
	}

	return x;
}

/**
 * @brief Calculates the binomial coefficient of n and k.
 *
 * We used to use gsl_sf_lnchoose(xx,kk) for this functionality.
 * After all, why implement something that has already been done?
 * Well, the reason is simplicity: GSL is used for only this one
 * function and the input (n<=20) is not even considered big.
 * Hence its much easier to have our own implementation and ditch
 * the GSL dependency even if that means our code is a tiny bit
 * less optimized and slower.
 *
 * @param n - The n part of the binomial coefficient.
 * @param k - analog.
 * @returns (n choose k)
 */
size_t binomial_coefficient(size_t n, size_t k)
{
	if (n <= 0 || k > n) {
		return 0;
	}

	if (k == 0 || k == n) {
		return 1;
	}

	if (k > n - k) {
		k = n - k;
	}

	size_t res = 1;

	for (size_t i = 1; i <= k; i++) {
		res *= n - k + i;
		res /= i;
	}

	return res;
}

/**
 * @brief Given `x` this function calculates the probability of a shustring
 * with a length less than `x`.
 *
 * Let X be the longest shortest unique substring (shustring) at any position.
 * Then this function computes P{X <= x} with respect to the given parameter
 * set. See Haubold et al. (2009).
 *
 * @param x - The maximum length of a shustring.
 * @param g - The the half of the relative amount of GC in the DNA.
 * @param l - The length of the subject.
 * @returns The probability of a certain shustring length.
 */
double shuprop(size_t x, double p, size_t l)
{
	double xx = (double)x;
	double ll = (double)l;
	size_t k;

	double s = 0.0;

	for (k = 0; k <= x; k++) {
		double kk = (double)k;
		double t = pow(p, kk) * pow(0.5 - p, xx - kk);

		s += pow(2, xx) * (t * pow(1 - t, ll)) *
			 (double)binomial_coefficient(x, k);
		if (s >= 1.0) {
			s = 1.0;
			break;
		}
	}

	return s;
}

/**
 * @brief Compute homologies between two sequences.
 */
std::vector<homology> anchor_homologies(const esa_s *C, double gc,
										const sequence *query_ptr)
{
	std::vector<homology> hv;

	size_t border = C->len / 2;

	const char *query = query_ptr->S.c_str();
	size_t query_length = query_ptr->len;
	lcp_inter_t inter;

	size_t last_pos_Q = 0;
	size_t last_pos_S = 0;
	size_t last_length = 0;
	// This variable indicates that the last anchor was the right anchor of a
	// pair.
	size_t last_was_right_anchor = 0;

	size_t this_pos_Q = 0;
	size_t this_pos_S;
	size_t this_length;

	size_t threshold =
		minAnchorLength(1 - sqrt(1 - RANDOM_ANCHOR_PROP), gc, C->len);

	homology current = {query_ptr};
	current.str = "";

	// Iterate over the complete query.
	while (this_pos_Q < query_length) {
		inter =
			get_match_cached(C, query + this_pos_Q, query_length - this_pos_Q);

		this_length = inter.l <= 0 ? 0 : inter.l;
		if (this_pos_Q + this_length > query_length) {
			this_length = query_length - this_pos_Q;
		}

		if (inter.i == inter.j && this_length >= threshold) {
			// anchor
			this_pos_S = C->SA[inter.i];

			if (this_pos_S > last_pos_S &&
				this_pos_Q - last_pos_Q == this_pos_S - last_pos_S &&
				(this_pos_S < border) == (last_pos_S < border)) {
				// right anchor

				current.length +=
					this_pos_Q - last_pos_Q - last_length + this_length;

				current.str += query_ptr->S.substr(
					last_pos_Q + last_length,
					this_pos_Q - last_pos_Q - last_length + this_length);

				last_was_right_anchor = 1;
			} else if (this_pos_S > last_pos_S &&
					   this_pos_Q - last_pos_Q < this_pos_S - last_pos_S &&
					   (this_pos_Q - last_pos_Q) + 200 >
						   this_pos_S - last_pos_S &&
					   (this_pos_S < border) == (last_pos_S < border)) {

				//
				auto gap_S = this_pos_S - last_pos_S - last_length;
				auto gap_Q = this_pos_Q - last_pos_Q - last_length;
				auto aln = asymalign(
					std::string{C->S + last_pos_S + last_length, gap_S},
					query_ptr->S.substr(last_pos_Q + last_length, gap_Q));

				// current.length +=
				// 	this_pos_Q - last_pos_Q - last_length + this_length;

				current.str += aln;
				current.str += query_ptr->S.substr(this_pos_Q, this_length);

				last_was_right_anchor = 1;
			} else {
				// no anchor pair, left anchor!
				if (last_was_right_anchor || last_length / 2 >= threshold) {
					// push last
					current.finalise(border);
					hv.push_back(std::move(current));
				}

				current = {query_ptr, this_pos_S, this_pos_Q, this_length,
						   query_ptr->S.substr(this_pos_Q, this_length)};

				// this is not right anchor
				last_was_right_anchor = 0;
			}

			// Cache values for later
			last_pos_Q = this_pos_Q;
			last_pos_S = this_pos_S;
			last_length = this_length;
		}

		// Advance
		this_pos_Q += this_length + 1;
	}

	// Very special case: The sequences are identical
	if (last_length >= query_length) {
		current = {query_ptr, last_pos_S, (size_t)0, query_length,
				   query_ptr->S};
	}

	if (last_was_right_anchor || last_length / 2 >= threshold) {
		current.finalise(border);
		hv.push_back(std::move(current));
	}

	return hv;
}

void reference_align(genome &reference_genome, std::vector<sequence> &sequences)
{
	printf("##maf version=1\n");

	esa_s E;

	clock_t t1 = clock();

	sequence subject = join(reference_genome.contigs);
	subject.name = reference_genome.name; // FIXME

	subject.prepare_subject();
	if (esa_init(&E, subject.RS.c_str(), subject.RSlen)) {
		errx(1, "Failed to create index for %s.", subject.name.c_str());
	}

	clock_t t2 = clock();

	std::vector<homology> all;
// now compare every other sequence to the subject
#pragma omp parallel for num_threads(THREADS)
	for (size_t j = 0; j < sequences.size(); j++) {
		sequences[j].prepare_query();

		// TODO: Provide a nicer progress indicator.
		if (FLAGS & F_EXTRA_VERBOSE) {
#pragma omp critical
			{
				fprintf(stderr, "comparing %zu and %zu\n", 0, j);
			}
		}

		std::vector<homology> hvlocal =
			anchor_homologies(&E, subject.gc, &sequences[j]);

#pragma omp critical
		{
			auto inserter = std::back_inserter(all);
			std::move(hvlocal.begin(), hvlocal.end(), inserter);
		}
	}

	esa_free(&E);

	clock_t t3 = clock();

	std::sort(begin(all), end(all), [](const homology &a, const homology &b) {
		return a.index_subject_relative < b.index_subject_relative;
	});

	clock_t t4 = clock();

	// walk the line
	std::vector<homology> pile;

	auto next_homo = begin(all);

	size_t curpos = 0;			 // w.r.t the joined genome
	size_t curpos_in_contig = 0; // w.r.t the current contig
	auto current_contig = reference_genome.contigs.begin();

	// loop while there are homologies waiting or on the pile
	while (next_homo < all.end() || pile.size()) {
		// the current alignment block can at most span the whole contig
		size_t length = current_contig->len - curpos_in_contig;

		// if there is a next homology, check its beginning
		if (next_homo < all.end()) {
			length =
				std::min(length, next_homo->index_subject_relative - curpos);
		}

		// the following loop could be optimized away if the pile were sorted by
		// end positions.
		for (const auto &entry : pile) {
			length = std::min(length, entry.index_subject_relative +
										  entry.length - curpos);
		}

		// print the whole pile from curpos to nextpos
		print_pile(*current_contig, curpos_in_contig, pile, curpos, length);

		size_t nextpos = curpos + length;

		// delete completely printed homologies from the pile
		auto is_done = [=](const homology &entry) {
			return nextpos >= entry.index_subject_relative + entry.length;
		};
		pile.erase(std::remove_if(pile.begin(), pile.end(), is_done),
				   pile.end());

		// add to the pile homologies which start at the current point
		auto split = next_homo;
		while (split < all.end() && nextpos == split->index_subject_relative) {
			split++;
		}

		// move from the `all` container into pile.
		auto inserter = std::back_inserter(pile);
		std::move(next_homo, split, inserter);
		next_homo = split;

		// jump to the next block
		curpos_in_contig += length;
		curpos += length;
		if (curpos_in_contig >= current_contig->len) {
			curpos++; // skip '!'
			current_contig++;
			curpos_in_contig = 0;
		}
	}

	clock_t t5 = clock();
	if (FLAGS & F_VERBOSE) {
		static const char str[] = {"index building: %lf\n"
								   "matching: %lf\n"
								   "sorting: %lf\n"
								   "output: %lf\n"};
		fprintf(stderr, str, ((double)(t2 - t1)) / CLOCKS_PER_SEC,
				((double)(t3 - t2)) / CLOCKS_PER_SEC,
				((double)(t4 - t3)) / CLOCKS_PER_SEC,
				((double)(t5 - t4)) / CLOCKS_PER_SEC);
	}
}

void print_pile(const sequence &subject, size_t curpos_in_contig,
				std::vector<homology> &pile, size_t curpos, size_t length)
{
	if (pile.empty() || !length) return;

	printf("a\ns %s %5zu %5zu + %5zu %.*s\n", subject.name.c_str(),
		   curpos_in_contig, length, subject.len, (int)length,
		   subject.S.c_str() + curpos_in_contig);

	for (auto &entry : pile) {
		const auto &query = *entry.query;

		size_t curpos_in_homo = curpos - entry.index_subject_relative;
		// TODO: don't copy the substring
		auto stripe = entry.str.substr(curpos_in_homo, length);
		auto gaps = count_gaps(stripe);

		// recount non_gaps in current shit
		if (entry.orientation == homology::FORWARD) {
			size_t offset =
				curpos_in_homo + entry.index_query - entry.gaps_so_far;
			printf("s %s %5zu %5zu + %5zu %s\n", query.name.c_str(), offset,
				   length - gaps, query.len, stripe.c_str());
		} else {
			size_t from_end_homo_end =
				query.len - (entry.index_query + entry.length_in_query);

			size_t revoffset =
				from_end_homo_end + curpos_in_homo - entry.gaps_so_far;

			printf("s %s %5zu %5zu - %5zu %s\n", query.name.c_str(), revoffset,
				   length - gaps, query.len, stripe.c_str());
		}

		entry.gaps_so_far += gaps;
	}

	// each multiple alignment ends with a blank line
	printf("\n");
}