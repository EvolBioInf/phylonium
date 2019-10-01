#include <string>
#include "catch.hpp"
#include "process.h"

void filter_overlaps_max(std::vector<homology> &pile);

bool operator==(const homology &a, const homology &b)
{
	if (a.start() != b.start()) return false;
	if (a.end() != b.end()) return false;
	if (a.start_query() != b.start_query()) return false;
	if (a.end_query() != b.end_query()) return false;

	return true;
}

TEST_CASE("Homology basics")
{
	auto A = homology{0, 0, 10};
	auto B = homology{1, 1, 10};

	REQUIRE(A.starts_left_of(B));
	REQUIRE(!A.ends_left_of(B));
	REQUIRE(A.overlaps(B));

	auto C = homology{10, 10, 10};
	REQUIRE(A.starts_left_of(C));
	REQUIRE(A.ends_left_of(C));
	REQUIRE(!A.overlaps(C));

	// Query coordinate doesn't matter
	A = homology{0, 23456, 10};
	B = homology{1, 678, 10};
	C = homology{10, 987, 10};

	REQUIRE(A.starts_left_of(B));
	REQUIRE(!A.ends_left_of(B));
	REQUIRE(A.overlaps(B));
	REQUIRE(A.starts_left_of(C));
	REQUIRE(A.ends_left_of(C));
	REQUIRE(!A.overlaps(C));

	auto D = homology{0, 0, 100}.trim(0, 10);
	A = homology{0, 0, 10};

	REQUIRE(D.start() == A.start());
	REQUIRE(D.end() == A.end());
	REQUIRE(D.start_query() == A.start_query());
	REQUIRE(D.end_query() == A.end_query());
}

TEST_CASE("Homology filtering")
{
	// two possible beginnings
	auto pile = std::vector<homology>{{0, 0, 10}, {1, 1, 3}};

	filter_overlaps_max(pile);
	REQUIRE(pile.size() == 1);
	REQUIRE(pile[0] == homology{0, 0, 10});

	// overlap in the middle
	pile = std::vector<homology>{
		{0, 0, 10}, {10, 10, 10}, {10, 10, 20}, {40, 40, 5}};

	auto expected =
		std::vector<homology>{{0, 0, 10}, {10, 10, 20}, {40, 40, 5}};

	filter_overlaps_max(pile);
	REQUIRE(pile == expected);

	// two possible endings
	pile = std::vector<homology>{
		{0, 0, 10}, {10, 10, 10}, {10, 10, 20}, {40, 40, 5}, {42, 42, 2}};

	filter_overlaps_max(pile);
	REQUIRE(pile == expected);

	// two chains
	pile = std::vector<homology>{{10, 10, 10}, {0, 0, 10},   {20, 20, 10},
								 {5, 5, 10},   {15, 15, 10}, {25, 25, 10},
								 {30, 30, 10}};

	expected = std::vector<homology>{
		{0, 0, 10}, {10, 10, 10}, {20, 20, 10}, {30, 30, 10}};

	std::sort(begin(pile), end(pile),
			  [](const homology &self, const homology &other) {
				  return self.starts_left_of(other);
			  });
	filter_overlaps_max(pile);
	REQUIRE(pile == expected);
}

TEST_CASE("Complete deletion")
{
	auto homologies = std::vector<std::vector<homology>>{
		{{10, 10, 10}, {110, 110, 20}, {220, 220, 10}, {260, 260, 10}},
		{{10, 10, 10}, {120, 120, 20}, {200, 200, 100}},
		{{0, 0, 300}, {300, 300, 100}}};

	auto expected = std::vector<std::vector<homology>>{
		{{10, 10, 10}, {120, 120, 10}, {220, 220, 10}, {260, 260, 10}},
		{{10, 10, 10}, {120, 120, 10}, {220, 220, 10}, {260, 260, 10}},
		{{10, 10, 10}, {120, 120, 10}, {220, 220, 10}, {260, 260, 10}}};

	REQUIRE(complete_delete(homologies) == expected);
	REQUIRE(complete_delete(expected) == expected);

	// query coordinates vary per sequence
	homologies = std::vector<std::vector<homology>>{
		{{10, 110, 10}, {110, 210, 20}, {220, 320, 10}, {260, 460, 10}},
		{{10, 510, 10}, {120, 620, 20}, {200, 700, 100}},
		{{0, 0, 300}, {300, 300, 100}}};

	expected = std::vector<std::vector<homology>>{
		{{10, 110, 10}, {120, 220, 10}, {220, 320, 10}, {260, 460, 10}},
		{{10, 510, 10}, {120, 620, 10}, {220, 720, 10}, {260, 760, 10}},
		{{10, 10, 10}, {120, 120, 10}, {220, 220, 10}, {260, 260, 10}}};

	REQUIRE(complete_delete(homologies) == expected);
}
