#include <string>
#include "catch.hpp"
#include "sequence.h"

TEST_CASE("Sequence basics")
{
	auto s = sequence("Name", "ACGTACGT");

	REQUIRE(s.get_name() == "Name");
	REQUIRE(s.get_nucl() == "ACGTACGT");
	REQUIRE(s.size() == 8);
}

TEST_CASE("Revcomp")
{
	REQUIRE(reverse("") == "");
	REQUIRE(reverse("A") == "T");
	REQUIRE(reverse("C") == "G");
	REQUIRE(reverse("G") == "C");
	REQUIRE(reverse("T") == "A");
	REQUIRE(reverse("ACGTACGT") == "ACGTACGT");

	auto str = std::string("TACGATCGATCGAAAGCTAGTTCGCCCCGAGATA");
	auto rc = "TATCTCGGGGCGAACTAGCTTTCGATCGATCGTA";
	REQUIRE(reverse(str) == rc);
	REQUIRE(reverse(reverse(str)) == str);
}

TEST_CASE("Filtering nucleotides")
{
	REQUIRE(filter_nucl("") == "");
	REQUIRE(filter_nucl("A") == "A");
	REQUIRE(filter_nucl("C") == "C");
	REQUIRE(filter_nucl("G") == "G");
	REQUIRE(filter_nucl("T") == "T");
	REQUIRE(filter_nucl("!") == "");

	auto str = std::string("TACGATCGATCGAAAGCTAGTTCGCCCCGAGATA");
	REQUIRE(filter_nucl(str) == str);

	REQUIRE(filter_nucl("tacgatc!gatc!gaa__agctagttcgcc#ccgagata") == str);
}
