#include <cinttypes>

class evo_model
{
protected:
	enum {
		AtoA,
		AtoC,
		AtoG,
		AtoT,
		CtoC,
		CtoG,
		CtoT,
		GtoG,
		GtoT,
		TtoT,
		MUTCOUNTS,
		CtoA = AtoC,
		GtoA = AtoG,
		GtoC = CtoG,
		TtoA = AtoT,
		TtoC = CtoT,
		TtoG = GtoT
	};


	int counts[16] = {0};

public:
	static int hash(char c) noexcept;

	evo_model() = default;

	void account(char a, char b) noexcept;
	void account(const char *stra, const char *strb, size_t length) noexcept;
	void account_rev(const char *stra, const char *strb, size_t b_offset,  size_t length) noexcept;
	evo_model &operator+=(const evo_model &other) noexcept;
	size_t total() const noexcept;
	double estimate_raw() const noexcept;
	double estimate_JC() const noexcept;
};
