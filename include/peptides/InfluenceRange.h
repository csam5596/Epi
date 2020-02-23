#ifndef EPITOPE_INFLUENCERANGE_H
#define EPITOPE_INFLUENCERANGE_H


#include <phasing/CodonChange.h>

struct InfluenceRange {

public:
	explicit InfluenceRange(const CodonChange& codonChange);

	void shift(int amount);
	void convertPositions(); // change coding sequence positions to protein positions

public:
	unsigned variantEffectId;
	VariantEffect::VariantEffectType veType;

	// in terms of base positions
	long firstStart;
	long lastStart;
};


#endif //EPITOPE_INFLUENCERANGE_H
