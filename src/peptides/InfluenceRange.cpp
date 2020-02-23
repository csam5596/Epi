#include "peptides/InfluenceRange.h"

#include <iostream>

InfluenceRange::InfluenceRange(const CodonChange &codonChange)
	: variantEffectId(codonChange.variantEffectId), veType(codonChange.veType){

	unsigned offset = 0;
	while(codonChange.codonRef.size() > offset &&
		  codonChange.codonAlt.size() > offset &&
		  codonChange.codonRef.at(offset) == codonChange.codonAlt.at(offset)){
		offset++;
	}

	// get index of first base that is different (0-based)
	firstStart = (codonChange.proteinPosition - 1) * 3 + offset;

	switch(veType){

		case VariantEffect::SNV:
		case VariantEffect::INFRAME_DELETION:
			lastStart = firstStart;
			break;
		case VariantEffect::INFRAME_INSERTION:
			lastStart = firstStart + (codonChange.getShift() - 1);
			break;
		case VariantEffect::FRAMESHIFT_VARIANT:
			lastStart = std::numeric_limits<long>::max()/2;
			break;

	}
}

void InfluenceRange::shift(int amount) {
	firstStart += amount;
	if(veType != VariantEffect::FRAMESHIFT_VARIANT) lastStart += amount;
}

void InfluenceRange::convertPositions() {
	firstStart /= 3;
	if(veType != VariantEffect::FRAMESHIFT_VARIANT){
		lastStart /= 3;
	}
}