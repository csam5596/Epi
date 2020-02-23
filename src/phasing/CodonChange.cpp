#include "phasing/CodonChange.h"
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <Logging.h>

CodonChange::CodonChange(const Variant &variant, const VariantEffect &ve)
	: variantId(variant.id), variantEffectId(ve.variantEffectId), proteinPosition(ve.proteinPosition),
	  codonRef(ve.codonRef), codonAlt(ve.codonAlt), vType(variant.type), veType(ve.type) {

	// if indel remove equal and unaffected bases at the end from both codons
	if(codonRef.length() != codonAlt.length()){
		while(*codonRef.rbegin() == *codonAlt.rbegin()){
			if(std::isupper(*codonRef.rbegin()) || std::isupper(*codonAlt.rbegin())) break;
			codonRef.pop_back();
			codonAlt.pop_back();
		}
	}
	codonReplace = codonAlt;
	boost::to_upper(codonReplace);
}

bool CodonChange::operator<(const CodonChange &rhs) const {
	return proteinPosition < rhs.proteinPosition;
}

bool CodonChange::operator()(std::string &codingSequence) const {

	unsigned long cdsPosition = (proteinPosition - 1) * 3;
	if(codingSequence.length() < cdsPosition){
		LOG_WARN("CDS length " + std::to_string(codingSequence.length()) + " < " + std::to_string(cdsPosition));
		return false;
	}
	std::string_view foundCodon = std::string_view(codingSequence).substr(cdsPosition, codonRef.length());
	if(!boost::iequals(foundCodon, codonRef)){
		LOG_WARN("Expected: [" + codonRef + "] \tBut found: [" + std::string(foundCodon) + "]");
		return false;
	}
	codingSequence.replace(cdsPosition, codonRef.length(), codonReplace);
	return true;
}

int CodonChange::getShift() const {
	return static_cast<int>(codonRef.length()) - codonAlt.length();
}