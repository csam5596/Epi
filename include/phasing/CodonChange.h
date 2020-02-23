#ifndef EPITOPE_CODONCHANGE_H
#define EPITOPE_CODONCHANGE_H

#include <string>
#include <vcf/Variant.h>
#include <vep/VariantEffect.h>

struct CodonChange {

public:
	CodonChange(const Variant& variant, const VariantEffect& ve);

	[[nodiscard]] int getShift() const;

	bool operator<(const CodonChange &rhs) const;
	bool operator()(std::string& codingSequence) const;

public:
	unsigned variantId;
	unsigned variantEffectId;
	unsigned proteinPosition;
	std::string codonRef;
	std::string codonAlt;
	std::string codonReplace;
	Variant::VariantType vType;
	VariantEffect::VariantEffectType veType;

};


#endif //EPITOPE_CODONCHANGE_H
