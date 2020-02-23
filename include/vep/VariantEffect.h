#ifndef EPITOPE_VARIANTEFFECT_H
#define EPITOPE_VARIANTEFFECT_H

#include <string>
#include <optional>

#include "vcf/Contig.h"

class VariantEffect {

public:
	enum VariantEffectType {
		SNV, INFRAME_INSERTION, INFRAME_DELETION, FRAMESHIFT_VARIANT
	};

	static std::optional<VariantEffect> parse(const std::string& line);
	static VariantEffectType getType(const std::string& refCodon, const std::string& altCodon);

	[[nodiscard]] const std::string& getTypeAsString() const;

private:

	VariantEffect(unsigned variantId, Contig &chr, unsigned pos, std::string &gene, std::string &transcript, int strand,
			      std::string &ensp, std::string &codonRef, std::string &codonAlt, unsigned proteinPos);

public:
	unsigned variantId;
	unsigned variantEffectId;
	Contig chromosome;
	unsigned position;
	std::string gene;
	std::string transcript; // Feature
	int strand;

	std::string proteinId;
	std::string codonRef;
	std::string codonAlt;
	unsigned proteinPosition;


	VariantEffectType type;
};


#endif //EPITOPE_VARIANTEFFECT_H
