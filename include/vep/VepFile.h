#ifndef EPITOPE_VEPFILE_H
#define EPITOPE_VEPFILE_H

#include <string>
#include <vector>
#include <set>
#include "vcf/VcfFile.h"
#include "VariantEffect.h"

class VepFile {

public:

	enum Column {
		UPLOADED_VARIATION, LOCATION, ALLELE, GENE, FEATURE, FEATURE_TYPE, CONSEQUENCE, CDNA_POSITION,
		CDS_POSITION, PROTEIN_POSITION, AMINO_ACIDS, CODONS, EXISTING_VARIATION, EXTRA
	};

	static VepFile parse(const std::string& vepFile);

	const std::set<std::string>& getTranscriptIds() const;
	const std::optional<std::vector<VariantEffect> const*> getVariantEffects(unsigned variantId) const;
	const VariantEffect& getVariantEffectById(unsigned id) const;

private:
	VepFile(std::unordered_map<unsigned, std::vector<VariantEffect>> &effects, std::set<std::string>& ids,
			std::unordered_map<unsigned, VariantEffect>& effectsById);

private:
	std::unordered_map<unsigned, std::vector<VariantEffect>> variantEffects;
	std::unordered_map<unsigned, VariantEffect> effectsById;
	std::set<std::string> transcriptIds;
};


#endif //EPITOPE_VEPFILE_H
