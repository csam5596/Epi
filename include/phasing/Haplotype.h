#ifndef EPITOPE_HAPLOTYPE_H
#define EPITOPE_HAPLOTYPE_H

#include <string>
#include <vector>
#include <vep/VepFile.h>
#include "Haploblock.h"

#include <Types.h>
#include <phasing/CodonChange.h>

class Haplotype {

	
public:
	Haplotype(const std::string& filePrefix, const std::string& rna);

	void computeProteinChanges(const VcfFile& vcf, const VepFile& vep, unsigned minRcAbs, float minRcRel);
	void addHomozygousVariants(const std::string& vcf, const VepFile& vep);
	[[nodiscard]] const std::unordered_map<ENSP, std::vector<std::set<CodonChange>>>& getProteinChanges() const;
	[[nodiscard]] std::optional<float> getRnaVaf(unsigned long variantId) const;
	[[nodiscard]] std::optional<unsigned> getRcRNA(unsigned long variantId) const;

private:
	void parseRnaCounts(const std::string& file, const std::string& rna);

private:
	std::vector<Haploblock> haploblocks;
	std::unordered_map<std::string, std::pair<unsigned, unsigned>> rnaPhaseCounts;
	std::unordered_map<ENSP, std::vector<std::set<CodonChange>>> changesByProtein;
	std::unordered_map<unsigned, float> vafRNAByVariant;
	std::unordered_map<unsigned, unsigned> rcRNAByVariant;
};


#endif //EPITOPE_HAPLOTYPE_H
