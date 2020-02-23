#ifndef EPITOPE_VCFFILE_H
#define EPITOPE_VCFFILE_H

#include <vector>
#include <string>
#include <vcf/VcfHeader.h>
#include <vcf/Variant.h>

#include "VcfLineIterator.h"
#include "vcf/filters/MutectLineFilter.h"


class VcfFile {

public:

	enum Column {
		CHROM = 0, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, NORMAL, TUMOR
	};

	static VcfFile parse(const std::string& path);
	const Variant& getVariantById(unsigned long id) const;
	const std::unordered_map<unsigned long, Variant>& getVariants();

private:
	VcfFile(VcfHeader& header, std::unordered_map<unsigned long, Variant>& variants);

private:
	VcfHeader vcfHeader;
	std::unordered_map<unsigned long, Variant> variants;

};
#endif //EPITOPE_VCFFILE_H
