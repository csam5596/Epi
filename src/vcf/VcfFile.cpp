#include <Logging.h>
#include "vcf/VcfFile.h"

VcfFile VcfFile::parse(const std::string &path) {

	LOG_DEBUG("Loading variants from " + path + ".");

	std::unordered_map<unsigned long, Variant> variants;

	VcfLineIterator it(path, "None", 0.0f, 0.0f);
	while(it.isGood()){
		Variant v = Variant::parse(it.lineBuffer);
		variants.emplace(v.id, std::move(v));
		it.advance();
	}

	LOG_DEBUG(std::to_string(variants.size()) + " variants loaded.");

	return VcfFile(it.header, variants);
}

const Variant& VcfFile::getVariantById(unsigned long id) const {
	return variants.at(id);
}

const std::unordered_map<unsigned long, Variant>& VcfFile::getVariants() {
	return variants;
}

VcfFile::VcfFile(VcfHeader &header, std::unordered_map<unsigned long, Variant> &variants)
	: vcfHeader(std::move(header)), variants(std::move(variants)){}