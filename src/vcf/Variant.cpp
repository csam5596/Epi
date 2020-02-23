#include <sstream>
#include <util/VariantUtil.h>
#include "vcf/Variant.h"

Variant Variant::parse(const std::string &line) {

	int vcCount = VariantUtil::getVariantCallerCount(line);
	float vafNormal = VariantUtil::getNormalVaf(line);
	float vafTumor = VariantUtil::getTumorVaf(line);

	std::istringstream ss(line);

	std::string chr;
	unsigned position;
	unsigned id;
	std::string ref;
	std::string alt;
	std::string info;

	ss >> chr >> position >> id >> ref >> alt;

	// ignore quality and filter
	ss.ignore(9999, '\t');
	ss.ignore(9999, '\t');
	ss.ignore(9999,'\t');

	ss >> info;

	VariantType t = VariantUtil::getVariantType(line);
	float rcNormal = VariantUtil::getReadCountNormal(line);
	float rcTumor = VariantUtil::getReadCountTumor(line);;

	return Variant(Contig(chr), position, id, ref, alt, t, vcCount, vafNormal, vafTumor, rcNormal, rcTumor);
}

Variant::Variant(Contig &&chr, unsigned position, unsigned id, std::string &ref, std::string &alt, VariantType type,
				 int count, float vafNormal, float vafTumor, float rcNormal, float rcTumor)
	: contig(std::move(chr)), position(position), id(id), ref(std::move(ref)), alt(std::move(alt)), type(type),
	  vcCount(count), vafNormal(vafNormal), vafTumor(vafTumor), rcNormal(rcNormal), rcTumor(rcTumor) {

	name = chr.getName() + ":" + std::to_string(position) + ":" + this->ref + ":" + this->alt;
}

const std::string& Variant::getName() const {
	return name;
}