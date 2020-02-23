#ifndef EPITOPE_VCFLINEFILTER_H
#define EPITOPE_VCFLINEFILTER_H

#include <string>
#include <memory>
#include <map>
#include "Contig.h"


class VcfLineFilter {

public:
	static std::unique_ptr<VcfLineFilter> createLineFilter(const std::string& variantCaller, float normalVaf, float tumorVaf);

	VcfLineFilter(float normalVaf, float tumorVaf);

	[[nodiscard]] virtual bool filterVariant(const std::string& line) const = 0;
	[[nodiscard]] virtual std::pair<float,float> getVaf(const std::string& line) const = 0;
	[[nodiscard]] virtual std::pair<float,float> getReadCounts(const std::string& line) const = 0;
	[[nodiscard]] bool chromosomeValid(const std::string& chromosome) const;

protected:
	[[nodiscard]] bool vafValid(float normalVaf, float tumorVaf) const;

protected:
	float maxNormalVaf;
	float minTumorVaf;

};


#endif //EPITOPE_VCFLINEFILTER_H
