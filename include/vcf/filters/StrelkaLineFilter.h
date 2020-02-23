#ifndef EPITOPE_STRELKALINEFILTER_H
#define EPITOPE_STRELKALINEFILTER_H


#include "vcf/VcfLineFilter.h"

class StrelkaLineFilter : public VcfLineFilter {

public:
	StrelkaLineFilter(float maxNormalVaf, float minTumorVaf);

	[[nodiscard]] bool filterVariant(const std::string &line) const override;
	[[nodiscard]] std::pair<float, float> getVaf(const std::string &line) const override;
	[[nodiscard]] std::pair<float, float> getReadCounts(const std::string &line) const override;
};


#endif //EPITOPE_STRELKALINEFILTER_H
