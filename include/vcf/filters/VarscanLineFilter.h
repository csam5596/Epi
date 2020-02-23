#ifndef EPITOPE_VARSCANLINEFILTER_H
#define EPITOPE_VARSCANLINEFILTER_H


#include "vcf/VcfLineFilter.h"

class VarscanLineFilter : public VcfLineFilter {

public:
	VarscanLineFilter(float maxNormalVaf, float minTumorVaf);

	[[nodiscard]] bool filterVariant(const std::string &line) const override;
	[[nodiscard]] std::pair<float, float> getVaf(const std::string &line) const override;
	[[nodiscard]] std::pair<float, float> getReadCounts(const std::string &line) const override;
};


#endif //EPITOPE_VARSCANLINEFILTER_H
