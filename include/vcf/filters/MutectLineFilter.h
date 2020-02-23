#ifndef EPITOPE_MUTECTLINEFILTER_H
#define EPITOPE_MUTECTLINEFILTER_H

#include "vcf/VcfLineFilter.h"

class MutectLineFilter : public VcfLineFilter {

public:
	MutectLineFilter(float maxNormalVaf, float minTumorVaf);

	[[nodiscard]] bool filterVariant(const std::string &line) const override;
	[[nodiscard]] std::pair<float, float> getVaf(const std::string &line) const override;
	[[nodiscard]] std::pair<float, float> getReadCounts(const std::string &line) const override;

};


#endif //EPITOPE_MUTECTLINEFILTER_H
