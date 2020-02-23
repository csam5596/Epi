#ifndef EPITOPE_PASSTHROUGHLINEFILTER_H
#define EPITOPE_PASSTHROUGHLINEFILTER_H


#include <vcf/VcfLineFilter.h>

class PassThroughLineFilter : public VcfLineFilter {

public:
	PassThroughLineFilter();

	[[nodiscard]] bool filterVariant(const std::string &line) const override;
	[[nodiscard]] std::pair<float, float> getVaf(const std::string &line) const override;
	[[nodiscard]] std::pair<float, float> getReadCounts(const std::string &line) const override;
};


#endif //EPITOPE_PASSTHROUGHLINEFILTER_H
