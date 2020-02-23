#ifndef EPITOPE_HAPLOTYPECALLERFILTER_H
#define EPITOPE_HAPLOTYPECALLERFILTER_H


#include <vcf/VcfLineFilter.h>

class HaplotypeCallerFilter : public VcfLineFilter {

public:
	HaplotypeCallerFilter();

	[[nodiscard]] bool filterVariant(const std::string &line) const override;
	[[nodiscard]] std::pair<float, float> getVaf(const std::string &line) const override;
	[[nodiscard]] std::pair<float, float> getReadCounts(const std::string &line) const override;
};


#endif //EPITOPE_HAPLOTYPECALLERFILTER_H
