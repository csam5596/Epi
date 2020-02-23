#include "vcf/filters/PassThroughLineFilter.h"

PassThroughLineFilter::PassThroughLineFilter()
	: VcfLineFilter(0.0f, 0.0f) {}

bool PassThroughLineFilter::filterVariant(const std::string &) const {
	return false;
}

std::pair<float, float> PassThroughLineFilter::getVaf(const std::string&) const {
	return std::pair<float, float>();
}

std::pair<float, float> PassThroughLineFilter::getReadCounts(const std::string&) const {
	return std::pair<float, float>();
}
