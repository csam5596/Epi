#include "vcf/VcfLineFilter.h"

#include <regex>
#include <vcf/filters/MutectLineFilter.h>
#include <vcf/filters/HaplotypeCallerFilter.h>
#include <vcf/filters/VarscanLineFilter.h>
#include <vcf/filters/StrelkaLineFilter.h>

#include <iostream>
#include <vcf/filters/PassThroughLineFilter.h>

std::unique_ptr<VcfLineFilter> VcfLineFilter::createLineFilter(const std::string &variantCaller, float normalVaf, float tumorVaf) {
	if(variantCaller == "Mutect2"){
		return std::make_unique<MutectLineFilter>(normalVaf, tumorVaf);
	}else if(variantCaller == "Varscan"){
		return std::make_unique<VarscanLineFilter>(normalVaf, tumorVaf);
	}else if(variantCaller == "Strelka") {
		return std::make_unique<StrelkaLineFilter>(normalVaf, tumorVaf);
	}else if(variantCaller == "HaplotypeCaller") {
		return std::make_unique<HaplotypeCallerFilter>();
	}else if(variantCaller == "None"){
		return std::make_unique<PassThroughLineFilter>();
	}else{
		throw std::runtime_error("[VcfLineFilter::createLineFilter] Unknown variant caller \"" + variantCaller + "\"");
	}
}

VcfLineFilter::VcfLineFilter(float normalVaf, float tumorVaf)
	: maxNormalVaf(normalVaf), minTumorVaf(tumorVaf) {}

bool VcfLineFilter::chromosomeValid(const std::string &chromosome) const {
	const static std::regex validChromosome("^chr([1-9]|1[0-9]|2[0-3]|X|Y)$");

	std::smatch match;
	return std::regex_search(chromosome, match, validChromosome);
}

bool VcfLineFilter::vafValid(float normalVaf, float tumorVaf) const {
	return tumorVaf >= minTumorVaf && normalVaf <= maxNormalVaf;
}