#include <vector>
#include <iterator>
#include <sstream>
#include <vcf/VcfFile.h>
#include "vcf/filters/HaplotypeCallerFilter.h"

HaplotypeCallerFilter::HaplotypeCallerFilter()
	: VcfLineFilter(0.0f, 0.0f) {}

bool HaplotypeCallerFilter::filterVariant(const std::string &line) const {

	std::istringstream ss(line);
	std::istream_iterator<std::string> beg(ss), end;
	std::vector<std::string> tokens(beg, end);

	// remove invalid chromosomes
	return !chromosomeValid(tokens[VcfFile::CHROM]);

}

std::pair<float, float> HaplotypeCallerFilter::getVaf(const std::string&) const {
	return std::pair<float, float>();
}

std::pair<float, float> HaplotypeCallerFilter::getReadCounts(const std::string &) const {
	return std::pair<float, float>();
}
