#include <sstream>
#include <iterator>
#include <vector>
#include "vcf/VcfFile.h"
#include "vcf/filters/MutectLineFilter.h"

MutectLineFilter::MutectLineFilter(float maxNormalVaf, float minTumorVaf)
	: VcfLineFilter(maxNormalVaf, minTumorVaf) {}

bool MutectLineFilter::filterVariant(const std::string &line) const {

	std::istringstream ss(line);
	std::istream_iterator<std::string> beg(ss), end;
	std::vector<std::string> tokens(beg, end);

	// remove invalid chromosomes
	if(!chromosomeValid(tokens[VcfFile::CHROM])) return true;

	// if alt allele in normal then ignore
	if(tokens[VcfFile::FILTER].find("alt_allele_in_normal") != std::string::npos) return true;

	std::pair<float,float> vaf = getVaf(line);

	return !vafValid(vaf.first, vaf.second);

}

std::pair<float, float> MutectLineFilter::getVaf(const std::string &line) const {

	std::istringstream ss(line);
	std::istream_iterator<std::string> beg(ss), end;
	std::vector<std::string> tokens(beg, end);

	std::vector<std::string> format;
	boost::split(format, tokens[VcfFile::FORMAT], boost::is_any_of(":"));

	std::vector<std::string> normal;
	boost::split(normal, tokens[VcfFile::NORMAL], boost::is_any_of(":"));

	std::vector<std::string> tumor;
	boost::split(tumor, tokens[VcfFile::TUMOR], boost::is_any_of(":"));

	std::unordered_map<std::string, std::string> normalValues;
	std::unordered_map<std::string, std::string> tumorValues;

	std::transform(format.begin(), format.end(), normal.begin(), std::inserter(normalValues, normalValues.end()),
			[](const std::string& a, const std::string& b) { return std::make_pair(a, b); }
	);

	std::transform(format.begin(), format.end(), tumor.begin(), std::inserter(tumorValues, tumorValues.end()),
			[](const std::string& a, const std::string& b) { return std::make_pair(a, b); }
	);

	float tumorVaf = std::stof(tumorValues.at("AF"));
	float normalVaf = std::stof(normalValues.at("AF"));

	return std::make_pair(normalVaf, tumorVaf);
}

std::pair<float, float> MutectLineFilter::getReadCounts(const std::string &line) const {

	std::istringstream ss(line);
	std::istream_iterator<std::string> beg(ss), end;
	std::vector<std::string> tokens(beg, end);

	std::vector<std::string> format;
	boost::split(format, tokens[VcfFile::FORMAT], boost::is_any_of(":"));

	std::vector<std::string> normal;
	boost::split(normal, tokens[VcfFile::NORMAL], boost::is_any_of(":"));

	std::vector<std::string> tumor;
	boost::split(tumor, tokens[VcfFile::TUMOR], boost::is_any_of(":"));

	std::unordered_map<std::string, std::string> normalValues;
	std::unordered_map<std::string, std::string> tumorValues;

	std::transform(format.begin(), format.end(), normal.begin(), std::inserter(normalValues, normalValues.end()),
	               [](const std::string& a, const std::string& b) { return std::make_pair(a, b); }
	);

	std::transform(format.begin(), format.end(), tumor.begin(), std::inserter(tumorValues, tumorValues.end()),
	               [](const std::string& a, const std::string& b) { return std::make_pair(a, b); }
	);

	float rcNormal = std::stof(normalValues.at("DP"));
	float rcTumor = std::stof(tumorValues.at("DP"));

	return std::make_pair(rcNormal, rcTumor);
}
