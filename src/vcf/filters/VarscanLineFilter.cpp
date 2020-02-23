#include <vector>
#include <sstream>
#include <iterator>
#include <vcf/VcfFile.h>
#include "vcf/filters/VarscanLineFilter.h"

VarscanLineFilter::VarscanLineFilter(float maxNormalVaf, float minTumorVaf)
	: VcfLineFilter(maxNormalVaf, minTumorVaf) {}

bool VarscanLineFilter::filterVariant(const std::string &line) const {

	std::istringstream ss(line);
	std::istream_iterator<std::string> beg(ss), end;
	std::vector<std::string> tokens(beg, end);

	// remove invalid chromosomes
	if(!chromosomeValid(tokens[VcfFile::CHROM])) return true;

	// check info
	std::unordered_map<std::string, std::string> info;
	std::istringstream infoStream(tokens[VcfFile::INFO]);
	std::string infoToken;

	while(getline(infoStream, infoToken, ';')){
		size_t pos = infoToken.find('=');
		info.emplace(infoToken.substr(0,pos), infoToken.substr(pos+1));
	}

	// Varscan: <ID=SS,Number=1,Type=String,Description="Somatic status of variant (0=Reference,1=Germline,2=Somatic,3=LOH, or 5=Unknown)"
	if(info.at("SS") != "2") return true;

	// check vaf
	std::pair<float,float> vaf = getVaf(line);
	return !vafValid(vaf.first, vaf.second);
}

std::pair<float, float> VarscanLineFilter::getVaf(const std::string &line) const {

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

	const std::string& normalVafString = normalValues.at("FREQ");
	const std::string& tumorVafString = tumorValues.at("FREQ");

	float normalVaf = std::stof(normalVafString.substr(0, normalVafString.find('%'))) / 100.f;
	float tumorVaf = std::stof(tumorVafString.substr(0, tumorVafString.find('%'))) / 100.f;

	return std::make_pair(normalVaf, tumorVaf);

}

std::pair<float, float> VarscanLineFilter::getReadCounts(const std::string &line) const {

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


	float rcNormal = std::stof(normalValues.at("RD")) + std::stof(normalValues.at("AD"));
	float rcTumor = std::stof(tumorValues.at("RD")) + std::stof(tumorValues.at("AD"));
	return std::make_pair(rcNormal, rcTumor);

}
