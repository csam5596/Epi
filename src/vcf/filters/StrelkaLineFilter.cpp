#include <sstream>
#include <iterator>
#include <vector>
#include <vcf/VcfFile.h>
#include "vcf/filters/StrelkaLineFilter.h"

StrelkaLineFilter::StrelkaLineFilter(float maxNormalVaf, float minTumorVaf)
	: VcfLineFilter(maxNormalVaf, minTumorVaf) {}

bool StrelkaLineFilter::filterVariant(const std::string &line) const {

	std::istringstream ss(line);
	std::istream_iterator<std::string> beg(ss), end;
	std::vector<std::string> tokens(beg, end);

	// remove invalid chromosomes
	if(!chromosomeValid(tokens[VcfFile::CHROM])) return true;

	std::pair<float, float> vaf = getVaf(line);
	return !vafValid(vaf.first, vaf.second);
}

/**
 ** computed according to:
 ** https://github.com/Illumina/strelka/tree/master/docs/userGuide#somatic
 **/
std::pair<float, float> StrelkaLineFilter::getVaf(const std::string &line) const {

	std::istringstream ss(line);
	std::istream_iterator<std::string> beg(ss), end;
	std::vector<std::string> tokens(beg, end);

	std::string key = std::string(tokens[VcfFile::REF] + "U");

	std::vector<std::string> format;
	boost::split(format, tokens[VcfFile::FORMAT], boost::is_any_of(":"));
	bool isSNV = std::find(format.begin(), format.end(), key) != format.end();

	const std::string refKey = isSNV ? (tokens[VcfFile::REF] + "U") : "TAR";
	const std::string altKey = isSNV ? (tokens[VcfFile::ALT] + "U") : "TIR";

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

	const std::string& normalRef = normalValues.at(refKey);
	const std::string& tumorRef = tumorValues.at(refKey);

	const std::string& normalAlt = normalValues.at(altKey);
	const std::string& tumorAlt = tumorValues.at(altKey);

	float tier1RefCountsNormal = std::stof(normalRef.substr(0, normalRef.find(',')));
	float tier1AltCountsNormal = std::stof(normalAlt.substr(0, normalAlt.find(',')));


	float tier1RefCountsTumor = std::stof(tumorRef.substr(0, tumorRef.find(',')));
	float tier1AltCountsTumor = std::stof(tumorAlt.substr(0, tumorAlt.find(',')));

	return std::make_pair(
			(tier1AltCountsNormal / (tier1AltCountsNormal + tier1RefCountsNormal)),
			(tier1AltCountsTumor / (tier1AltCountsTumor + tier1RefCountsTumor))
	);
}

std::pair<float, float> StrelkaLineFilter::getReadCounts(const std::string &line) const {
	std::istringstream ss(line);
	std::istream_iterator<std::string> beg(ss), end;
	std::vector<std::string> tokens(beg, end);

	std::string key = std::string(tokens[VcfFile::REF] + "U");

	std::vector<std::string> format;
	boost::split(format, tokens[VcfFile::FORMAT], boost::is_any_of(":"));
	bool isSNV = std::find(format.begin(), format.end(), key) != format.end();

	const std::string refKey = isSNV ? (tokens[VcfFile::REF] + "U") : "TAR";
	const std::string altKey = isSNV ? (tokens[VcfFile::ALT] + "U") : "TIR";

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

	const std::string& normalRef = normalValues.at(refKey);
	const std::string& tumorRef = tumorValues.at(refKey);

	const std::string& normalAlt = normalValues.at(altKey);
	const std::string& tumorAlt = tumorValues.at(altKey);

	float tier1RefCountsNormal = std::stof(normalRef.substr(0, normalRef.find(',')));
	float tier1AltCountsNormal = std::stof(normalAlt.substr(0, normalAlt.find(',')));

	float tier1RefCountsTumor = std::stof(tumorRef.substr(0, tumorRef.find(',')));
	float tier1AltCountsTumor = std::stof(tumorAlt.substr(0, tumorAlt.find(',')));

	float rcNormal = tier1RefCountsNormal + tier1AltCountsNormal;
	float rcTumor = tier1RefCountsTumor + tier1AltCountsTumor;

	return std::make_pair(rcNormal, rcTumor);
}
