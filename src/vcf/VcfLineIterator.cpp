#include <iostream>
#include <regex>
#include "vcf/VcfLineIterator.h"
#include <boost/algorithm/string/find.hpp>
#include <boost/iostreams/filter/gzip.hpp>

VcfLineIterator::VcfLineIterator(const std::string &vcf, const std::string& vc, float maxVafNormal, float minVafTumor)
	: filter(VcfLineFilter::createLineFilter(vc, maxVafNormal, minVafTumor)) {

	file = std::make_unique<std::ifstream>(vcf);
	if(file->fail()) throw std::runtime_error("Could not open file: " + vcf);

	in = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
	if(vcf.substr(vcf.find_last_of('.') + 1) == "gz"){
		in->push(boost::iostreams::gzip_decompressor());
	}
	in->push(*file);

	iStream = std::make_unique<std::istream>(in.get());

	parseHeader();
	advance();
}

void VcfLineIterator::parseHeader() {

	getline(*iStream, lineBuffer);
	header.fileFormat = std::move(lineBuffer);

	while(getline(*iStream, lineBuffer)){
		if(lineBuffer.at(1) != '#'){
			header.columnHeaders = std::move(lineBuffer);
			break;
		}
		if(lineBuffer.find("##contig=") != std::string::npos){
			std::string chromosome = getHeaderId(lineBuffer);
			if(filter->chromosomeValid(chromosome)) {
				header.contigs.emplace(Contig(chromosome), std::move(lineBuffer));
			}
		}else if(lineBuffer.find("##FILTER=") != std::string::npos){
			header.filters.emplace(getHeaderId(lineBuffer), std::move(lineBuffer));
		}else if(lineBuffer.find("##INFO=") != std::string::npos){
			header.info.emplace(getHeaderId(lineBuffer), std::move(lineBuffer));
		}else if(lineBuffer.find("##FORMAT=") != std::string::npos){
			header.format.emplace(getHeaderId(lineBuffer), std::move(lineBuffer));
		}
	}

}

bool VcfLineIterator::isGood() const {
	return iStream->good();
}

void VcfLineIterator::advance() {
	while(getline(*iStream, lineBuffer) && filter->filterVariant(lineBuffer));
	if(isGood()) updateComparisonBuffer();
}


void VcfLineIterator::updateComparisonBuffer() {

	std::vector<int> tabs;
	tabs.reserve(5);
	size_t pos = lineBuffer.find('\t');
	tabs.emplace_back(pos);

	while( pos != std::string::npos && tabs.size() != 5){
		pos =lineBuffer.find('\t', pos + 1);
		tabs.emplace_back(pos);
	}

	contigBuffer = Contig(std::string_view(lineBuffer).substr(0, tabs[0]));
	positionBuffer = std::stoul(std::string_view(lineBuffer).substr(tabs[0]+1, tabs[1]-(tabs[0]+1)).data());
	refBuffer = std::string_view(lineBuffer).substr(tabs[2]+1, tabs[3]-(tabs[2]+1));
	altBuffer = std::string_view(lineBuffer).substr(tabs[3]+1, tabs[4]-(tabs[3]+1));

	std::pair<float, float> vaf = filter->getVaf(lineBuffer);
	vafNormal = vaf.first;
	vafTumor = vaf.second;

	std::pair<float, float> readCounts = filter->getReadCounts(lineBuffer);
	rcNormal = readCounts.first;
	rcTumor = readCounts.second;

}

bool VcfLineIterator::operator<(const VcfLineIterator &rhs) const {
	return std::tie(contigBuffer, positionBuffer, refBuffer, altBuffer) <
	       std::tie(rhs.contigBuffer, rhs.positionBuffer, rhs.refBuffer, rhs.altBuffer);
}

bool VcfLineIterator::operator==(const VcfLineIterator &rhs) const {
	return std::tie(contigBuffer, positionBuffer, refBuffer, altBuffer) ==
	       std::tie(rhs.contigBuffer, rhs.positionBuffer, rhs.refBuffer, rhs.altBuffer);
}

std::string VcfLineIterator::getHeaderId(const std::string &line) const {
	static const std::regex idRegex = std::regex("ID=([^,])+,");

	size_t start = line.find("ID=");
	size_t end = line.find_first_of(",>", start+3);
	return line.substr(start+3, end-(start+3));
}

