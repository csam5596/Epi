#include <algorithm>
#include <boost/format.hpp>
#include "vcf/VcfHeader.h"

VcfHeader VcfHeader::merge(const std::vector<VcfHeader> &headers) {
	VcfHeader header;

	std::for_each(headers.cbegin(), headers.cend(), [&header](const VcfHeader& h){
		if(header.fileFormat.empty()) header.fileFormat = h.fileFormat;
		header.filters.insert(h.filters.cbegin(), h.filters.cend());
		header.info.insert(h.info.cbegin(), h.info.cend());
		header.format.insert(h.format.cbegin(), h.format.cend());
		header.contigs.insert(h.contigs.cbegin(), h.contigs.cend());
		if(header.columnHeaders.empty()) header.columnHeaders = h.columnHeaders;
	});

	return header;
}

void VcfHeader::addInfo(const std::string &id, const std::string &number, const std::string &type,
                        const std::string &desc) {

	static auto f = boost::format("##INFO=<ID=%1%,Number=%2%,Type=%3%,Description=\"%4%\">");
	info.emplace(id, boost::str(f % id % number % type % desc));
}


std::ostream &operator<<(std::ostream &os, const VcfHeader &header) {

	os << header.fileFormat << std::endl;

	std::for_each(header.filters.cbegin(), header.filters.cend(), [&os](const auto& kv){
		os << kv.second << std::endl;
	});
	std::for_each(header.info.cbegin(), header.info.cend(), [&os](const auto& kv){
		os << kv.second << std::endl;
	});
	std::for_each(header.format.cbegin(), header.format.cend(), [&os](const auto& kv){
		os << kv.second << std::endl;
	});
	std::for_each(header.contigs.cbegin(), header.contigs.cend(), [&os](const auto& kv){
		os << kv.second << std::endl;
	});

	os << header.columnHeaders;

	return os;
}
