#ifndef EPITOPE_VCFHEADER_H
#define EPITOPE_VCFHEADER_H

#include <unordered_map>
#include <map>
#include <vcf/Contig.h>
#include <ostream>
#include <vector>

struct VcfHeader {

public:
	static VcfHeader merge(const std::vector<VcfHeader>& headers);

	VcfHeader() = default;
	void addInfo(const std::string& id, const std::string& number, const std::string& type, const std::string& desc);

public:
	std::string fileFormat;

	std::unordered_map<std::string, std::string> filters;
	std::unordered_map<std::string, std::string> info;
	std::unordered_map<std::string, std::string> format;
	std::map<Contig, std::string> contigs;

	std::string columnHeaders;

	friend std::ostream &operator<<(std::ostream &os, const VcfHeader &header);

};


#endif //EPITOPE_VCFHEADER_H
