#ifndef EPITOPE_VCFLINEITERATOR_H
#define EPITOPE_VCFLINEITERATOR_H

#include <string>
#include <boost/process.hpp>
#include "VcfLineFilter.h"
#include <vcf/Contig.h>
#include <vcf/VcfHeader.h>
#include <boost/iostreams/filtering_streambuf.hpp>

class VcfLineIterator {

public:
	VcfLineIterator(const std::string& vcf, const std::string& vc, float maxVafNormal, float minVafTumor);

	[[nodiscard]] bool isGood() const;
	void advance();

	bool operator<(const VcfLineIterator &rhs) const;
	bool operator==(const VcfLineIterator &rhs) const;

private:
	void updateComparisonBuffer();
	void parseHeader();
	std::string getHeaderId(const std::string& line) const;


public:
	std::string lineBuffer;
	VcfHeader header;

	float vafNormal = -1.f;
	float vafTumor = -1.f;
	float rcNormal = -1.f;
	float rcTumor = -1.f;

private:
	std::unique_ptr<VcfLineFilter> filter;

	std::unique_ptr<std::ifstream> file;
	std::unique_ptr<boost::iostreams::filtering_streambuf<boost::iostreams::input>> in;
	std::unique_ptr<std::istream> iStream;

	// used for comparison
	Contig contigBuffer;
	unsigned long positionBuffer = 0;
	std::string_view refBuffer;
	std::string_view altBuffer;

};


#endif //EPITOPE_VCFLINEITERATOR_H
