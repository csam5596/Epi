#include "index/Exon.h"

#include <boost/format.hpp>

Exon::Exon(unsigned number, std::string chr, unsigned start, unsigned end)
	: exonNumber(number), chromosome(std::move(chr)), startPosition(start), endPosition(end){

	location = boost::str(boost::format("chr%1%:%2%-%3%") % chromosome % start % end);
}

bool Exon::operator<(const Exon &rhs) const {
	return exonNumber < rhs.exonNumber;
}