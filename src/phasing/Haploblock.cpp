#include <sstream>
#include "phasing/Haploblock.h"

#include <vector>
#include <iterator>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>

Haploblock Haploblock::parse(const std::string &line) {

	std::istringstream is(line);

	is.ignore(9999, '\t'); // contig
	is.ignore(9999, '\t'); // start
	is.ignore(9999, '\t'); // stop
	is.ignore(9999, '\t'); // length

	unsigned variantCount;
	std::string variantIds;
	std::string alleles;
	unsigned countA;
	unsigned countB;

	is >> variantCount >> variantIds >> alleles >> countA >> countB;

	// parse variant ids
	std::replace(variantIds.begin(), variantIds.end(), ',', ' ');  // replace ',' by ' '
	std::istringstream iss(variantIds);
	std::vector<unsigned> ids(std::istream_iterator<unsigned>{iss}, std::istream_iterator<unsigned>());


	// parse both allele configurations
	size_t splitIndex = alleles.find('|');
	std::vector<std::string> allelesA;
	std::vector<std::string> allelesB;

	allelesA.reserve(variantCount);
	allelesB.reserve(variantCount);

	std::string phaseA = alleles.substr(0, splitIndex);
	boost::split(allelesA, phaseA, boost::is_any_of(","));

	std::string phaseB = alleles.substr(splitIndex+1);
	boost::split(allelesB, phaseB, boost::is_any_of(","));

	return Haploblock(ids, allelesA, allelesB, countA, countB);
}

Haploblock::Haploblock(std::vector<unsigned int> &variantIds, std::vector<std::string> &allelesA,
                       std::vector<std::string> &allelesB, unsigned countA, unsigned countB)
       : variantIds(std::move(variantIds)), allelesA(std::move(allelesA)),
         allelesB(std::move(allelesB)), countA(countA), countB(countB) {}
