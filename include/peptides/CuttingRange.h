#ifndef EPITOPE_CUTTINGRANGE_H
#define EPITOPE_CUTTINGRANGE_H

#include <unordered_map>
#include <map>
#include <string>
#include <deque>
#include "InfluenceRange.h"
#include "Peptide.h"
#include <set>

class CuttingRange {

public:
	CuttingRange(unsigned length, std::deque<InfluenceRange>& ir, const std::string& protein);

	std::map<std::string, Peptide> operator()(const std::string& protein);

private:
	unsigned length;
	std::map<std::pair<unsigned, unsigned>, std::set<unsigned>> ranges;
};


#endif //EPITOPE_CUTTINGRANGE_H
