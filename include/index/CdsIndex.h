#ifndef EPITOPE_CDSINDEX_H
#define EPITOPE_CDSINDEX_H


#include <string>
#include <unordered_map>
#include <Types.h>
#include <set>
#include "Exon.h"
#include "StartCodon.h"

class CdsIndex {

public:
	CdsIndex(const std::string& gtf, const std::string& ref, const std::set<ENST>& transcripts);

	[[nodiscard]] std::optional<std::string> getCodingSequence(const ENST& enst) const;

private:
	void parseGtf(const std::string& file, const std::set<ENST> &transcripts);
	void parsePositions(const std::string& refFile);
	void buildCodingSequences(const std::set<ENST>& transcripts);


	std::unordered_map<ENST, std::string> codingSequences;

	std::unordered_map<ENST, std::set<Exon>> exons;
	std::unordered_map<ENST, StartCodon> startCodons;
	std::unordered_map<std::string, std::string> sequencesByLocation;
};


#endif //EPITOPE_CDSINDEX_H
