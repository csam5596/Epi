#include "index/Protein.h"

Protein Protein::parse(const std::string &line) {

	std::string proteinId = line.substr(1, line.find_first_of(". ")-1);
	size_t transcriptStart = line.find("transcript:");
	size_t transcriptEnd = line.find_first_of(". ", transcriptStart);

	std::string transcriptId = line.substr(transcriptStart + 11, transcriptEnd - transcriptStart - 11);

	size_t geneSymbolStart = line.find("gene_symbol:");
	size_t geneSymbolEnd = line.find_first_of(". ", geneSymbolStart);
	std::string geneSymbol = line.substr(geneSymbolStart + 12, geneSymbolEnd - geneSymbolStart - 12);

	return Protein(proteinId, transcriptId, geneSymbol);
}

Protein::Protein(ENSP &proteinId, ENST &transcriptId, std::string& geneSymbol)
	: proteinId(std::move(proteinId)), transcriptId(std::move(transcriptId)), geneSymbol(std::move(geneSymbol)){}
