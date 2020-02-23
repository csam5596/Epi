#ifndef EPITOPE_PROTEIN_H
#define EPITOPE_PROTEIN_H

#include <string>
#include <Types.h>

class Protein {

public:
	static Protein parse(const std::string& line);

private:
	Protein(ENSP &proteinId, ENST &transcriptId, std::string& geneSymbol);

public:
	ENSP proteinId;
	ENST transcriptId;
	std::string geneSymbol;
	std::string aminoAcids;

};


#endif //EPITOPE_PROTEIN_H
