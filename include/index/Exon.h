#ifndef EPITOPE_EXON_H
#define EPITOPE_EXON_H

#include <string>

struct Exon {

public:
	Exon(unsigned number, std::string chr ,unsigned start, unsigned end);

	bool operator<(const Exon &rhs) const;

public:
	unsigned exonNumber;
	std::string chromosome;
	unsigned startPosition;
	unsigned endPosition;
	std::string location;
};


#endif //EPITOPE_EXON_H
