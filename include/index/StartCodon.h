#ifndef EPITOPE_STARTCODON_H
#define EPITOPE_STARTCODON_H

#include <string>
#include <ostream>

struct StartCodon {

public:
	StartCodon(unsigned exonNumber, std::string chr, unsigned start, unsigned end, char strand);

	[[nodiscard]] std::string toString() const;

public:
	unsigned exonNumber;
	std::string chromosome;
	unsigned start;
	unsigned end;
	int strand;
};


#endif //EPITOPE_STARTCODON_H
