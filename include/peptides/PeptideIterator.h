#ifndef EPITOPE_PEPTIDEITERATOR_H
#define EPITOPE_PEPTIDEITERATOR_H

#include <string>
#include <memory>
#include <fstream>
#include "Peptide.h"

class PeptideIterator {

public:
	explicit PeptideIterator(const std::string& file);

	void advance();
	void rewind();
	void close();

	[[nodiscard]] bool isGood() const;

public:
	Peptide peptide;

private:
	std::unique_ptr<std::ifstream> in;
	std::string lineBuffer;


};


#endif //EPITOPE_PEPTIDEITERATOR_H
