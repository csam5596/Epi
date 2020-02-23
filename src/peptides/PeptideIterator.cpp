#include "peptides/PeptideIterator.h"
#include <sstream>
#include <iostream>
PeptideIterator::PeptideIterator(const std::string &file){

	in = std::make_unique<std::ifstream>(file);
	if(in->fail()) throw std::runtime_error("Could not open file " + file);
	if(in->peek() == '#') getline(*in, lineBuffer);
	advance();

}

void PeptideIterator::advance() {
	getline(*in, lineBuffer);
	if(isGood()){
		peptide = Peptide::parse(lineBuffer);
	}
}

bool PeptideIterator::isGood() const {
	return in->good();
}

void PeptideIterator::rewind() {
	in->clear();
	in->seekg(0);
	if(in->peek() == '#') getline(*in, lineBuffer);
	advance();
}

void PeptideIterator::close() {
	in->close();
}