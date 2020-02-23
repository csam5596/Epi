#include "peptides/LineIterator.h"
#include <fstream>

LineIterator::LineIterator(const std::string &file) {
	in = std::make_unique<std::ifstream>(file);
	if(in->fail()) throw std::runtime_error("Could not open file " + file);

	advance();
}

void LineIterator::advance() {
	getline(*in, lineBuffer);
}

bool LineIterator::isGood() const {
	return in->good();
}

void LineIterator::close() {
	in->close();
}