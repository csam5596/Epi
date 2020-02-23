#include "index/StartCodon.h"

StartCodon::StartCodon(unsigned exonNumber, std::string chr, unsigned start, unsigned end, char strand)
	: exonNumber(exonNumber), chromosome(std::move(chr)), start(start), end(end), strand(strand == '+' ? 1 : -1){}


std::string StartCodon::toString() const {
	return "chr" + chromosome + ":" + std::to_string(start) + ":" + std::to_string(end);
}
