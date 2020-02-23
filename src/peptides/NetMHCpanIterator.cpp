#include "peptides/NetMHCpanIterator.h"
#include <iostream>
#include <sstream>

NetMHCpanIterator::NetMHCpanIterator(const std::string &file) {

	in = std::make_unique<std::ifstream>(file);
	if(in->fail()) throw std::runtime_error("Could not open file " + file);

	advance();
}

void NetMHCpanIterator::advance() {

	getline(*in, lineBuffer);
	size_t pos = lineBuffer.find_first_not_of(' ');

	// skip lines until a non empty line with a number at the beginning is found (and not eof reached)
	while(lineBuffer.empty() || pos == std::string::npos || !std::isdigit(lineBuffer.at(pos))){
		getline(*in, lineBuffer);
		if(!isGood()) break;
		pos = lineBuffer.find_first_not_of(' ');
	}

	if(isGood()) parseLine();

}

bool NetMHCpanIterator::isGood() const {
	return in->good();
}

void NetMHCpanIterator::parseLine() {

	std::istringstream iss(lineBuffer);

	// ignore values
	std::string i;
	unsigned ignore = 0;

	iss >> ignore >> hla >> peptide >> core;
	iss >> ignore >> ignore >> ignore >> ignore >> ignore;
	iss >> icore >> i >>score >> affinity >> rank >> exp >> i;

	if(i.empty() || i != "<="){
		bindingLevel = "";
	}else{
		iss >> bindingLevel;
	}

}