#include <fstream>
#include <boost/algorithm/string/join.hpp>
#include <sstream>
#include "hla/HlaTypes.h"
#include <iostream>
#include <Logging.h>

HlaTypes HlaTypes::parse(const std::string &file) {


	std::ifstream in(file);
	if(in.fail()) throw std::runtime_error("Could not open file: " + file);

	LOG_DEBUG("Parsing HLA-Types from " + file);

	const std::string hlaPrefix = "HLA-";
	std::set<std::string> types;

	std::string line;
	getline(in, line);	// ignore header
	getline(in, line);
	std::istringstream iss(line);

	iss.ignore(9999, '\t'); // ignore line number

	for(int i = 0; i<6; i++){
		iss >> line;
		line.erase (std::remove(line.begin(), line.end(), '*'), line.end());
		types.emplace(hlaPrefix + line);
	}

	return HlaTypes(types);
}

std::string HlaTypes::toString() const {
	return boost::join(types, ",");
}

HlaTypes::HlaTypes(std::set<std::string> &types)
	: types(std::move(types)){}