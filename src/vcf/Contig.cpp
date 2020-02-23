#include "vcf/Contig.h"
#include <iostream>

Contig::Contig(const std::string_view &chromosome) {

	size_t chrNumberIdx = chromosome.find_first_not_of("chr");
	std::string_view chr = chromosome.substr(chrNumberIdx);

	// if chromosome is a number
	if(chr.find_first_not_of("0123456789") == std::string::npos){
		chrNumber = std::stoul(chr.data());
	}else{
		chrNumber = chr[0];
	}

}

std::string Contig::getName() const {

	std::string chr = "chr";
	if(chrNumber.index() == 1){
		chr += std::get<1>(chrNumber);
	}else{
		chr += std::to_string(std::get<0>(chrNumber));
	}
	return chr;
}

bool Contig::operator<(const Contig &rhs) const {
	// number < X < Y
	return (chrNumber.index() == rhs.chrNumber.index()) ?
		chrNumber < rhs.chrNumber
		: chrNumber.index() < rhs.chrNumber.index();
}

bool Contig::operator==(const Contig &rhs) const {
	return chrNumber == rhs.chrNumber;
}
