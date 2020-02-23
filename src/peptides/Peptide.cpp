#include "peptides/Peptide.h"
#include <algorithm>
#include <vector>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <sstream>
#include <sstream>

Peptide Peptide::parse(const std::string &line) {

	std::istringstream iss(line);
	std::string sequence, closestMatch, variantEffectString;
	double immunogenicity;

	iss >> sequence >> variantEffectString >> immunogenicity;

	std::set<std::set<unsigned>> effects;
	std::vector<std::string> variantEffectBlocks;
	boost::split(variantEffectBlocks, variantEffectString, boost::is_any_of("|"));

	std::for_each(variantEffectBlocks.cbegin(), variantEffectBlocks.cend(), [&effects](const std::string& s){

		std::set<unsigned> current;

		std::stringstream ss(s);
		for (int i; ss >> i;) {
			current.emplace(i);
			if (ss.peek() == ',')
				ss.ignore();
		}

		effects.emplace(std::move(current));
	});

	return Peptide(sequence, effects, immunogenicity);
}


Peptide::Peptide(std::string sequence, const std::set<unsigned> &variants)
		: sequence(std::move(sequence)), effectIds() {
	effectIds.emplace(variants);
}

Peptide::Peptide(std::string sequence, std::set<std::set<unsigned>> variants, double immuno)
	: sequence(std::move(sequence)), effectIds(std::move(variants)), immunogenicityScore(immuno){}

void Peptide::addVariantEffects(const std::set<unsigned> &effects){
	effectIds.emplace(effects);
}

void Peptide::combine(const Peptide &other) {
	std::for_each(other.effectIds.cbegin(), other.effectIds.cend(), [this](const std::set<unsigned>& set){
		effectIds.emplace(set);
	});
}

const std::string& Peptide::getSequence() const {
	return sequence;
}

const std::set<std::set<unsigned>>& Peptide::getEffectIds() const {
	return effectIds;
}

void Peptide::setScore(double score) {
	immunogenicityScore = score;
}

double Peptide::getScore() const {
	return immunogenicityScore;
}

std::ostream &operator<<(std::ostream &os, const Peptide &peptide) {
	os << peptide.sequence << '\t';

	for(auto it = peptide.effectIds.cbegin(); it != peptide.effectIds.cend(); it++){
		for(auto it2 = it->cbegin(); it2 != it->cend(); it2++){
			os << *it2;
			if(std::next(it2) != it->cend()) os << ",";
		}
		if(std::next(it) != peptide.effectIds.cend()) os << "|";
	}
	os << '\t' << peptide.immunogenicityScore;

	return os;
}

