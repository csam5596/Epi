#ifndef EPITOPE_PEPTIDE_H
#define EPITOPE_PEPTIDE_H

#include <string>
#include <ostream>
#include <unordered_map>
#include <set>

class Peptide {

public:
	static Peptide parse(const std::string& line);

	Peptide() = default;
	Peptide(std::string sequence, const std::set<unsigned>& variants);

	[[nodiscard]] const std::string& getSequence() const;
	[[nodiscard]] const std::set<std::set<unsigned>>& getEffectIds() const;
	[[nodiscard]] double getScore() const;

	void addVariantEffects(const std::set<unsigned>& effects);
	void combine(const Peptide& other);
	void setScore(double score);

	friend std::ostream &operator<<(std::ostream &os, const Peptide &peptide);

private:
	Peptide(std::string sequence, std::set<std::set<unsigned>> variants, double immuno);

private:
	std::string sequence;
	std::set<std::set<unsigned>> effectIds;
	double immunogenicityScore = -1.;
};


#endif //EPITOPE_PEPTIDE_H
