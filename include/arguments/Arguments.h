#ifndef EPITOPE_ARGUMENTS_H
#define EPITOPE_ARGUMENTS_H

#include <vector>
#include <string>

struct Arguments {

public:
	Arguments(int argc, char* argv[]);

public:
	std::string dna;
	std::string rna;
	std::vector<std::string> somaticVariants;
	std::vector<std::string> variantCallers;
	std::string germlineVariants;
	std::string proteins;
	std::string gtf;
	std::string ref;
	std::string hla;
	std::string tpm;
	std::string phaser;
	std::string vepCacheDir;

	int minVcCount;
	float minTumorVaf;
	float maxNormalVaf;
	unsigned minReadCountAbs;
	float minReadCountRel;
	unsigned minTpmCount;
	unsigned dnaMapQ;
	unsigned rnaMapQ;

	bool dirty;

private:
	void printValues();

};

#endif //EPITOPE_ARGUMENTS_H
