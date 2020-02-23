#ifndef EPITOPE_VARIANT_H
#define EPITOPE_VARIANT_H

#include <string>
#include <vcf/Contig.h>

struct Variant {

public:

	enum VariantType {
		GERMLINE, SOMATIC
	};

	static Variant parse(const std::string& line);

	[[nodiscard]] const std::string& getName() const;

private:
	Variant(Contig &&chr, unsigned position, unsigned id, std::string &ref, std::string &alt, VariantType type,
			int count, float vafNormal, float vafTumor, float rcNormal, float rcTumor);

public:
	Contig contig;
	unsigned position;
	unsigned id;
	std::string ref;
	std::string alt;
	VariantType type;
	std::string name;
	int vcCount;
	float vafNormal;
	float vafTumor;
	float rcNormal;
	float rcTumor;

};


#endif //EPITOPE_VARIANT_H
