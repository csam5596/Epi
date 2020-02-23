#include "vep/VariantEffect.h"
#include <vector>
#include <boost/algorithm/string.hpp>
#include <vep/VepFile.h>
#include <util/VepUtil.h>

VariantEffect::VariantEffectType VariantEffect::getType(const std::string &refCodon, const std::string &altCodon) {

	unsigned refLength = refCodon.length();
	unsigned altLength = altCodon.length();

	if(refLength == altLength) return VariantEffectType::SNV;
	if(refLength % 3 != 0 || altLength % 3 != 0) return VariantEffectType::FRAMESHIFT_VARIANT;
	return (refLength < altLength) ? VariantEffectType::INFRAME_INSERTION : VariantEffectType::INFRAME_DELETION;

}

std::optional<VariantEffect> VariantEffect::parse(const std::string &line) {

	std::vector<std::string> tokens;
	boost::split(tokens, line, boost::is_any_of("\t"));

	if(VepUtil::filter(tokens)) return std::nullopt;

	unsigned id = std::stoul(tokens[VepFile::UPLOADED_VARIATION]);

	// parse location
	auto location = std::string_view(tokens[VepFile::LOCATION]);
	size_t colonPos = location.find(':');
	size_t positionEnd = location.find(colonPos, '-');

	Contig contig(location.substr(0, colonPos));
	unsigned position = std::stoul(location.substr(colonPos+1, positionEnd-colonPos-1).data());


	std::string& gene = tokens[VepFile::GENE];
	std::string& transcript = tokens[VepFile::FEATURE];


	auto proteinPosStr = std::string_view(tokens[VepFile::PROTEIN_POSITION]);
	unsigned proteinPosition = std::stoul(proteinPosStr.substr(0, proteinPosStr.find('-')).data());


	// Parse extras
	auto extras = std::string_view(tokens[VepFile::EXTRA]);
	size_t enspStart = extras.find("ENSP=");
	size_t enspEnd = extras.find(enspStart, ';');

	size_t strandStart = extras.find("STRAND=");
	size_t strandEnd = extras.find(strandStart, ';');

	std::string ensp = extras.substr(enspStart+5, enspEnd-enspStart-5).data();
	int strand = std::stoi(extras.substr(strandStart+7, strandEnd - strandStart - 7).data());


	// parse codons
	boost::erase_all(tokens[VepFile::CODONS], "-");
	auto codons = std::string_view(tokens[VepFile::CODONS]);
	size_t split = codons.find('/');
	auto ref = std::string(codons.substr(0, split));
	auto alt = std::string(codons.substr(split+1));

	return VariantEffect(id, contig, position, gene, transcript, strand, ensp, ref, alt, proteinPosition);

}

VariantEffect::VariantEffect(unsigned id, Contig &chr, unsigned pos, std::string &gene, std::string &transcript,
                             int strand, std::string &ensp, std::string &cRef, std::string &cAlt, unsigned proteinPos)
     : variantId(id), chromosome(std::move(chr)), position(pos), gene(std::move(gene)),
       transcript(std::move(transcript)), strand(strand), proteinId(std::move(ensp)), codonRef(std::move(cRef)),
       codonAlt(std::move(cAlt)), proteinPosition(proteinPos){

	static unsigned counter = 0;
	type = getType(codonRef, codonAlt);
	variantEffectId = counter++;

}

const std::string& VariantEffect::getTypeAsString() const {
	static const std::string snv = "SNV";
	static const std::string ins = "Insertion";
	static const std::string del = "Deletion";
	static const std::string fs = "Frameshift";
	static const std::string empty;

	switch(type){
		case SNV:
			return snv;
		case INFRAME_INSERTION:
			return ins;
		case INFRAME_DELETION:
			return del;
		case FRAMESHIFT_VARIANT:
			return fs;
		default:
			return empty;
	}
}