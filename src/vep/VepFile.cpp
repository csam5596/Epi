#include "vep/VepFile.h"
#include "vep/VariantEffect.h"
#include <fstream>
#include <Logging.h>

VepFile VepFile::parse(const std::string &vepFile) {

	LOG_DEBUG("Loading variant effects from " + vepFile + ".");

	std::ifstream in(vepFile);
	if(in.fail()) throw std::runtime_error("Could not open file: " + vepFile);

	std::string line;

	std::unordered_map<unsigned, std::vector<VariantEffect>> effects;
	std::unordered_map<unsigned, VariantEffect> effectsById;
	std::set<std::string> transcriptIds;

	// skip header
	while(getline(in, line)){
		if(line.at(1) != '#') break;
	}

	unsigned long veCount = 0;

	while(getline(in, line)){
		std::optional<VariantEffect> ve = VariantEffect::parse(line);
		if(ve.has_value()){
			transcriptIds.emplace(ve.value().transcript);
			effectsById.emplace(ve.value().variantEffectId, ve.value());
			if(effects.find(ve.value().variantId) == effects.end()){
				effects.emplace(ve.value().variantId, std::vector<VariantEffect>{});
			}
			effects.at(ve.value().variantId).emplace_back(std::move(ve.value()));
			veCount++;
		}
	}

	LOG_DEBUG(std::to_string(veCount) + " variant effects loaded.");

	return VepFile(effects, transcriptIds, effectsById);
}

VepFile::VepFile(std::unordered_map<unsigned, std::vector<VariantEffect>> &effects, std::set<std::string> &ids,
                 std::unordered_map<unsigned, VariantEffect> &effectsById)
		: variantEffects(std::move(effects)), effectsById(std::move(effectsById)), transcriptIds(std::move(ids)){}

const std::set<std::string>& VepFile::getTranscriptIds() const {
	return transcriptIds;
}

const std::optional<std::vector<VariantEffect> const*> VepFile::getVariantEffects(unsigned variantId) const {
	if(variantEffects.find(variantId) == variantEffects.end()) return std::nullopt;
	return &variantEffects.at(variantId);
}

const VariantEffect& VepFile::getVariantEffectById(unsigned id) const {
	return effectsById.at(id);
}