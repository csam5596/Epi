#include "phasing/Haplotype.h"

#include <fstream>
#include <iostream>
#include <Logging.h>
#include <boost/format.hpp>

Haplotype::Haplotype(const std::string &filePrefix, const std::string &rna){

	if(!rna.empty()){
		parseRnaCounts(filePrefix, rna);
	}

	std::ifstream in(filePrefix + ".haplotypes.txt");
	if(in.fail()) throw std::runtime_error("Could not open file: " + filePrefix + ".haplotypes.txt");

	std::string line;

	// ignore header
	getline(in, line);

	while(getline(in, line)){
		haploblocks.emplace_back(Haploblock::parse(line));
	}
}

void Haplotype::parseRnaCounts(const std::string &file, const std::string &rna) {

	if(rna.empty()) return;

	std::ifstream in(file + ".haplotypic_counts.txt");
	if(in.fail()) throw std::runtime_error("Could not open " + file + ".haplotypic_counts.txt");

	unsigned phaseA;
	unsigned phaseB;
	std::string variants;

	std::string line;

	size_t slashIndex = rna.find_last_of('/');
	size_t dotIndex = rna.find_last_of('.');
	std::string_view rnaFileName = std::string_view(rna).substr(slashIndex+1, dotIndex-slashIndex-1);

	getline(in, line); //ignore header

	while(getline(in, line)){
		if(line.find(rnaFileName) != std::string::npos){
			std::istringstream iss(line);
			iss.ignore(9999, '\t'); // contig
			iss.ignore(9999, '\t'); // start
			iss.ignore(9999, '\t'); // end
			iss >> variants;
			iss.ignore(9999, '\t'); // tab
			iss.ignore(9999, '\t'); // variant count
			iss.ignore(9999, '\t'); // blacklisted
			iss.ignore(9999, '\t'); // blacklisted count
			iss.ignore(9999, '\t'); // haplotype A
			iss.ignore(9999, '\t'); // haplotype B
			iss >> phaseA >> phaseB;

			std::vector<std::string> variantNames;
			boost::split(variantNames, variants, boost::is_any_of(","));

			std::for_each(variantNames.begin(), variantNames.end(), [this, phaseA, phaseB](std::string& name){
				rnaPhaseCounts.emplace(std::move(name), std::make_pair(phaseA, phaseB));
			});

		}
	}
}

void Haplotype::computeProteinChanges(const VcfFile &vcf, const VepFile &vep, unsigned minRcAbs, float minRcRel){

	for(const Haploblock& block : haploblocks){

		std::unordered_map<ENSP, std::set<CodonChange>> phaseA;
		std::unordered_map<ENSP, std::set<CodonChange>> phaseB;

		long totalCount = block.countA + block.countB;
		if(totalCount == 0) continue; // should not happen

		for(unsigned long index = 0; index < block.variantIds.size(); index++){
			unsigned variantId = block.variantIds[index];

			const Variant& variant = vcf.getVariantById(variantId);
			const std::optional<std::vector<VariantEffect> const*> effects = vep.getVariantEffects(variantId);

			if(!effects.has_value()) continue;
			std::string variantName = variant.getName();
			std::replace( variantName.begin(), variantName.end(), ':', '_');

			const std::pair<unsigned, unsigned>& rnaCounts = [this,&variantName]{
				if(rnaPhaseCounts.find(variantName) != rnaPhaseCounts.end()){
					return rnaPhaseCounts.at(variantName);
				}else{
					LOG_WARN("No rna counts found for " + variantName);
					return std::make_pair(0u,0u);
				}
			}();
			float rnaCountsSum = static_cast<float>(std::max(rnaCounts.first + rnaCounts.second, 1u));

			bool isPhaseA = block.allelesA[index] == variant.alt;
			unsigned readCount = (isPhaseA) ? block.countA : block.countB;
			float rcRel = static_cast<float>(readCount) / static_cast<float>(totalCount);

			if(readCount < minRcAbs || rcRel < minRcRel){
				static auto reasonAbs = boost::format("Absolute Read Count: %1% < %2%");
				static auto reasonRel = boost::format("Relative Read Count: %1% < %2%");
				std::string reason = (readCount < minRcAbs) ? boost::str(reasonAbs % readCount % minRcAbs) :
						boost::str(reasonRel % rcRel % minRcRel);

				LOG_WARN("Dropped variant " + std::to_string(variantId) + "\t" + variant.getName() + " - " + reason);
				continue;
			}

			std::unordered_map<ENSP, std::set<CodonChange>>* phase = (isPhaseA) ? &phaseA : &phaseB;

			float rnaVaf = isPhaseA ? (rnaCounts.first / rnaCountsSum) : (rnaCounts.second / rnaCountsSum);
			vafRNAByVariant.emplace(variantId, rnaVaf);
			rcRNAByVariant.emplace(variantId, rnaCountsSum);

			for(const VariantEffect& ve : *effects.value()){
				if(phase->find(ve.proteinId) == phase->end()){
					phase->emplace(ve.proteinId, std::set<CodonChange>{});
				}
				phase->at(ve.proteinId).emplace(variant, ve);
			}
		}

		// merge
		std::for_each(phaseA.begin(), phaseA.end(), [this](auto& kv){
			if(changesByProtein.find(kv.first) == changesByProtein.end()){
				changesByProtein.emplace(kv.first, std::vector<std::set<CodonChange>>{});
			}
			changesByProtein.at(kv.first).emplace_back(std::move(kv.second));
		});

		std::for_each(phaseB.begin(), phaseB.end(), [this](auto& kv){
			if(changesByProtein.find(kv.first) == changesByProtein.end()){
				changesByProtein.emplace(kv.first, std::vector<std::set<CodonChange>>{});
			}
			changesByProtein.at(kv.first).emplace_back(std::move(kv.second));
		});
	}

}


void Haplotype::addHomozygousVariants(const std::string &filename, const VepFile &vep) {

	LOG_DEBUG("Adding Homozygous Variants from " + filename);

	VcfFile vcf = VcfFile::parse(filename);
	const std::unordered_map<unsigned long, Variant>& variants = vcf.getVariants();

	std::for_each(variants.cbegin(), variants.cend(), [this, &vep](const auto& kv){
		const std::optional<std::vector<VariantEffect> const*> effects = vep.getVariantEffects(kv.first);

		if(effects.has_value()){

			for(const VariantEffect& ve : *effects.value()){
				if(changesByProtein.find(ve.proteinId) == changesByProtein.end()){

					std::vector<std::set<CodonChange>> changeVector;
					std::set<CodonChange> changes{};
					changes.emplace(kv.second, ve);
					changeVector.emplace_back(std::move(changes));

					changesByProtein.emplace(ve.proteinId, std::move(changeVector));

				}else{

					std::vector<std::set<CodonChange>>& proteinChanges = changesByProtein.at(ve.proteinId);
					std::for_each(proteinChanges.begin(), proteinChanges.end(), [&](std::set<CodonChange>& changes){
						changes.emplace(kv.second, ve);
					});

				}
			}
		}

	});

}

const std::unordered_map<ENSP, std::vector<std::set<CodonChange>>>& Haplotype::getProteinChanges() const {
	return changesByProtein;
}

std::optional<float> Haplotype::getRnaVaf(unsigned long variantId) const {
	if(vafRNAByVariant.find(variantId) != vafRNAByVariant.end()){
		return vafRNAByVariant.at(variantId);
	}
	return std::nullopt;
}

std::optional<unsigned> Haplotype::getRcRNA(unsigned long variantId) const {
	if(rcRNAByVariant.find(variantId) != rcRNAByVariant.end()){
		return rcRNAByVariant.at(variantId);
	}
	return std::nullopt;
}