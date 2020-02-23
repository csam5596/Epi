#ifndef EPITOPE_RESULTS_H
#define EPITOPE_RESULTS_H

#include <string>
#include <peptides/NetMHCpanIterator.h>
#include <fstream>
#include <peptides/PeptideIterator.h>

namespace Results {

	static inline void generate(const std::string& netMhc, const std::string& peptides, const VcfFile& vcf,
							    const VepFile& vep, const Haplotype& haplotype, const TpmCount& tpm,
							    const ProteinIndex& proteinIndex, unsigned minTpmCount, const std::string& out){

		static const std::string header =
				"Peptide\tVariants\tVariantTypes\tVariantCallers\t"
				"DNA VAF Normal\tDNA VAF Tumor\tRNA VAF\t"
				"Read Count DNA Normal\tRead Count DNA Tumor\tRead Count RNA Tumor\tgene\tProtein\t"
				"HLA\tBindingCore\tICore\tRawPredictionScore\tAffinity(nM)\t%Rank\tExp\tStrength\t"
		        "TPM\texpressed\tImmunogenicityScore\tImmunogenic";


		LOG_DEBUG("Printing results to " + out);

		NetMHCpanIterator iterator(netMhc);
		PeptideIterator peptideIterator(peptides);

		std::ofstream o(out);
		if(o.fail()) throw std::runtime_error("Could not open file " + out);

		o << header << std::endl;

		while(iterator.isGood()){

			bool valid = true;

			if(iterator.peptide != peptideIterator.peptide.getSequence()){
				LOG_ERROR("Peptides do not match " + iterator.peptide + "\t" + peptideIterator.peptide.getSequence());
				valid = false;
			}


			const std::set<std::set<unsigned>>& effectGroups = peptideIterator.peptide.getEffectIds();
			std::for_each(effectGroups.cbegin(), effectGroups.cend(), [&](const std::set<unsigned>& ids){

				std::string variants;
				std::string variantCallers;
				std::string variantTypes;
				std::string normalVaf;
				std::string tumorVaf;
				std::string rnaVaf;
				std::string rcDnaN;
				std::string rcDnaT;
				std::string rcRNA;
				std::string protein;
				std::string transcript;

				std::for_each(ids.cbegin(), ids.cend(), [&](const unsigned& effectId){
					const VariantEffect& ve = vep.getVariantEffectById(effectId);
					const Variant& v = vcf.getVariantById(ve.variantId);

					if(!variants.empty()) variants += ",";
					if(!variantCallers.empty()) variantCallers += ",";
					if(!variantTypes.empty()) variantTypes += ",";
					if(!normalVaf.empty()) normalVaf += ",";
					if(!tumorVaf.empty()) tumorVaf += ",";
					if(!rnaVaf.empty()) rnaVaf += ",";
					if(!rcRNA.empty()) rcRNA += ",";
					if(!rcDnaN.empty()) rcDnaN += ",";
					if(!rcDnaT.empty()) rcDnaT += ",";

					variants += v.getName();
					variantCallers += std::to_string(v.vcCount);
					variantTypes += ve.getTypeAsString();
					normalVaf += std::to_string(v.vafNormal);
					tumorVaf += std::to_string(v.vafTumor);
					rnaVaf += std::to_string(haplotype.getRnaVaf(ve.variantId).value_or(0.f));
					rcRNA += std::to_string(haplotype.getRcRNA(ve.variantId).value_or(0));
					rcDnaN += std::to_string(v.rcNormal);
					rcDnaT += std::to_string(v.rcTumor);

					if(protein.empty()) {
						protein = ve.proteinId;
					}else if(protein != ve.proteinId){
						LOG_ERROR("Protein Ids do not match: [Expected]: " + protein + " [Actual]: " + ve.proteinId);
						valid = false;
					}
					if(transcript.empty()) {
						transcript = ve.transcript;
					}else if(transcript != ve.transcript){
						LOG_ERROR("Transcripts do not match: [Expected]: " + transcript + " [Actual]: " + ve.transcript);
						valid = false;
					}

				});

				if(valid){
					std::optional<Protein const*> affectedProtein = proteinIndex.getProtein(protein);
					std::string gene = affectedProtein.has_value() ? affectedProtein.value()->geneSymbol : "";
					const std::string& proteinId = affectedProtein.has_value() ? affectedProtein.value()->proteinId : "";

					const Peptide& peptide = peptideIterator.peptide;

					o << peptide.getSequence() << '\t' << variants << '\t' << variantTypes << '\t' << variantCallers << '\t';
					o << normalVaf << '\t' << tumorVaf << '\t' << rnaVaf << '\t';
					o << rcDnaN << '\t' << rcDnaT << '\t' << rcRNA << '\t' << gene << '\t' << proteinId << '\t';

					o << iterator.hla << '\t' << iterator.core << "\t" << iterator.icore << '\t' << iterator.score << '\t';
					o << iterator.affinity << '\t' << iterator.rank << '\t' << iterator.exp << '\t' << iterator.bindingLevel << '\t';

					float tpmCount = tpm.getTpm(transcript).value_or(0);
					char expressed = tpmCount > static_cast<float>(minTpmCount) ? 'Y' : 'N';
					o << std::to_string(tpmCount) << '\t' << expressed << '\t';
					o << peptide.getScore() << '\t' << (peptide.getScore() > 0.5 ? "yes" : "no") << std::endl;
				}
			});

			iterator.advance();
			peptideIterator.advance();

			if(!peptideIterator.isGood()) peptideIterator.rewind();
		}

		LOG_DEBUG("Finished");
	}

}

#endif //EPITOPE_RESULTS_H
