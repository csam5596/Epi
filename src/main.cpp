#include <vep/VepFile.h>
#include <phasing/Haplotype.h>
#include <peptides/NetMhcPan.h>
#include <hla/HlaTypes.h>
#include <Logging.h>
#include <peptides/Peptidome.h>
#include <tpm/TpmCount.h>
#include <Results.h>
#include <index/ProteinIndex.h>
#include "vcf/VcfFile.h"
#include "util/VepUtil.h"
#include "arguments/Arguments.h"
#include "phasing/Phaser.h"
#include "util/VcfUtil.h"
#include <experimental/filesystem>
#include <Prerequisites.h>
#include <immunogenicity/RandomForest.h>
#include <immunogenicity/BlomapEncoder.h>

#define LOAD_RESOURCE(x) ([]() {                                            \
        extern const char _binary_##x##_start, _binary_##x##_end;           \
        return std::istringstream(std::string(&_binary_##x##_start, &_binary_##x##_end - &_binary_##x##_start)); \
    })()

int main(int argc, char* argv[]) {

	Arguments args(argc, argv);

	/** SET UP **/
	mkdir("tmp", 0777);
	mkdir("log", 0777);
	Logging::init();

	/** CHECK PREREQUISITES **/
	if(!std::experimental::filesystem::exists(args.dna + ".bai")){
		LOG_ERROR("No index file found for " + args.dna + "\t[expected: " + args.dna + ".bai]");
		exit(1);
	}

	if(!args.rna.empty() && ! std::experimental::filesystem::exists(args.rna + ".bai")){
		LOG_ERROR("No index file found for " + args.rna + "\t[expected: " + args.rna + ".bai]");
		exit(1);
	}

	if(!Prerequisites::check()) exit(1);


	/** COMBINE AND FILTER VARIANTS **/
	VcfUtil::filterAndCombine(args.somaticVariants, args.variantCallers, args.minVcCount, args.maxNormalVaf, args.minTumorVaf, "tmp/somatic.combined.vcf");
	VcfUtil::filter(args.germlineVariants, "HaplotypeCaller", 0.0f, 0.0f,"tmp/germline.filtered.vcf");
	VcfUtil::difference("tmp/somatic.combined.vcf", "tmp/germline.filtered.vcf", "tmp/somatic.filtered.vcf");
	VcfUtil::mergeGermlineSomatic("tmp/germline.filtered.vcf", "tmp/somatic.filtered.vcf", "tmp/combined.vcf");
	VcfFile vcf = VcfFile::parse("tmp/combined.vcf");

	/** CALL VEP **/
	VepUtil::predictEffects(args.vepCacheDir, "tmp/combined.vcf", "tmp/effects.vep.txt");
	VepFile vepFile = VepFile::parse("tmp/effects.vep.txt");

	/** PHASE VARIANTS **/
	VcfUtil::splitZygosity("tmp/combined.vcf", "tmp/combined");

	mkdir("tmp/phaser", 0777);
	Phaser::phaseVariants(args.phaser, "tmp/combined.heterozygous.vcf", args.dna, args.rna, args.dnaMapQ, args.rnaMapQ, "tmp/phaser/phaser");
	Haplotype haplotype("tmp/phaser/phaser", args.rna);
	haplotype.computeProteinChanges(vcf, vepFile, args.minReadCountAbs, args.minReadCountRel);
	haplotype.addHomozygousVariants("tmp/combined.homozygous.vcf", vepFile);

	/** BUILD PROTEIN INDEX **/
	ProteinIndex proteinIndex(args.proteins);

	/** BUILD PEPTIDOME **/
	Peptidome peptidome;
	peptidome.createPeptides(args, haplotype, vepFile, proteinIndex);
	peptidome.filterPeptides();
	peptidome.mergeNeoPeptides("tmp/neo.pep");

	/** IMMUNOGENICITY PREDICTION **/
	auto rfIss = LOAD_RESOURCE(resources_rfBm_txt_gz);
	RandomForest rf(rfIss, 500);
	rf.predictAll("tmp/neo.pep", BlomapEncoder::encode);

	/** PARSE HLA **/
	HlaTypes hla = HlaTypes::parse(args.hla);

	/** NET MHC PAN **/
	NetMhcPan::call( "tmp/neo.pep", hla);

	/** PARSE KALLISTO **/
	TpmCount tpm = TpmCount::parse(args.tpm);

	/** PRINT RESULTS **/
	Results::generate("log/netMHCpan.out", "tmp/neo.pep", vcf, vepFile,
			          haplotype, tpm, proteinIndex, args.minTpmCount, "epitopes.tsv");

	if(!args.dirty) rmdir("tmp");

	return 0;
}