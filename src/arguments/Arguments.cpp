#include "arguments/Arguments.h"
#include <iostream>
#include <Logging.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/join.hpp>

Arguments::Arguments(int argc, char *argv[]) {

	namespace po = boost::program_options;

	po::options_description desc("Program options");
	desc.add_options()
			("help,h", "produce help message")
			("dna", po::value<std::string>(&dna),
			 "Tumor DNA - aligned, sorted and indexed (.bam)")
			("rna", po::value<std::string>(&rna)->required(),
			 "Tumor RNA - aligned, sorted and indexed (.bam)")
			("somatic", po::value<std::vector<std::string>>(&somaticVariants)->multitoken()->required(),
			 "vcf files containing somatic variants (.vcf | .vcf.gz) - space separated list")
			("vc", po::value<std::vector<std::string>>(&variantCallers)->multitoken()->required(),
			 "name of variant callers by vcf file [Mutect2, Strelka, Varscan]")
			("germline", po::value<std::string>(&germlineVariants)->required(),
			 "vcf file containing variants called by GATK Haplotypecaller")
			("proteins", po::value<std::string>(&proteins)->required(),
			 ("File containing all protein ENSP numbers and corresponding amino acid sequences"))
			("gtf", po::value<std::string>(&gtf)->required(),
			 ("GTF file containing positions"))
			("ref", po::value<std::string>(&ref)->required(),
			 ("reference genome file"))
			("hla", po::value<std::string>(&hla)->required(),
			 ("hla output as derived from Optitype"))
			("tpm", po::value<std::string>(&tpm)->required(),
			 ("Kallisto output file containing TPM (.tsv)"))
			("phaser", po::value<std::string>(&phaser)->required(),
			 "Path to phaser.py")
			("vepCacheDir", po::value<std::string>(&vepCacheDir)->required(),
			 "Path to cache for Variant effect Predictor")

			// Optional Parameters
			("minVcCount", po::value<int>(&minVcCount)->default_value(2),
			 "Minimum Number of variant callers that need to find a specific variant "
			 "for it to be considered for further analysis")
			("maxNormalVaf", po::value<float>(&maxNormalVaf)->default_value(0.05),
			 "Maximum Variant Allele Frequency for Normal sample")
			("minTumorVaf", po::value<float>(&minTumorVaf)->default_value(0.05),
			 "Minimum Variant Allele Frequency for Tumor Sample")
			("minRcAbs", po::value<unsigned>(&minReadCountAbs)->default_value(1),
			 "Minimum absolute number of reads a phased Variant needs to appear on")
			("minRcRel", po::value<float>(&minReadCountRel)->default_value(0.2f),
			 "Minimum Fraction of reads a phased Variants needs to be found on")
			("minTpmCount", po::value<unsigned>(&minTpmCount)->default_value(1),
			 "Minimum Tpm Count for a Transcript to be considered expressed")
			("dnaMapQ", po::value<unsigned>(&dnaMapQ)->default_value(1),
			 "Minimum DNA mapping quality for reads to be considered for variant phasing")
			("rnaMapQ", po::value<unsigned>(&rnaMapQ)->default_value(255),
			 "Minimum RNA mapping quality for reads to be considered for variant phasing")
			("dirty", po::bool_switch(&dirty)->default_value(false),
			"Keep temporary files");

	po::variables_map vm;

	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);

		if (vm.count("help")) {
			std::cout << desc << std::endl;
			exit(0);
		}

		po::notify(vm);

		printValues();

	} catch (po::error &e) {
		LOG_ERROR(e.what());
		std::cout << desc << std::endl;
		exit(1);
	}
}

void Arguments::printValues() {

	std::cout << "Arguments: " << std::endl;
	std::cout << std::setw(15) << std::left << "--dna" << '\t' << dna << std::endl;
	std::cout << std::setw(15) << std::left << "--rna" << '\t' << rna << std::endl;

	std::cout << std::setw(15) << std::left << "--somatic" << '\t' << somaticVariants[0] << std::endl;
	for(unsigned i = 1; i<somaticVariants.size(); i++){
		std::cout << std::setw(15) << ' ' << '\t' << somaticVariants[i] << std::endl;
	}

	std::cout << std::setw(15) << std::left << "--vc" << '\t' << variantCallers[0] << std::endl;
	for(unsigned i = 1; i<variantCallers.size(); i++){
		std::cout << std::setw(15) << ' ' << '\t' << variantCallers[i] << std::endl;
	}

	std::cout << std::setw(15) << std::left << "--germline" << '\t' << germlineVariants << std::endl;
	std::cout << std::setw(15) << std::left << "--proteins" << '\t' << proteins << std::endl;
	std::cout << std::setw(15) << std::left << "--gtf" << '\t' << gtf << std::endl;
	std::cout << std::setw(15) << std::left << "--ref" << '\t' << ref << std::endl;
	std::cout << std::setw(15) << std::left << "--hla" << '\t' << hla << std::endl;
	std::cout << std::setw(15) << std::left << "--tpm" << '\t' << tpm << std::endl;
	std::cout << std::setw(15) << std::left << "--phaser" << '\t' << phaser << std::endl;
	std::cout << std::setw(15) << std::left << "--vepCacheDir" << '\t' << vepCacheDir << std::endl;

	std::cout << std::endl;
	std::cout << "Cutoff Parameters: " << std::endl;

	std::cout << std::setw(15) << std::left << "--minVcCount" << '\t' << minVcCount << std::endl;
	std::cout << std::setw(15) << std::left << "--maxNormalVaf" << '\t' << maxNormalVaf << std::endl;
	std::cout << std::setw(15) << std::left << "--minTumorVaf" << '\t' << minTumorVaf << std::endl;
	std::cout << std::setw(15) << std::left << "--minRcAbs" << '\t' << minReadCountAbs << std::endl;
	std::cout << std::setw(15) << std::left << "--minRcRel" << '\t' << minReadCountRel << std::endl;
	std::cout << std::setw(15) << std::left << "--minTpmCount" << '\t' << minTpmCount << std::endl;
	std::cout << std::setw(15) << std::left << "--dnaMapQ" << '\t' << dnaMapQ << std::endl;
	std::cout << std::setw(15) << std::left << "--rnaMapQ" << '\t' << rnaMapQ << std::endl;
	std::cout << std::endl;

	if(dirty){
		std::cout << std::setw(15) << std::left << "--dirty" << std::endl << std::endl;
	}

}