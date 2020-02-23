#ifndef EPITOPE_PHASER_H
#define EPITOPE_PHASER_H

#include <string>
#include <boost/process.hpp>
#include <boost/format.hpp>
#include "Logging.h"

namespace Phaser {

	static inline void phaseVariants(const std::string& phaserPath, const std::string& vcf, const std::string& dna,
									 const std::string& rna, unsigned dnaMapQ, unsigned rnaMapQ, const std::string& out){

		namespace bp = boost::process;

		static auto phaseCommand = boost::format(
				"python %1%/phaser.py --vcf %2% --bam %3% --pass_only 0 --sample Combined --baseq 10 --mapq %4% "
				"--paired_end 1 --process_slow 1 --threads 8 --include_indels 1 "
				"--gw_phase_vcf 1 --o %5%"
		);

		const std::string vcfGz = vcf + ".gz";
		const std::string bam = dna + (rna.empty() ? "" : ","+rna);
		const std::string mapq = std::to_string(dnaMapQ) + (rna.empty() ? "" : "," + std::to_string(rnaMapQ));


		// BGZIP
		const std::string bgzipCommand = "bgzip -f " + vcf;
		bp::child bgzip(bgzipCommand, bp::std_out > "log/bgzip.out", bp::std_err > "log/bgzip.err");
		bgzip.wait();
		if(bgzip.exit_code() != 0){
			LOG_ERROR("[" + bgzipCommand + "] returned " + std::to_string(bgzip.exit_code()));
			exit(0);
		}

		// TABIX
		const std::string tabixCommand = "tabix -p vcf " + vcfGz;
		bp::child tabix(tabixCommand, bp::std_out > "log/tabix.out", bp::std_err > "log/tabix.err");
		tabix.wait();
		if(tabix.exit_code() != 0){
			LOG_ERROR("[" + tabixCommand + "] returned " + std::to_string(tabix.exit_code()));
			exit(1);
		}

		// PHASER
		LOG_DEBUG("Phasing variants from " + vcf);

		const std::string command = boost::str(phaseCommand % phaserPath % vcfGz % bam % mapq % out);
		std::ofstream phaserStdout("log/phaser.out");
		if(phaserStdout.fail()) throw std::runtime_error("Could not open file log/phaser.out");

		std::string line;

		bp::ipstream outStream;
		bp::child phaser(command, bp::std_out > outStream, bp::std_err > "log/phaser.err");

		while(outStream && getline(outStream, line)){

			phaserStdout << line << std::endl;
			if(line.find("processing chromosome") != std::string::npos){
				LOG_DEBUG(line);
			}
		}
		phaser.wait();
		if(phaser.exit_code() != 0){
			LOG_ERROR("[" + command + "] returned " + std::to_string(tabix.exit_code()));
			exit(1);
		}

		LOG_DEBUG("Phasing completed");
	}

}

#endif //EPITOPE_PHASER_H
