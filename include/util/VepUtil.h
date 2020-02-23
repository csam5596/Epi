#ifndef EPITOPE_VEPUTIL_H
#define EPITOPE_VEPUTIL_H

#include <iostream>
#include <string>
#include <boost/process.hpp>
#include <boost/format.hpp>
#include "Logging.h"

namespace VepUtil {

	static inline void predictEffects(const std::string& cache, const std::string& vcf, const std::string& out){

		namespace bp = boost::process;

		static auto command = boost::format("vep --offline --dir %1% -i %2% -o %3% --force_overwrite --fork 4 --protein");

		LOG_DEBUG("Predicting variant effects for " + vcf + " > " + out + ".");

		const std::string vepCommand = boost::str(command % cache % vcf % out);
		bp::child vep(vepCommand, bp::std_out > "log/vep.out", bp::std_err > "log/vep.err");
		vep.wait();
		if(vep.exit_code() != 0){
			LOG_ERROR("[" + vepCommand + "] returned " + std::to_string(vep.exit_code()));
			exit(1);
		}

		LOG_DEBUG("Calling VEP done.");
	}

	static inline bool filter(const std::vector<std::string>& tokens) {
		return  tokens[VepFile::CDS_POSITION] == "-" || tokens[VepFile::PROTEIN_POSITION] == "-" ||
				tokens[VepFile::AMINO_ACIDS] == "-" || tokens[VepFile::CODONS] == "-" ||
				tokens[VepFile::EXTRA].find("ENSP=") == std::string::npos ||
				tokens[VepFile::EXTRA].find("STRAND=") == std::string::npos ||
		        tokens[VepFile::EXTRA].find("cds_end_NF") != std::string::npos ||
		        tokens[VepFile::EXTRA].find("cds_start_NF") != std::string::npos;
	}
}

#endif //EPITOPE_VEPUTIL_H
