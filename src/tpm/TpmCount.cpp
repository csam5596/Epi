#include "tpm/TpmCount.h"

#include <fstream>
#include <sstream>
#include <Logging.h>

TpmCount TpmCount::parse(const std::string &file) {

	std::ifstream in(file);
	if(in.fail()) throw std::runtime_error("Could not open " + file);

	LOG_DEBUG("Parsing TPM from " + file + ".");

	std::unordered_map<ENST, float> tpmCounts;

	std::string line;
	getline(in, line); // skip header

	std::string targetId;
	unsigned length;
	float effLength, estCounts, tpm;

	while(getline(in, line)){

		std::istringstream iss(line);
		iss >> targetId >> length >> effLength >> estCounts >> tpm;

		// clean targetId - remove .version
		targetId = targetId.substr(0, targetId.find('.'));
		tpmCounts.emplace(std::move(targetId), tpm);
	}

	LOG_DEBUG(std::to_string(tpmCounts.size()) + " transcripts loaded.");
	return TpmCount(tpmCounts);

}


TpmCount::TpmCount(std::unordered_map<ENST, float> &counts)
	: tpmCounts(std::move(counts)){}

std::optional<float> TpmCount::getTpm(const ENST& transcriptId) const {
	if(tpmCounts.find(transcriptId) != tpmCounts.end()){
		return tpmCounts.at(transcriptId);
	}
	return std::nullopt;
}