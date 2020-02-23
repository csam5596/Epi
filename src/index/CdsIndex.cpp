#include "index/CdsIndex.h"
#include "index/Exon.h"
#include <boost/process.hpp>
#include <iostream>

#include <sstream>
#include <algorithm>
#include <util/GtfUtil.h>
#include <Logging.h>
#include <util/DNAUtil.h>

CdsIndex::CdsIndex(const std::string &gtf, const std::string &ref, const std::set<ENST> &transcripts) {

	LOG_DEBUG("Building CDS-Index for " + std::to_string(transcripts.size()) + " transcripts");
	parseGtf(gtf, transcripts);
	parsePositions(ref);
	buildCodingSequences(transcripts);
	LOG_DEBUG("CDS Index built for " + std::to_string(codingSequences.size()) + " transcripts");

}

void CdsIndex::parseGtf(const std::string &file, const std::set<ENST> &transcripts) {

	namespace bp = boost::process;
	LOG_DEBUG("Loading exons and start codons from " + file);

	std::ofstream gtfArgs("tmp/gtfArgs.tmp");
	std::for_each(transcripts.cbegin(), transcripts.cend(), [&gtfArgs](const auto& s){ gtfArgs << s << std::endl; });

	bp::ipstream pipeStreamPositions;
	const std::string cPositionsCommand = "grep -F -f tmp/gtfArgs.tmp " + file;
	bp::child cPositions(cPositionsCommand, bp::std_out > pipeStreamPositions);

	std::string line;

	std::string chr, source, feature, quality;
	unsigned startPosition, endPosition;
	char strand;

	while(pipeStreamPositions && getline(pipeStreamPositions, line)){

		if(line.empty() || line.at(0) == '#') continue;
		std::istringstream iss(line);
		iss >> chr >> source >> feature  >> startPosition >> endPosition;

		if(feature != "exon" && feature != "start_codon") continue;

		const std::string transcriptId = GtfUtil::getTranscriptId(line);
		unsigned exonNumber = GtfUtil::getExonNumber(line);

		if(feature == "exon"){
			Exon e(exonNumber, chr, startPosition, endPosition);
			if(exons.find(transcriptId) == exons.end()){
				exons.emplace(transcriptId, std::set<Exon>{});
			}
			exons.at(transcriptId).emplace(std::move(e));

		}else if(feature == "start_codon"){
			iss >> quality >> strand;
			StartCodon c(exonNumber, chr, startPosition, endPosition, strand);
			if(startCodons.find(transcriptId) != startCodons.end()){
				LOG_WARN("Additional Start Codon for " + transcriptId + " [" + c.toString() + "] - using [" + startCodons.at(transcriptId).toString() + "]");
			}
			startCodons.emplace(transcriptId, std::move(c));
		}
	}
	cPositions.wait();
	if(cPositions.exit_code() != 0){
		LOG_ERROR("[" + cPositionsCommand + "] returned " + std::to_string(cPositions.exit_code()));
		exit(1);
	}
	LOG_DEBUG("Finished parsing exons and start codons");

}


void CdsIndex::parsePositions(const std::string &refFile) {

	namespace bp = boost::process;

	// write all positions to file
	std::ofstream posFile("tmp/positions.tmp");
	if(posFile.fail()) throw std::runtime_error("Could not open file [positions.tmp]");

	std::set<std::string> pos;

	std::for_each(exons.cbegin(), exons.cend(), [&pos](const auto& kv){
		std::transform(kv.second.cbegin(), kv.second.cend(), std::inserter(pos, pos.begin()), [](const Exon& e){
			return e.location;
		});
	});
	std::for_each(pos.cbegin(), pos.cend(), [&posFile](const auto& sv){ posFile << sv << std::endl; });

	LOG_DEBUG("Found " + std::to_string(pos.size()) + " unique positions");

	// get sequences from reference genome!
	bp::ipstream pipeStreamSequences;
	const std::string cCommand = "xargs --arg-file=tmp/positions.tmp samtools faidx " + refFile;
	bp::child c(cCommand, bp::std_out > pipeStreamSequences);

	std::string line;
	std::string current;

	while (pipeStreamSequences && std::getline(pipeStreamSequences, line) && !line.empty()){
		// new position
		if(line.at(0) == '>'){
			const std::string currentPosition = line.substr(1);
			sequencesByLocation.emplace(currentPosition, "");
			current = std::string_view(currentPosition);
			continue;
		}
		sequencesByLocation.at(current).append(line);
	}
	c.wait();
	if(c.exit_code() != 0){
		LOG_ERROR("[" + cCommand + "] returned " + std::to_string(c.exit_code()));
		exit(1);
	}

	LOG_DEBUG("Retrieved " + std::to_string(sequencesByLocation.size()) + " positions from reference genome");
}


void CdsIndex::buildCodingSequences(const std::set<ENST>& transcripts) {

	std::for_each(transcripts.cbegin(), transcripts.cend(), [this](const ENST& enst){

		if(startCodons.find(enst) == startCodons.end() || exons.find(enst) == exons.end()) return;

		const StartCodon& start = startCodons.at(enst);
		unsigned leading = 0;

		std::string cds;

		for(const Exon& e : exons.at(enst)){
			if(e.exonNumber < start.exonNumber) continue;
			if(e.exonNumber == start.exonNumber) {
				leading = (start.strand == 1) ? start.start - e.startPosition : e.endPosition - start.end;
			}

			std::string sequence = sequencesByLocation.at(e.location);
			if(start.strand == -1){
				std::reverse(sequence.begin(), sequence.end());
				DNAUtil::complementarySequence(sequence);
			}
			cds += sequence;
		}
		cds.erase(0, leading);
		codingSequences.emplace(enst, std::move(cds));
	});

}

std::optional<std::string> CdsIndex::getCodingSequence(const ENST &enst) const {
	if(codingSequences.find(enst) != codingSequences.end()){
		return codingSequences.at(enst);
	}else{
		return std::nullopt;
	}
}
