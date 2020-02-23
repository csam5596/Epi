#include <Logging.h>
#include <index/CdsIndex.h>
#include <util/DNAUtil.h>
#include "peptides/Peptidome.h"
#include "peptides/CuttingRange.h"
#include "peptides/Peptide.h"
#include "peptides/PeptideIterator.h"
#include <boost/fiber/barrier.hpp>

#include <boost/asio.hpp>
#include <boost/thread/latch.hpp>
#include <peptides/LineIterator.h>

Peptidome::Peptidome()
	: pool(4), selfPeptideFiles(), neoPeptideFiles(), selfPeptideStreams() {

	for(auto i = 0; i<4; i++){
		selfPeptideFiles.at(i) = "tmp/pep" + std::to_string(i+8) + ".txt";
		selfPeptideStreams.at(i) = std::make_unique<std::ofstream>(selfPeptideFiles.at(i));
		neoPeptideFiles.at(i) = "tmp/neopep" + std::to_string(i+8) + ".txt";
		std::ofstream o(neoPeptideFiles.at(i)); // create the files
	}
}

Peptidome::~Peptidome() {
	pool.close();
	pool.join();
}

void Peptidome::createPeptides(const Arguments &args, const Haplotype &hap, const VepFile& vep, const ProteinIndex &pi) {

	LOG_DEBUG("Creating peptides");
	CdsIndex cdsIndex(args.gtf, args.ref, vep.getTranscriptIds());
	const auto& changesByProtein = hap.getProteinChanges();

	for(auto&[proteinId, changeList] : changesByProtein){

		std::optional<Protein const*> refProtein = pi.getProtein(proteinId);
		if(!refProtein.has_value()){
			LOG_WARN("No reference protein found for " + proteinId);
			continue;
		}

		std::optional<std::string> codingSequence = cdsIndex.getCodingSequence(refProtein.value()->transcriptId);
		if(!codingSequence.has_value()){
			LOG_WARN("No coding sequence found for " + proteinId);
			continue;
		}

		const std::string& referenceProtein = refProtein.value()->aminoAcids;
		std::string translated = DNAUtil::transcriptToProtein(codingSequence.value());

		SANITY_CHECK(proteinsMatch(referenceProtein, translated, proteinId));

		std::for_each(changeList.cbegin(), changeList.cend(), [this, &codingSequence, &proteinId](const auto& changeSet){
			applyChanges(proteinId, codingSequence.value(), changeSet);
		});
	}

	LOG_DEBUG("Printing reference peptides");
	const std::unordered_map<ENSP, Protein>& refProteins = pi.getAllProteins();
	std::for_each(refProteins.cbegin(), refProteins.cend(), [this, &changesByProtein](const auto& kv){
		if(changesByProtein.find(kv.first) == changesByProtein.end()) {
			printSelfPeptides(kv.second.aminoAcids);
		}
	});

	std::for_each(selfPeptideStreams.begin(), selfPeptideStreams.end(), [](auto& file){ file->close(); });
	sortUniqueSelfPeptides();

}

bool Peptidome::proteinsMatch(const std::string &expected, const std::string &actual, const ENSP &ensp) {
	if(expected != actual){
		LOG_WARN("Proteins differ for " + ensp);
		LOG_WARN("\t[Expected]: " + expected);
		LOG_WARN("\t[Actual]:   " + expected);
		return false;
	}
	return true;
}

void Peptidome::sortUniqueSelfPeptides(){

	LOG_DEBUG("Sorting self peptides");

	boost::latch latch(4);

	static auto sortUnique = [](unsigned length, boost::latch& l){
		std::string filename = "tmp/pep" + std::to_string(length) + ".txt";
		const std::string sortCommand = "sort -u " + filename + " -o " + filename;
		boost::process::child sort(sortCommand);
		sort.wait();
		if(sort.exit_code() != 0){
			LOG_ERROR("[" + sortCommand + "] returned " + std::to_string(sort.exit_code()));
			exit(1);
		}

		LOG_DEBUG("Sorted " + filename);
		l.count_down();
	};

	for(auto i = 8; i < 12; i++){
		pool.submit([i, &latch]{ return sortUnique(i, latch); });
	}
	latch.wait();
}

void Peptidome::printSelfPeptides(const std::string &protein) {
	for(auto length = 8ul; length < 12; length++) {

		if (protein.length() < length) return;

		std::unique_ptr<std::ofstream>& out = selfPeptideStreams.at(length - 8);
		for (unsigned long i = 0; i < protein.length() - length; i++) {
			*out << std::string_view(protein).substr(i, length) << '\n';
		}
	}
}

void Peptidome::applyChanges(const std::string& id, const std::string &cds, const std::set<CodonChange> &changes) {

	unsigned somaticCount = std::count_if(changes.cbegin(), changes.cend(), [](const CodonChange& c){
		return c.vType == Variant::SOMATIC;
	});

	// only germline mutations
	if(somaticCount == 0){

		std::string selfCds = cds;
		bool check = std::all_of(changes.rbegin(), changes.rend(), [&selfCds](const CodonChange& cc){
			return cc(selfCds);
		});
		if(!check){
			LOG_WARN("Dropped Protein " + id);
			return;
		}
		printSelfPeptides(DNAUtil::transcriptToProtein(selfCds));

	}else{

		std::string selfCds = cds;
		std::string nonSelfCds = cds;
		std::deque<InfluenceRange> ir;

		bool check = std::all_of(changes.rbegin(), changes.rend(), [&selfCds, &nonSelfCds, &ir](const auto& cc){
			if(!cc(nonSelfCds)){
				return false;
			}
			if(cc.getShift() != 0) std::for_each(ir.begin(), ir.end(), [&cc](auto& i){ i.shift(cc.getShift()); });
			if(cc.vType == Variant::GERMLINE){
				cc(selfCds);
			}else{
				ir.emplace_front(cc);
			}
			return true;
		});

		if(!check){
			LOG_WARN("Dropped Protein " + id);
			return;
		}

		printSelfPeptides(DNAUtil::transcriptToProtein(selfCds));
		cutAndPrintNeoPeptides(DNAUtil::transcriptToProtein(nonSelfCds), ir);
	}


}

void Peptidome::cutAndPrintNeoPeptides(const std::string &protein, std::deque<InfluenceRange> &ir) {

	std::for_each(ir.begin(), ir.end(), [](auto& range){ range.convertPositions(); });

	// remove ranges after stop codon
	while(!ir.empty() && static_cast<unsigned long>(ir.back().firstStart) >= protein.length()){
		ir.pop_back();
	}

	for(unsigned pepLength = 8; pepLength < 12; pepLength++){
		if(protein.length() < pepLength) break;

		CuttingRange cr(pepLength, ir, protein);
		std::map<std::string, Peptide> peptides = cr(protein);

		const std::string& inFile = neoPeptideFiles.at(pepLength-8);
		std::ofstream out(inFile+".tmp");

		PeptideIterator pepIt(inFile);
		auto it = peptides.begin();

		while(pepIt.isGood() && it != peptides.end()){

			if(pepIt.peptide.getSequence() < it->first){
				out << pepIt.peptide << '\n';
				pepIt.advance();
			}else if(pepIt.peptide.getSequence() == it->first){
				pepIt.peptide.combine(it->second);
				it++;
			}else{
				out << it->second << '\n';
				it++;
			}
		}

		while(pepIt.isGood()){
			out << pepIt.peptide << '\n';
			pepIt.advance();
		}

		while(it != peptides.end()){
			out << it->second << '\n';
			it++;
		}

		pepIt.close();
		out.close();

		remove(inFile.c_str());
		rename((inFile + ".tmp").c_str(), inFile.c_str());
	}
}

void Peptidome::filterPeptides(){

	LOG_DEBUG("Filtering found peptides");

	boost::latch latch(4);

	static auto filter = [](const std::string& inFile, const std::string& maskFile, boost::latch& l){

		PeptideIterator it(inFile);
		LineIterator mask(maskFile);
		std::ofstream out(inFile + ".tmp");

		while(it.isGood() && mask.isGood()){
			if(mask.lineBuffer < it.peptide.getSequence()){
				mask.advance();
			}else if(mask.lineBuffer == it.peptide.getSequence()){
				it.advance();
			}else{
				out << it.peptide << '\n';
				it.advance();
			}
		}

		while(it.isGood()){
			out << it.peptide << '\n';
			it.advance();
		}

		it.close();
		out.close();

		remove(inFile.c_str());
		rename((inFile+".tmp").c_str(), inFile.c_str());
		LOG_DEBUG("Filtering " + inFile + " done");
		l.count_down();
	};

	for(auto i = 8; i < 12; i++){
		pool.submit([this, i, &latch]{ filter(neoPeptideFiles.at(i-8), selfPeptideFiles.at(i-8), latch); });
	}

	latch.wait();
	LOG_DEBUG("Filtering complete");
}

void Peptidome::mergeNeoPeptides(const std::string &outFile) const {

	std::ofstream out(outFile);
	if(out.fail()) throw std::runtime_error("Could not open " + outFile);

	unsigned long count = 0;
	out << "#Peptide\teffectIds\tImmunogenicity" << '\n';
	std::for_each(neoPeptideFiles.cbegin(), neoPeptideFiles.cend(), [&out, &count](const std::string& fileName){
		std::ifstream in(fileName);
		if(in.fail()) throw std::runtime_error("Could not open file " + fileName);
		std::string line;
		while(getline(in, line)){
			out << line << '\n';
			count++;
		}
	});

	LOG_DEBUG("Found " + std::to_string(count) + " neopeptides");
}