#ifndef EPITOPE_PEPTIDOME_H
#define EPITOPE_PEPTIDOME_H


#include <boost/thread/executors/basic_thread_pool.hpp>
#include <Types.h>
#include <index/ProteinIndex.h>
#include <arguments/Arguments.h>
#include <phasing/Haplotype.h>
#include "InfluenceRange.h"
#include <vep/VepFile.h>

#define SANITY_CHECK(x) if(!(x)) continue

class Peptidome {

public:
	Peptidome();
	~Peptidome();

	void createPeptides(const Arguments &args, const Haplotype &hap, const VepFile& vep, const ProteinIndex &pi);
	void filterPeptides();
	void mergeNeoPeptides(const std::string& outFile) const;


private:
	static bool proteinsMatch(const std::string& expected, const std::string& actual, const ENSP& ensp);

	void printSelfPeptides(const std::string& protein);
	void cutAndPrintNeoPeptides(const std::string& protein, std::deque<InfluenceRange>& ir);
	void sortUniqueSelfPeptides();
	void applyChanges(const std::string& id, const std::string& cds, const std::set<CodonChange>& changes);

private:
	boost::executors::basic_thread_pool pool;
	std::array<std::string, 4> selfPeptideFiles;
	std::array<std::string, 4> neoPeptideFiles;
	std::array<std::unique_ptr<std::ofstream>, 4> selfPeptideStreams;

};


#endif //EPITOPE_PEPTIDOME_H
