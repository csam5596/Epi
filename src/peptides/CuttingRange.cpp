#include "peptides/CuttingRange.h"
#include <algorithm>
#include <map>

CuttingRange::CuttingRange(unsigned length, std::deque<InfluenceRange> &ir, const std::string &protein)
	: length(length) {

	std::map<long, unsigned> effectsByStart;
	std::map<long, unsigned> effectsByEnd;
	std::set<unsigned> effectIdBuffer;

	std::for_each(ir.cbegin(), ir.cend(), [&effectsByStart, &effectsByEnd, &protein, length](const InfluenceRange& i){
		effectsByStart.emplace(std::max(i.firstStart - (length-1), 0l), i.variantEffectId);
		effectsByEnd.emplace(std::min(static_cast<long>(protein.length()) - length, i.lastStart), i.variantEffectId);
	});

	auto itStart = effectsByStart.begin();
	auto itEnd = effectsByEnd.begin();

	while(itStart != effectsByStart.end() || itEnd != effectsByEnd.end()){

		unsigned s;
		unsigned e;

		if(itStart == effectsByStart.end() || itEnd->first < itStart->first){
			s = itEnd->first;
			effectIdBuffer.erase(itEnd->second);
			std::advance(itEnd, 1);
		}else{
			s = itStart->first;
			effectIdBuffer.emplace(itStart->second);
			std::advance(itStart, 1);
		}

		if(itStart == effectsByStart.end() && itEnd == effectsByEnd.end()) break;
		if(effectIdBuffer.empty()) continue;

		e = (itStart == effectsByStart.end() || itEnd->first < itStart->first) ? itEnd->first : itStart->first-1;

		ranges.emplace(std::make_pair(s,e), effectIdBuffer);
	}
}

std::map<std::string, Peptide> CuttingRange::operator()(const std::string &protein) {

	std::map<std::string, Peptide> peptides;

	for(auto& cr : ranges){

		long start = cr.first.first;
		long end = cr.first.second;

		for(long i = start; i<=end; i++){
			std::string  p = protein.substr(i, length);
			if(peptides.find(p) == peptides.end()){
				Peptide peptide(p, cr.second);
				peptides.emplace(peptide.getSequence(), std::move(peptide));
			}else{
				peptides.at(p).addVariantEffects(cr.second);
			}
		}
	}
	return peptides;
}