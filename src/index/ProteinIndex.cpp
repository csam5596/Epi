#include "index/ProteinIndex.h"
#include <fstream>
#include <Logging.h>

ProteinIndex::ProteinIndex(const std::string &file) {

	LOG_DEBUG("Loading proteins from " + file);
	std::ifstream in(file);
	if(in.fail()) throw std::runtime_error("Could not open file: " + file);

	std::string currentProtein;
	std::string line;

	while(getline(in, line)){

		// new protein
		if(line.at(0) == '>'){
			Protein p = Protein::parse(line);
			currentProtein = p.proteinId;
			proteinsById.emplace(p.proteinId, std::move(p));
		// amino acids
		}else{
			if(proteinsById.find(currentProtein) != proteinsById.end()){
				proteinsById.at(currentProtein).aminoAcids.append(line);
			}
		}
	}

	unsigned loadedProteins = proteinsById.size();
	LOG_DEBUG("Loaded " + std::to_string(loadedProteins) + " proteins");

	// remove invalid proteins
	for(auto it = proteinsById.begin(); it != proteinsById.end(); ){
		if(it->second.aminoAcids.find_first_of("XU*") != std::string::npos){
			proteinsById.erase(it++);
		}else{
			it++;
		}
	}

	LOG_DEBUG("Removed " + std::to_string(loadedProteins-proteinsById.size()) + " invalid proteins");

}

std::optional<Protein const*> ProteinIndex::getProtein(const ENSP &ensp) const {
	if(proteinsById.find(ensp) != proteinsById.end()){
		return &proteinsById.at(ensp);
	}else{
		return std::nullopt;
	}
}

const std::unordered_map<ENSP, Protein>& ProteinIndex::getAllProteins() const {
	return proteinsById;
}