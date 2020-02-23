#include <algorithm>
#include <sstream>
#include <Logging.h>
#include <peptides/Peptide.h>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include "immunogenicity/RandomForest.h"


RandomForest::RandomForest(std::istringstream& iss, unsigned int size) {

	boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
	in.push(boost::iostreams::gzip_decompressor());
	in.push(iss);

	std::istream iStream(&in);

	LOG_DEBUG("Loading Random Forest");

	std::string line;
	unsigned level;

	getline(iStream, line); // skip first line

	getline(iStream, line);
	std::istringstream issLine(line);

	while(issLine >> level){
		parameters.emplace_back(level);
	}

	if(size != 0) forest.reserve(size);

	while(getline(iStream, line)){
		forest.emplace_back(parameters, iStream);
	}

	LOG_DEBUG(std::to_string(forest.size()) + " trees loaded");

}

void RandomForest::predictAll(const std::string &pepFile, const Encoder &encoder) const {

	std::ifstream in(pepFile);
	if(in.fail()) throw std::runtime_error("Could not open " + pepFile);

	std::ofstream out(pepFile + ".tmp");
	std::string line;

	getline(in, line); // skip header
	out << "#Peptide\teffectIds\tImmunogenicity" << '\n';

	while(getline(in, line)){

		Peptide peptide = Peptide::parse(line);
		double score = predict(encoder(peptide.getSequence()));
		peptide.setScore(score);
		out << peptide << '\n';
	}

	remove(pepFile.c_str());
	rename((pepFile + ".tmp").c_str(), pepFile.c_str());

}

double RandomForest::predict(const std::vector<Predictor> &predictors) const {

	long posCount = std::count_if(forest.cbegin(), forest.cend(), [&predictors](const auto& tree){
		return tree.predict(predictors);
	});

	return static_cast<double>(posCount) / static_cast<double>(forest.size());
}