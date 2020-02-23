#include <sstream>
#include "immunogenicity/DecisionTree.h"

DecisionTree::DecisionTree(const std::vector<Parameter> &parameters, std::istream &in)
	: parameters(parameters){

	std::string line;

	unsigned node, l, r, var, pred;
	double point;
	int status;

	while(getline(in, line)){
		if(line.empty() || line.find("##") != std::string::npos){
			break;
		}

		std::istringstream iss(line);
		iss >> node >> l >> r >> var >> point >> status >> pred;

		left.emplace_back(l - 1);
		right.emplace_back(r - 1);

		// inner node
		if(status == 1){
			splitVars.emplace_back(var - 1);
			switch(parameters.at(var - 1).type){


				case ParameterType::FACTOR:
					splitPoints.emplace_back(static_cast<unsigned>(point));
					break;

				case ParameterType::NUMERIC:
					splitPoints.emplace_back(point);
					break;

				default:
					throw std::runtime_error("Invalid parameter type");

			}
		// leaf node
		}else if(status == -1){
			splitVars.emplace_back(std::nullopt);
			splitPoints.emplace_back(pred == 2);
		}else{
			throw std::runtime_error("Invalid node type: " + std::to_string(status));
		}
	}

}

bool DecisionTree::predict(const std::vector<Predictor> &predictors) const {

	int index = 0;

	while(index != -1){
		if(splitVars.at(index).has_value()){

			int splitParIndex = splitVars.at(index).value();
			parameters.at(splitParIndex).compare(splitPoints.at(index), predictors.at(splitParIndex)) ?
				index = left.at(index) :
				index = right.at(index);
		}else{
			return std::get<2>(splitPoints.at(index));
		}
	}
	throw std::runtime_error("Decision Tree Prediction failed");
}