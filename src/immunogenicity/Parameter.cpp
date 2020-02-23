#include "immunogenicity/Parameter.h"

#include <stdexcept>


Parameter::Parameter(unsigned nLevels)
	: type(nLevels == 1 ? ParameterType::NUMERIC : ParameterType::FACTOR), nLevels(nLevels){}

bool Parameter::compare(const std::variant<unsigned, double, bool> &splitPoint, const Predictor &predictor) const {

	switch(type){

		case ParameterType::FACTOR:
			// see: https://www.rdocumentation.org/packages/randomForest/versions/4.6-14/topics/getTree
			return (std::get<0>(splitPoint) & (1u << std::get<0>(predictor)));

		case ParameterType::NUMERIC:
			return std::get<1>(predictor) <= std::get<1>(splitPoint);

		default:
			throw std::runtime_error("Invalid Parameter type");
	}

}