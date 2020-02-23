#ifndef EPI_RANDOMFOREST_H
#define EPI_RANDOMFOREST_H

#include <vector>
#include <string>
#include "DecisionTree.h"
#include <functional>

using Encoder = std::function<std::vector<Predictor>(const std::string&)>;

class RandomForest {

public:
	explicit RandomForest(std::istringstream& iss, unsigned size = 0);
	void predictAll(const std::string& pepFile, const Encoder& encoder) const;

	[[nodiscard]] double predict(const std::vector<Predictor> &predictors) const;

private:
	std::vector<DecisionTree> forest;
	std::vector<Parameter> parameters;

};


#endif //EPI_RANDOMFOREST_H
