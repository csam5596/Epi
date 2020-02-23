#ifndef EPI_DECISIONTREE_H
#define EPI_DECISIONTREE_H

#include <vector>
#include <optional>
#include <fstream>
#include "Parameter.h"

class DecisionTree {

public:
	DecisionTree(const std::vector<Parameter>& parameters, std::istream& in);
	[[nodiscard]] bool predict(const std::vector<Predictor>& predictors) const;

private:
	const std::vector<Parameter>& parameters;

	std::vector<std::optional<int>> splitVars;
	std::vector<std::variant<unsigned, double, bool>> splitPoints;
	std::vector<int> left;
	std::vector<int> right;

};


#endif //EPI_DECISIONTREE_H
