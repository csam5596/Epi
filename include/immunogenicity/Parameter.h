#ifndef EPI_PARAMETER_H
#define EPI_PARAMETER_H

#include <variant>

#include "ParameterType.h"

using Predictor = std::variant<unsigned, double>;

struct Parameter {
	explicit Parameter(unsigned nLevels);

	[[nodiscard]] bool compare(const std::variant<unsigned, double, bool>& splitPoint, const Predictor& predictor) const;

	const ParameterType type;
	const unsigned nLevels;

};


#endif //EPI_PARAMETER_H
