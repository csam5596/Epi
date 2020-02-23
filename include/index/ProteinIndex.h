#ifndef EPITOPE_PROTEININDEX_H
#define EPITOPE_PROTEININDEX_H

#include <string>
#include <optional>
#include <unordered_map>
#include <Types.h>
#include "Protein.h"

class ProteinIndex {

public:
	explicit ProteinIndex(const std::string& file);

	[[nodiscard]] std::optional<Protein const*> getProtein(const ENSP& ensp) const;
	[[nodiscard]] const std::unordered_map<ENSP, Protein>& getAllProteins() const;

private:
	std::unordered_map<ENSP, Protein> proteinsById;

};


#endif //EPITOPE_PROTEININDEX_H
