#ifndef EPITOPE_TPMCOUNT_H
#define EPITOPE_TPMCOUNT_H

#include <string>
#include <unordered_map>
#include <optional>
#include <Types.h>

class TpmCount {

public:
	static TpmCount parse(const std::string& file);

	[[nodiscard]] std::optional<float> getTpm(const ENST& transcriptId) const;

private:
	explicit TpmCount(std::unordered_map<ENST, float>& counts);

private:
	std::unordered_map<ENST, float> tpmCounts;

};


#endif //EPITOPE_TPMCOUNT_H
