#ifndef EPITOPE_HAPLOBLOCK_H
#define EPITOPE_HAPLOBLOCK_H

#include <string>
#include <vector>

struct Haploblock {

public:
	static Haploblock parse(const std::string& line);

private:
	Haploblock(std::vector<unsigned> &variantIds, std::vector<std::string> &allelesA,
	           std::vector<std::string> &allelesB, unsigned countA, unsigned countB);

public:
	std::vector<unsigned> variantIds;
	std::vector<std::string> allelesA;
	std::vector<std::string> allelesB;
	unsigned countA;
	unsigned countB;
};


#endif //EPITOPE_HAPLOBLOCK_H
