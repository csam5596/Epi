#ifndef EPITOPE_NETMHCPANITERATOR_H
#define EPITOPE_NETMHCPANITERATOR_H

#include <string>
#include <memory>
#include <fstream>

class NetMHCpanIterator {


public:
	explicit NetMHCpanIterator(const std::string& file);

	void advance();
	[[nodiscard]] bool isGood() const;

private:
	void parseLine();

public:

	// line values
	std::string hla;
	std::string peptide;
	std::string core;
	std::string icore;
	float score = 0.f;
	float affinity = 0.f;
	float rank = 0.f;
	float exp = 0.f;
	std::string bindingLevel;

private:
	std::unique_ptr<std::ifstream> in;
	std::string lineBuffer;
};


#endif //EPITOPE_NETMHCPANITERATOR_H
