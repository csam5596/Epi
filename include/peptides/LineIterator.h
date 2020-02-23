#ifndef EPITOPE_LINEITERATOR_H
#define EPITOPE_LINEITERATOR_H

#include <string>
#include <memory>

class LineIterator {

public:
	explicit LineIterator(const std::string& file);

	void advance();
	void close();
	[[nodiscard]] bool isGood() const;

public:
	std::string lineBuffer;

private:
	std::unique_ptr<std::ifstream> in;

};


#endif //EPITOPE_LINEITERATOR_H
