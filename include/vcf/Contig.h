#ifndef EPITOPE_CONTIG_H
#define EPITOPE_CONTIG_H


#include <string_view>
#include <variant>

struct Contig {

public:
	Contig() = default;
	explicit Contig(const std::string_view& chromosome);

	[[nodiscard]] std::string getName() const;
	bool operator<(const Contig &rhs) const;
	bool operator==(const Contig &rhs) const;

private:
	std::variant<unsigned long, char> chrNumber;

};


#endif //EPITOPE_CONTIG_H
