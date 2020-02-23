#ifndef EPITOPE_HLATYPES_H
#define EPITOPE_HLATYPES_H

#include <string>
#include <set>

struct HlaTypes {

public:
	static HlaTypes parse(const std::string& file);

	[[nodiscard]] std::string toString() const;

private:
	explicit HlaTypes(std::set<std::string>& types);

private:
	std::set<std::string> types;

};


#endif //EPITOPE_HLATYPES_H
