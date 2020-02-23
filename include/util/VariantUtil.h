#ifndef EPITOPE_VARIANTUTIL_H
#define EPITOPE_VARIANTUTIL_H

#include <string>
#include "vcf/Variant.h"

namespace VariantUtil{


	static inline float getFloatValue(const std::string& variantString, const std::string& key, float def){

		size_t pos = variantString.find(key);
		if(pos == std::string::npos) return def;

		size_t end = variantString.find_first_of(";\t", pos + key.length());
		std::string_view value = std::string_view(variantString).substr(pos + key.length(), end - pos - key.length());
		return std::stof(value.data());
	}

	static inline float getNormalVaf(const std::string& variantString){
		return getFloatValue(variantString, "NAF=", 1.f);
	}

	static inline float getTumorVaf(const std::string& variantString){
		return getFloatValue(variantString, "TAF=", 1.f);
	}

	static inline void addInfoValue(std::string& variantString, const std::string& name, const std::string& value){

		size_t pos = variantString.find('\t');
		short count = 1;

		while( pos != std::string::npos && count != 8){
			pos =variantString.find('\t', pos + 1);
			count++;
		}

		std::string kvString = ";" + name + "=" + value + '\t';
		variantString.replace(pos, 1, kvString);

	}

	static inline int getVariantCallerCount(const std::string& variantString){
		static const std::string key = "VC_COUNT=";

		size_t pos = variantString.find(key);
		if(pos == std::string::npos) return 1; // at least one vc has to have found it

		size_t end = variantString.find_first_of(";\t", pos + key.length());
		std::string_view value = std::string_view(variantString).substr(pos + key.length(), end - pos - key.length());
		return std::stoi(value.data());
	}

	static inline Variant::VariantType getVariantType(const std::string& variantString){

		static const std::string key = "VT=";

		size_t pos = variantString.find(key);

		size_t end = variantString.find_first_of(";\t", pos + key.length());
		std::string_view value = std::string_view(variantString).substr(pos + key.length(), end - pos - key.length());
		return (value == "SOMATIC") ? Variant::VariantType::SOMATIC : Variant::VariantType::GERMLINE;

	}

	static inline float getReadCountNormal(const std::string& variantString){
		return getFloatValue(variantString, "RCN=", -1.f);
	}

	static inline float getReadCountTumor(const std::string& variantString){
		return getFloatValue(variantString, "RCT=", -1.f);
	}

}

#endif //EPITOPE_VARIANTUTIL_H
