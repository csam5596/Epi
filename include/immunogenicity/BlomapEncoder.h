#ifndef EPI_BLOMAPENCODER_H
#define EPI_BLOMAPENCODER_H

#include <array>
#include <variant>
#include <vector>
#include "Parameter.h"

namespace BlomapEncoder {

	static inline std::array<double, 5> encodeAminoAcid(char aminoAcid){
		switch(aminoAcid){
			case 'A': return {-0.57, 0.39, -0.96, -0.61, -0.69};
			case 'R': return {-0.4, -0.83, -0.61, 1.26, -0.28};
			case 'N': return {-0.7, -0.63, -1.47, 1.02, 1.06};
			case 'D': return {-1.62, -0.52, -0.67, 1.02, 1.47};
			case 'C': return {0.07, 2.04, 0.65, -1.13, -0.39};

			case 'Q': return {-0.05, -1.50, -0.67, 0.49, 0.21};
			case 'E': return {-0.64, -1.59, -0.39, 0.69, 1.04};
			case 'G': return {-0.90, 0.87, -0.36, 1.08, 1.95};
			case 'H': return {0.73, -0.67, -0.42, 1.13, 0.99};
			case 'I': return {0.59, 0.79, 1.44, -1.90, -0.93};

			case 'L': return {0.65, 0.84, 1.25, -0.99, -1.90};
			case 'K': return {-0.64, -1.19, -0.65, 0.68, -0.13};
			case 'M': return {0.76, 0.05, 0.06, -0.62, -1.59};
			case 'F': return {1.87, 1.04, 1.28, -0.61, -0.16};
			case 'P': return {-1.82, -0.63, 0.32, 0.03, 0.68};

			case 'S': return {-0.39, -0.27, -1.51, -0.25, 0.31};
			case 'T': return {-0.04, -0.3, -0.82, -1.02, -0.04};
			case 'W': return {1.38, 1.69, 1.91, 1.07, -0.05};
			case 'Y': return {1.75, 0.11, 0.65, 0.21, -0.41};
			case 'V': return {-0.02, 0.30, 0.97, -1.55, -1.16};
			default:
				return {0.0, 0.0, 0.0, 0.0, 0.0};
		}
	}

	static inline std::vector<Predictor> encode(const std::string& peptide){

		std::vector<Predictor> encoding;
		encoding.reserve(55);

		for(unsigned i = 0; i<11; i++){

			char aa = peptide.length() > i ? peptide.at(i) : 'X';
			std::array<double, 5> enc = encodeAminoAcid(aa);

			for(unsigned j = 0; j<5; j++){
				encoding.emplace_back(enc.at(j));
			}
		}

		return encoding;
	}

}


#endif //EPI_BLOMAPENCODER_H
