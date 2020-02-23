#ifndef EPITOPE_DNAUTIL_H
#define EPITOPE_DNAUTIL_H

#include <unordered_map>
#include <iostream>
#include <algorithm>

namespace DNAUtil {

	static inline char complementaryBase(const char base){
		static const std::unordered_map<char,char> bases = {
				{'A', 'T'},
				{'T', 'A'},
				{'G', 'C'},
				{'C', 'G'}
		};
		if(bases.find(base) == bases.end()){
			std::cerr << "Unknown Base: " << base << std::endl;
			return '?';
		}
		return bases.at(base);
	}

	static inline void complementarySequence(std::string& sequence){
		std::transform(sequence.begin(), sequence.end(), sequence.begin(), complementaryBase);
	}


	static inline char translateDNACodon(const std::string_view& codon){
		static std::unordered_map<std::string_view, char> dnaCodons = {
				{"TTT", 'F'}, {"TCT", 'S'}, {"TAT", 'Y'}, {"TGT", 'C'},
				{"TTC", 'F'}, {"TCC", 'S'}, {"TAC", 'Y'}, {"TGC", 'C'},
				{"TTA", 'L'}, {"TCA", 'S'}, {"TAA", '*'}, {"TGA", '*'},
				{"TTG", 'L'}, {"TCG", 'S'}, {"TAG", '*'}, {"TGG", 'W'},

				{"CTT", 'L'}, {"CCT", 'P'}, {"CAT", 'H'}, {"CGT", 'R'},
				{"CTC", 'L'}, {"CCC", 'P'}, {"CAC", 'H'}, {"CGC", 'R'},
				{"CTA", 'L'}, {"CCA", 'P'}, {"CAA", 'Q'}, {"CGA", 'R'},
				{"CTG", 'L'}, {"CCG", 'P'}, {"CAG", 'Q'}, {"CGG", 'R'},

				{"ATT", 'I'}, {"ACT", 'T'}, {"AAT", 'N'}, {"AGT", 'S'},
				{"ATC", 'I'}, {"ACC", 'T'}, {"AAC", 'N'}, {"AGC", 'S'},
				{"ATA", 'I'}, {"ACA", 'T'}, {"AAA", 'K'}, {"AGA", 'R'},
				{"ATG", 'M'}, {"ACG", 'T'}, {"AAG", 'K'}, {"AGG", 'R'},

				{"GTT", 'V'}, {"GCT", 'A'}, {"GAT", 'D'}, {"GGT", 'G'},
				{"GTC", 'V'}, {"GCC", 'A'}, {"GAC", 'D'}, {"GGC", 'G'},
				{"GTA", 'V'}, {"GCA", 'A'}, {"GAA", 'E'}, {"GGA", 'G'},
				{"GTG", 'V'}, {"GCG", 'A'}, {"GAG", 'E'}, {"GGG", 'G'},
		};
		if(codon.size() != 3 || dnaCodons.find(codon) == dnaCodons.end()){
			std::cerr << "Invalid codon: " << codon << std::endl;
			return '?';
		}
		return dnaCodons.at(codon);
	}

	static inline std::string transcriptToProtein(const std::string& cds){
		std::vector<char> aminoAcids;
		aminoAcids.reserve(cds.length() / 3);
		for(unsigned long idx = 0; idx < cds.length()-2; idx+=3){
			std::string_view view(cds);
			char aa = translateDNACodon(view.substr(idx, 3));
			if(aa == '*') break;
			aminoAcids.emplace_back(aa);
		}
		aminoAcids.shrink_to_fit();
		return std::string(aminoAcids.begin(), aminoAcids.end());
	}

}

#endif //EPITOPE_DNAUTIL_H
