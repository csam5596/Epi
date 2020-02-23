#ifndef EPITOPE_GTFUTIL_H
#define EPITOPE_GTFUTIL_H

#include <string>
#include <regex>

namespace GtfUtil {

	static inline unsigned getExonNumber(const std::string& gtfLine){

		size_t exonStart = gtfLine.find("exon_number \"");
		size_t exonEnd = gtfLine.find('\"', exonStart+13);

		return std::stoul(gtfLine.substr(exonStart+13, exonEnd-exonStart-13));
	}

	static inline std::string getTranscriptId(const std::string& gtfLine){

		size_t transcriptStart = gtfLine.find("transcript_id \"");
		size_t transcriptEnd = gtfLine.find('\"', transcriptStart+15);

		return gtfLine.substr(transcriptStart+15, transcriptEnd-transcriptStart-15);
	}
}

#endif //EPITOPE_GTFUTIL_H
