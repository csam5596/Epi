#ifndef EPITOPE_VCFUTIL_H
#define EPITOPE_VCFUTIL_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <map>

#include "Logging.h"
#include "vcf/VcfLineIterator.h"
#include "VariantUtil.h"

namespace VcfUtil {

	static inline void filterAndCombine(const std::vector<std::string>& files, const std::vector<std::string>& vcs,
	                                    unsigned minCount, float maxNormalVaf, float minTumorVaf, const std::string& outFile){

		std::ofstream out(outFile);
		if(out.fail()) throw std::runtime_error("Could not open " + outFile);

		LOG_DEBUG("Filtering and combining files ...");
		for(unsigned long i = 0; i<files.size(); i++){
			LOG_DEBUG("\tParsing " + files[i] + " - VariantCaller: " + vcs[i]);
		}

		std::vector<VcfLineIterator> iterators;
		std::transform(files.cbegin(), files.cend(), vcs.begin(), std::back_inserter(iterators),
               [maxNormalVaf, minTumorVaf](const std::string& file, const std::string& vc) {
							return VcfLineIterator(file, vc, maxNormalVaf, minTumorVaf);
				}
		);

		// combine headers
		std::vector<VcfHeader> headers;
		std::for_each(iterators.cbegin(), iterators.cend(), [&headers](const VcfLineIterator& i){
			headers.emplace_back(i.header);
		});
		VcfHeader header = VcfHeader::merge(headers);
		header.addInfo("NAF", "1", "Float", "Average normal allele frequency");
		header.addInfo("TAF", "1", "Float", "Average tumor allele frequency");
		header.addInfo("VC_COUNT", "1", "Integer", "Number of variant confirming the mutation");
		header.addInfo("RCN", "1", "Float", "Average read count for Normal Sample");
		header.addInfo("RCT", "1", "Float", "Average read count for Tumor Sample");
		out << header << std::endl;

		// combine variants
		std::vector<VcfLineIterator*> matching;

		while(!iterators.empty()){
			auto it = std::min_element(iterators.begin(), iterators.end());
			for(auto& lineIterator : iterators){
				if(lineIterator == *it) matching.emplace_back(&lineIterator);
			}

			if(matching.size() >= minCount){

				float normalVaf = 0.f;
				float tumorVaf = 0.f;
				float rcNormal = 0.f;
				float rcTumor = 0.f;

				std::for_each(matching.cbegin(), matching.cend(), [&](const VcfLineIterator* it){
					normalVaf += it->vafNormal;
					tumorVaf += it->vafTumor;
					rcNormal += it->rcNormal;
					rcTumor += it->rcTumor;
				});

				normalVaf /= matching.size();
				tumorVaf /= matching.size();
				rcNormal /= matching.size();
				rcTumor /= matching.size();

				VariantUtil::addInfoValue(it->lineBuffer, "NAF", std::to_string(normalVaf));
				VariantUtil::addInfoValue(it->lineBuffer, "TAF", std::to_string(tumorVaf));
				VariantUtil::addInfoValue(it->lineBuffer, "VC_COUNT", std::to_string(matching.size()));
				VariantUtil::addInfoValue(it->lineBuffer, "RCN", std::to_string(rcNormal));
				VariantUtil::addInfoValue(it->lineBuffer, "RCT", std::to_string(rcTumor));

				out << it->lineBuffer << std::endl;
			}

			std::for_each(matching.begin(), matching.end(), [](auto* iterator){ iterator->advance(); });

			matching.clear();
			iterators.erase(std::remove_if(iterators.begin(), iterators.end(), [](const auto& item){
				return !item.isGood();
			}), iterators.end());

		}
	}

	static inline void filter(const std::string& file, const std::string& vc, float maxVafNormal, float minVafTumor,
							  const std::string& outFile){

		std::ofstream out(outFile);
		if(out.fail()) throw std::runtime_error("Could not open " + outFile);

		LOG_DEBUG("Filtering " + file);

		VcfLineIterator vcf(file, vc, maxVafNormal, minVafTumor);

		out << vcf.header << std::endl;

		while(vcf.isGood()){
			out << vcf.lineBuffer << std::endl;
			vcf.advance();
		}
	}




	static inline void difference(const std::string& file1, const std::string& file2, const std::string& outFile){

		std::ofstream out(outFile);
		if(out.fail()) throw std::runtime_error("Could not open " + outFile);

		LOG_DEBUG("Creating difference of " + file1 + " and " + file2);

		VcfLineIterator origin(file1, "None", 0.0f, 0.0f);
		VcfLineIterator excluded(file2, "None", 0.0f, 0.0f);

		// Print header
		out << origin.header << std::endl;

		// Combine Variants
		while(origin.isGood() && !excluded.isGood()){
			if(origin < excluded){
				out << origin.lineBuffer << std::endl;
				origin.advance();
			}else if(excluded < origin){
				excluded.advance();
			}else{
				excluded.advance();
				origin.advance();
			}
		}

		while(origin.isGood()){
			out << origin.lineBuffer << std::endl;
			origin.advance();
		}
	}

	static inline void mergeGermlineSomatic(const std::string& germline, const std::string& somatic, const std::string& out) {

		static unsigned long variantId = 1;

		std::ofstream outFile(out);
		if(outFile.fail()) throw std::runtime_error("Could not open " + out);

		VcfLineIterator gl(germline, "None", 0.0f, 0.0f);
		VcfLineIterator s(somatic, "None", 0.0f, 0.0f);

		// Print header
		outFile << "##fileformat=VCFv4.3" << std::endl;
		VcfHeader header = VcfHeader::merge({s.header, gl.header});
		header.info.emplace("VT", "##INFO=<ID=VT,Number=1,Type=String,Description=\"Type of Variant\">");
		header.columnHeaders = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCombined";
		outFile << header << std::endl;

		// Print variants
		while(gl.isGood() && s.isGood()){

			std::string variant = s < gl ? s.lineBuffer : gl.lineBuffer;

			std::vector<size_t> tabs;
			size_t pos = variant.find('\t');
			while( pos != std::string::npos) {
				tabs.push_back(pos);
				pos = variant.find('\t', pos + 1);
			}


			if(s < gl){
				// remove normal sample
				variant.erase(tabs[8]+2, tabs[9]-tabs[8]);
				variant.insert(tabs[7], ";VT=SOMATIC");
				s.advance();
			}else{
				variant.insert(tabs[7], ";VT=GERMLINE");
				gl.advance();
			}

			// add id
			variant.replace(tabs[1]+1, tabs[2]-tabs[1]-1, std::to_string(variantId++));
			outFile << variant << std::endl;
		}

		// add rest
		while(gl.isGood()){
			std::string variant = gl.lineBuffer;

			std::vector<size_t> tabs;
			size_t pos = variant.find('\t');
			while( pos != std::string::npos) {
				tabs.push_back(pos);
				pos = variant.find('\t', pos + 1);
			}

			// add id
			variant.insert(tabs[7], ";VT=GERMLINE");
			variant.replace(tabs[1]+1, tabs[2]-tabs[1]-1, std::to_string(variantId++));

			gl.advance();

			outFile << variant << std::endl;
		}

		while(s.isGood()){
			std::string variant = s.lineBuffer;

			std::vector<size_t> tabs;
			size_t pos = variant.find('\t');
			while( pos != std::string::npos) {
				tabs.push_back(pos);
				pos = variant.find('\t', pos + 1);
			}

			variant.erase(tabs[8]+2, tabs[9]-tabs[8]);
			variant.insert(tabs[7], ";VT=SOMATIC");
			variant.replace(tabs[1]+1, tabs[2]-tabs[1]-1, std::to_string(variantId++));
			s.advance();

			outFile << variant << std::endl;
		}

	}

	static inline void splitZygosity(const std::string& combined, const std::string& prefix){

		const std::string fileHomozygous = prefix + ".homozygous.vcf";
		const std::string fileHeterozygous = prefix + ".heterozygous.vcf";

		std::ofstream outHomozygous(fileHomozygous);
		std::ofstream outHeterozygous(fileHeterozygous);

		if(outHomozygous.fail()) throw std::runtime_error("Could not open " + fileHomozygous);
		if(outHeterozygous.fail()) throw std::runtime_error("Could not open " + fileHeterozygous);

		LOG_DEBUG("Splitting " + combined + " into " + fileHomozygous + " | " + fileHeterozygous);

		VcfLineIterator vcf(combined, "None", 0.0f, 0.0f);

		outHomozygous << vcf.header << '\n';
		outHeterozygous << vcf.header << '\n';

		while(vcf.isGood()){

			bool isHomozygous = vcf.lineBuffer.find("\t1/1") != std::string::npos ||
								vcf.lineBuffer.find("\t1|1") != std::string::npos;

			if(isHomozygous){
				outHomozygous << vcf.lineBuffer << '\n';
			}else{
				outHeterozygous << vcf.lineBuffer << '\n';
			}

			vcf.advance();
		}


	}
}

#endif //EPITOPE_VCFUTIL_H
