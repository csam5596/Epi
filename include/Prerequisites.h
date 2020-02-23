#ifndef EPI_PREREQUISITES_H
#define EPI_PREREQUISITES_H

#include <string>

namespace Prerequisites {


	static inline bool check(){

		const static std::array<std::string, 8> programs = {"netMHCpan", "python", "bgzip", "tabix", "vep", "grep", "samtools", "xargs"};

		bool check = true;
		std::for_each(programs.cbegin(), programs.cend(), [&check](const std::string& program){

			std::string command = "which " + program + " > /dev/null";

			bool isInstalled = (system(command.c_str()) == 0);
			check = (check && isInstalled);

			if(!isInstalled){
				LOG_ERROR(program + " executable not found");
			}
		});

		return check;
	}

}

#endif //EPI_PREREQUISITES_H
