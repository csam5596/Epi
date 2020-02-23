#ifndef EPITOPE_NETMHCPAN_H
#define EPITOPE_NETMHCPAN_H

#include <string>
#include <boost/process.hpp>
#include <boost/format.hpp>
#include <hla/HlaTypes.h>
#include <Logging.h>

namespace NetMhcPan{

	static inline void call(const std::string& peptideFile, const HlaTypes& hla){

		namespace bp = boost::process;

		LOG_DEBUG("Calling netMHCpan...");

		std::ofstream netMhcPanOut("log/netMHCpan.out");
		if(netMhcPanOut.fail()) throw std::runtime_error("Could not open file log/netMHCpan.out");

		bp::ipstream outStream;
		static auto commandFormat = boost::format("netMHCpan -f %1% -p -inptype %2% -BA -a %3%");
		const std::string netMHCpanCommand = boost::str(commandFormat % peptideFile % 1 % hla.toString());
		bp::child netMHCPan(netMHCpanCommand, bp::std_out > outStream, bp::std_err > "log/netMHCpan.err");

		std::string line;
		while(outStream && getline(outStream, line)){
			netMhcPanOut << line << std::endl;
			if(line.find("Number of high binders") != std::string::npos) {
				LOG_DEBUG(line);
			}
		}
		netMHCPan.wait();
		if(netMHCPan.exit_code() != 0){
			LOG_ERROR("[" + netMHCpanCommand + "] returned " + std::to_string(netMHCPan.exit_code()));
			exit(1);
		}
	}


}

#endif //EPITOPE_NETMHCPAN_H
