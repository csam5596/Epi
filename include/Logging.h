#ifndef EPITOPE_LOGGING_H
#define EPITOPE_LOGGING_H

#include <boost/log/trivial.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/attributes/mutable_constant.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/attributes/mutable_constant.hpp>
#include <boost/log/utility/manipulators/add_value.hpp>


#define LOG_DEBUG(MSG) \
	LOG(boost::log::trivial::debug, MSG)

#define LOG_WARN(MSG) \
	LOG(boost::log::trivial::warning, MSG)

#define LOG_ERROR(MSG) \
	LOG(boost::log::trivial::error, MSG)

#define LOG(LEVEL, MSG) \
  BOOST_LOG_SEV(boost::log::trivial::logger::get(), LEVEL) \
    << boost::log::add_value("Line", __LINE__) \
    << boost::log::add_value("File", Logging::filename(__FILE__)) \
    << MSG


namespace Logging {

	static inline void init(){
		boost::log::add_common_attributes();

		boost::log::register_simple_filter_factory<boost::log::trivial::severity_level, char>("Severity");
		boost::log::register_simple_formatter_factory<boost::log::trivial::severity_level, char>("Severity");

		auto syslog_format(
			boost::log::expressions::stream <<
				boost::log::expressions::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S") <<
                "  " << std::setw(8) << std::left << boost::log::trivial::severity <<
                " " << std::setw(20) << std::left << boost::log::expressions::attr<std::string>("File") <<
				" | " << boost::log::expressions::smessage
		);

		auto coutSink = boost::log::add_console_log(std::cout,
                boost::log::keywords::auto_flush = true,
				boost::log::keywords::format = syslog_format
		);

		auto cerrSink = boost::log::add_console_log(std::cerr,
	            boost::log::keywords::auto_flush = true,
				boost::log::keywords::format = syslog_format
		);

		coutSink->set_filter(boost::log::trivial::severity < boost::log::trivial::warning);
		cerrSink->set_filter(boost::log::trivial::severity >= boost::log::trivial::warning);

		auto outFileSink = boost::log::add_file_log(
				boost::log::keywords::auto_flush = true,
				boost::log::keywords::file_name  = "log/epi.log",
				boost::log::keywords::format = syslog_format
		);

		auto errFileSink = boost::log::add_file_log(
				boost::log::keywords::auto_flush = true,
				boost::log::keywords::file_name  = "log/epi.err",
				boost::log::keywords::format = syslog_format
		);

		outFileSink->set_filter(boost::log::trivial::severity < boost::log::trivial::warning);
		errFileSink->set_filter(boost::log::trivial::severity >= boost::log::trivial::warning);

	}

	static inline std::string filename(const std::string& path){
		return path.substr(path.find_last_of("/\\")+1);
	}

}

#endif //EPITOPE_LOGGING_H
