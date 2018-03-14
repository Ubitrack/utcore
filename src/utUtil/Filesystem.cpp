#include <utUtil/Filesystem.h>

#include <boost/regex.hpp>

namespace Ubitrack { namespace Util {

std::string matchToEnv(boost::smatch const& match) {
  const char * s = getenv( match[1].str().c_str() );
  return s == NULL ? std::string("") : std::string(s);
}

std::string expandEnvironmentVariables(const std::string& in) {
  static boost::regex re( "\\$\\{([^}]+)\\}" );
  return boost::regex_replace(in, re, matchToEnv);
}

boost::filesystem::path getFilesystemPath(const std::string& in) {
	std::string expanded(expandEnvironmentVariables(in));
	boost::filesystem::path p(expanded);
	return p;
}

}} // Ubitrack::Util