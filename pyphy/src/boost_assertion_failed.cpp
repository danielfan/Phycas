#if !defined(NDEBUG)
#include <boost/assert.hpp>

void boost::assertion_failed(char const *expr, char const *function, char const *file, long line)
	{
	std::cerr << "\nBoost assertion failed:";
	std::cerr << "\n  expr: " << expr;
	std::cerr << "\n  func: " << function;
	std::cerr << "\n  file: " << file;
	std::cerr << "\n  line: " << line;
	std::cerr << std::endl;
	}
#endif
