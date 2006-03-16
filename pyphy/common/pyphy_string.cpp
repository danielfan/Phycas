#include "pyphy/common/pyphy_string.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|	Appends a string version of the supplied unsigned value `v' to the string `s' and returns a reference to `s'.
*/
std::string &append_unsigned(std::string &s, unsigned v)
	{
	static char tmp[128];
	sprintf(tmp, "%d", v);
	s << tmp;
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Utility function useful for reporting strings in exceptions. If supplied string `s' is longer than `max_chars', it 
|	is shortened to `max_chars' characters and an ellipsis ("...") is appended to the returned value. If `s' is less
|	than `max_chars', it does not need to be abbreviated and is thus returned in its entirety.
*/
std::string abbreviate(const std::string &s, unsigned max_chars)
	{
	if (s.length() < max_chars)
		return s;
	else
		{
		std::string ss(s.substr(0, max_chars));
		ss << "...";
		return ss;
		}
	}

