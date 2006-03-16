#if !defined (NCL_NXS_INPUT_HPP)
#define NCL_NXS_INPUT_HPP
#include <deque>
#include "ncl/nxs_defs.hpp"
namespace NxsIO {

class NxsInput
	{
	public:
		void AppendNextLine(std::string *nextLine, std::deque<std::string> *cmdHist);
	};

} // namespace NxsInput
#endif

