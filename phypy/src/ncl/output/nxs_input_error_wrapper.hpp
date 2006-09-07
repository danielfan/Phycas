#ifndef NXS_INPUT_ERROR_WRAPPER_HPP
#define NXS_INPUT_ERROR_WRAPPER_HPP
#include <string>
#include <stack>
#include "ncl/nxs_defs.hpp"
#include "ncl/misc/nxs_file_path.hpp"
#include "ncl/output/temp_basic_output_operators.hpp"

class NxsInputErrorWrapper
	{
	public:
		NxsInputErrorWrapper(
		  const std::string inMsg, 
		  const std::string & inCurrBlockID, 
		  const std::stack<NxsInFilePath> & inFileStack,  
		  file_pos inPos, 
		  unsigned inLineNum, 
		  unsigned inColNum)
		  	:message(inMsg),
		  	blockID(inCurrBlockID),
		  	fileStack(inFileStack),
		  	pos(inPos),
		  	line(inLineNum),
		  	col(inColNum)
		  	{
		  	}
	
		std::string message;
		std::string blockID; 
		std::stack<NxsInFilePath> fileStack;  
		file_pos pos; 
		unsigned line;
		unsigned col;
	};

template<int FORMAT_HINT, class OUT_STREAM>
class GenericPrinterClass<FORMAT_HINT, NxsInputErrorWrapper, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & errorStream, const NxsInputErrorWrapper & ies)
			{ 
			errorStream << "Error";
			if (!ies.fileStack.empty())
				{
				if (!ies.blockID.empty())
					errorStream << " [in a " << ies.blockID << " block";
				if (!ies.fileStack.empty())
					{
					errorStream << " in the file " << ies.fileStack.top().GetFileName();
					if (ies.line > 0)
						errorStream << " at line " << ies.line << ", column " << ies.col;
					if (ies.fileStack.size() > 1)
						{
						std::stack<NxsInFilePath> temp = ies.fileStack;
						temp.pop();
						while (!temp.empty())
							{
							errorStream << "\n\t(which was loaded by the file " << temp.top().GetFileName() << ')';
							temp.pop();
							}
						}
					}
					errorStream << ']';
				}
			if (ies.message.empty())
				errorStream << '.';
			else
				errorStream << ":\n\t" << ies.message;
			}
	};
#endif //NXS_GENERIC_OUTPUT_CLASSES_HPP