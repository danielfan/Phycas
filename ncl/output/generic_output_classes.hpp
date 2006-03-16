#ifndef NXS_GENERIC_OUTPUT_CLASSES_HPP
#define NXS_GENERIC_OUTPUT_CLASSES_HPP
#if defined(MTH_PHYCAS_MODULARIZING_OUTPUT) && MTH_PHYCAS_MODULARIZING_OUTPUT

#include <string>
#include <stack>
#include "ncl/nxs_defs.hpp"
#include "ncl/misc/nxs_file_path.hpp"
#include "ncl/output/temp_basic_output_operators.hpp"

class InputErrorSpec
	{
	public:
		InputErrorSpec(
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
class GenericPrinterClass<FORMAT_HINT, InputErrorSpec, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & errorStream, const InputErrorSpec & ies)
			{ 
			errorStream << "Error:\n";
			if (!ies.fileStack.empty())
				{
				if (!ies.blockID.empty())
					errorStream << " in a " << ies.blockID << " block";
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
				}
			if (ies.message.empty())
				errorStream << '.';
			else
				errorStream << ":\n" << ies.message;
			}
	};
# if 0
class Tree;	
NxsErrorStream & operator<<(NxsErrorStream & outS, const InputErrorSpec & s);

enum NCLSetOutStyle
		{
		kNoTranslateLabels	= 0x00,
		kTranslateLabels	= 0x00
		};

class NxsBasicListManager;
void PrintNamedSet(NxsWritableStream * outStream, const NxsIndexSet & desertTaxa, const NxsBasicListManager * const mgr = NULL, NCLSetOutputType setType = kGenericSetOutput, NCLOutputFormatHint fmt = kExpanded, NCLSetOutStyle style = kTranslateLabels);

enum NCLTreeOutStyle
	{
	kGraphicalTreeRep,
	kNewickTreeRep
	};

template<class OUTSTREAM>
void PrintTree(OUTSTREAM *, const Tree & , NCLTreeOutStyle = kNewickTreeRep);
#endif
#endif// defined(MTH_PHYCAS_MODULARIZING_OUTPUT) && MTH_PHYCAS_MODULARIZING_OUTPUT
#endif //NXS_GENERIC_OUTPUT_CLASSES_HPP

