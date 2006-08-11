//#include "phycas/force_include.h"
#include "pyphy/src/ncl/output/generic_output_classes.hpp"
#include "pyphy/src/ncl/output/temp_basic_output_operators.hpp"
#if defined(MTH_PHYCAS_MODULARIZING_OUTPUT) && MTH_PHYCAS_MODULARIZING_OUTPUT && 0
using std::stack;

class NxsXMLSocketOutputStream;

template<>
class GenericPrinter<NxsXMLSocketOutputStream, NxsInputErrorWrapper>
	{
	public:
		NxsXMLSocketOutputStream & Print(NxsXMLSocketOutputStream & stream, const NxsInputErrorWrapper & ies);
	};

	

template<class OUT_STREAM>
class  GenericPrinter<OUT_STREAM, NxsInputErrorWrapper>
	{
	public:
	OUT_STREAM & Print(OUT_STREAM & errorStream, const NxsInputErrorWrapper & ies)
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
					stack<NxsInFilePath> temp = ies.fileStack;
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
		return errorStream;
		}
	};


NxsErrorStream & operator<<(NxsErrorStream & errorStream, const NxsInputErrorWrapper & ies)
	{
	GenericPrinter<NxsErrorStream, NxsInputErrorWrapper> printer;
	return printer.Print(errorStream, ies);
	}

#include "ncl/output/nxs_xml_socket_output_stream.hpp"

NxsXMLSocketOutputStream & GenericPrinter<NxsXMLSocketOutputStream, NxsInputErrorWrapper>::Print(NxsXMLSocketOutputStream & stream, const NxsInputErrorWrapper & ies)
	{
	stream << "\nError found";
	if (!ies.fileStack.empty)
		{
		if (!ies.currBlockID.empty())
			stream << " in a " << ies.currBlockID << " block";
			{
			stream << " in the file " << ies.fileStack.top().GetFileName();
			if (ies.line > 0)
				stream << " at line " << ies.line << ", column " << ies.col;
			if (ies.fileStack.size() > 1)
				{
				stack<NxsInFilePath> temp = fileStack;
				temp.pop();
				while (!temp.empty())
					{
					stream << "\n\t(which was loaded by the file " << temp.top().GetFileName() << ')';
					temp.pop();
					}
				}
			}
		}
	if (ies.msg.empty())
		stream << '.';
	else
		stream << ":\n" << ies.msg;
	return stream;
	}


#endif //if defined(MTH_PHYCAS_MODULARIZING_OUTPUT) && MTH_PHYCAS_MODULARIZING_OUTPUT
