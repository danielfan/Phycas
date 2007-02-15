/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//#include "phycas/force_include.h"
#include "phycas/src/ncl/output/generic_output_classes.hpp"
#include "phycas/src/ncl/output/temp_basic_output_operators.hpp"
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
