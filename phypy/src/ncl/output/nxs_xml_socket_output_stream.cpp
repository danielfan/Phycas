/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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
#if defined (NCL_SOCKET_IO)
#	include "gui/phycasGUI/network/PhycasSocket.h"
#	if defined (WIN_PHOREST)
#		include <process.h>
#	endif
#	include <string>
#   include <list>
#   include <ctype.h> 
#   include "ncl/nxs_defs.hpp"
#	include "ncl/misc/algorithm_extensions.hpp"
#	include "ncl/nxs_token.hpp"
#	include "ncl/nxs_basic_manager.hpp"
#	include "ncl/misc/nxs_index_set.hpp"
#	include "ncl/output/nxs_output_stream.hpp"
#   include "ncl/output/nxs_output.hpp"
	
	using std::string;
	using std::pair;
	using std::vector;
	using ncl::endl;
	std::string EscapeXMLAttribute(const std::string &att);
	/// \bug not implemented correctly (currently just returns the argument, unescaped)
	std::string EscapeXMLAttribute(const std::string &att)
		{
		return att;
		}
	
	/// Flushes previous message if not HasContent() is true and then appends the attributes to the outstreams
	/// base level tag (does NOT check for duplication of attributes).
	/// \bug attribute values are not escaped - see 	EscapeXMLAttribute()
	void SetStreamAttributes(NxsXMLSocketOutputStream & outStream, const VecNxsSAXAttribute & attributes)
		{
		if (outStream.HasContent())
			outStream.PrintMessage();
		for (VecNxsSAXAttribute::const_iterator attPair = attributes.begin(); attPair != attributes.end(); ++attPair)
			outStream.currOpenTag << ' ' << attPair->first << "=\"" << EscapeXMLAttribute(attPair->second) << '\"';
		outStream.ResetMessage();
		}

	NxsXMLSocketOutputStream & NxsXMLSocketOutputStream::StartElementQName(const std::string &  qName, const VecNxsSAXAttribute & atts)
		{
		string e = qName;
		for (VecNxsSAXAttribute::const_iterator aIt = atts.begin(); aIt != atts.end(); ++aIt)
			e << ' ' << aIt->first << "=\"" << aIt->second << '\"';
		return OpenElement(e);
		}

	
	NxsXMLSocketOutputStream & NxsXMLSocketOutputStream::AddObject(const std::string & inOpenTag, const std::string & content, const std::string & closeTagName)
		{
		OpenElement(inOpenTag);
		*this << content;
		return CloseElement(closeTagName);
		}
		
	/*
	void NxsXMLSocketOutputStream::SetQueryOrPromptType(NCLQueryOrPromptType t)
			{
			if (t == NCLQueryOrPromptType(kGenericOutput))
				{
				if (GetQueryOrPromptType(*this) != NCLQueryOrPromptType(kGenericOutput))
					PrintMessage();
				msgType = t;
				}
			if (GetQueryOrPromptType(*this) != t)
				{
				if (HasContent())
					PrintMessage();
				// everything that is not kGenericOutput is expected to be xml
				msgType = t;
				SetIsXML();
				}
			ResetMessage();
			}
	*/				
/*		if (msgType == kGenericOutput)
			return ;
		std::string retStr = startOpenTag + " type=\"";
		switch (msgType)
			{
			case kPromptForFileChooser: 
				retStr.append("file");
				break;
			case kPromptForAnyString:
				retStr.append("string");
				break;
			case kPromptAlert:
				retStr.append("alert");
				break;
			case kPromptChoices:
				retStr.append("choices");
				break;
			case kPromptCancelOK:
				retStr.append("cancel_ok");
				break;
			case kPromptNoYes:
				retStr.append("no_yes");
				break;
			default: 
				break;
			}
		retStr << "\">";
		return retStr;
		}
*/
	NxsXMLSocketOutputStream::NxsXMLSocketOutputStream(const std::string & t, const std::string prePromptOutput)
		:dblFormat(UINT_MAX, 6),
		startOpening("<" + t),
		currOpenTag("<" + t),
		prePromptString(prePromptOutput),
		closeTag("</" + t + ">"),
		_isXML(false),
		needsCloseTag(false)
		{
		ResetMessage();
		SetWidth(kDefaultPrintWidth);
		}

	const std::string NxsXMLSocketOutputStream::GetOpenTag() const
		{
		return currOpenTag + '>';
		}
	void NxsXMLSocketOutputStream::DisplaySet(const NxsIndexSet &s, bool useNumbers, const NxsIndexInfo &indInfo)
		{
		*this << s.name << " = ";
		if (useNumbers)
				*this << s.GetNexusDescription(true);
		else
			{
			for (NxsIndexSet::const_iterator sIt = s.begin(); sIt != s.end(); ++sIt)
				*this << indInfo.GetLabel(*sIt) << ' ';	
			}
		*this << ncl::flush;		
		}

	void NxsXMLSocketOutputStream::ConvertMessageToXML()
		{
		if (!_isXML)
			{
			const string em = GetEscapedMessage();
			ResetMessage();
			storedMessage += em;
			_isXML = true;
			}
		}
		
	void NxsXMLSocketOutputStream::PrintMessage()
		{
		ConvertMessageToXML();
		storedMessage << GetCloseTag();
		RawPrint();	
		}
		
	void NxsXMLSocketOutputStream::RawPrint()
		{
		NxsOutputManager & outMgr = NxsOutputManager::GetInstance();
		if(outMgr.IsOpen()) 
			{
			std::cout << storedMessage << std::endl;
			outMgr.SendLine(storedMessage);
			}
		else
			std::cout << "socket null" << std::flush;
		ResetMessage();
		_isXML = false;
		needsCloseTag = false;
		currOpenTag = startOpening;
		}
		
	void NxsXMLSocketOutputStream::Prompt()
		{
		ConvertMessageToXML();
		PrintMessage();
		storedMessage = prePromptString;
		RawPrint();
		ResetMessage();
		}
		
	const std::string & EscapeXMLMessage(std::string & message)
		{
		replace_all_substr(message, "&", "&amp;");
		replace_all_substr(message, "<", "&lt;");
		replace_all_substr(message, ">", "&gt;");
		return message;
		}
		
	void NxsXMLSocketOutputStream::PrintTable(NxsTable *p)
		{
		PHYCAS_ASSERT(p == &table);	//you should be getting the Table pointer from NxsXMLSocketOutputStream::GetTablePtr so that it has the correct width etc
		try 
			{
			p->Show(*this);
			p->Reset();
			}
		catch(NxsTable::NxsX_InsufficientWidth &)
			{
			*this << "output width is too small" << ncl::endl;
			}
		}

		
#endif // defined(NCL_SOCKET_IO)
