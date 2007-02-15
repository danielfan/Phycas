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

#ifndef NCL_NXS_XML_SOCKET_OUTSTREAM_H
#define NCL_NXS_XML_SOCKET_OUTSTREAM_H

#include "phycas/src/ncl/nxs_defs.hpp"
#include <iostream>
#include <boost/noncopyable.hpp>
#include "phycas/src/ncl/output/nxs_table_cell.hpp"
#include "phycas/src/ncl/output/nxs_table.hpp"
#include "phycas/src/ncl/misc/string_extensions.hpp"

#if defined (SUPPORT_TREE_OUTPUT)
	class Tree;
#endif

class NxsIndexInfo;
class NxsInputErrorWrapper;
/*----------------------------------------------------------------------------------------------------------------------
|	Provides "put-to" operator interface for output of a temporary (and not very rich) xml format to a socket.
|   Under development for interaction with GUI
|	calling NxsXMLSocketOutputStream << endl causes the output stream to be flused.   
*/
class NxsXMLSocketOutputStream NON_COPYABLE
	{
	public :
		STATIC_CONST const unsigned kDefaultPrintWidth = 80;
		
		
		explicit NxsXMLSocketOutputStream(const std::string &xmlTagName, const std::string prePrompterOutput = std::string());
		
		NxsXMLSocketOutputStream & AddObject(const std::string & openTag, const std::string & content, const std::string & closeTagName);
		NxsXMLSocketOutputStream & OpenElement(const std::string & openTag);
		NxsXMLSocketOutputStream & CloseElement(const std::string & closeTag);

		unsigned	GetWidth() const;
		NxsTable  * GetTablePtr()
			{
			return &table;
			}
		void 		DirectOutputToScreen(bool displayOutput = true);
		void 		DirectOutputToFile(std::ostream *newLogFile);
		void		DisplaySet(const NxsIndexSet &s, bool useNumbers, const NxsIndexInfo &);
		void		Prompt() ;
		void 		PrintMessage();
		void 		PrintTable(NxsTable * t);
		void 		SetNDigitsAfterDecimal(unsigned i);
		void		SetWidth(unsigned outputWidth);
		void 		StopLoggingToFile();
		void SetIsXML(bool isXML = true)
			{
			if (!_isXML && isXML && HasContent())
				ConvertMessageToXML();
			else
				_isXML = isXML;
			}

		std::string			storedMessage; 			/* current stream (output that has not been printed) */
		DblFormatter		dblFormat;			/* class used to format doubles in a consistent manner */
		
		//  Basic (subset of the)  ContentHandler functionality from the SAX API
		//
		NxsXMLSocketOutputStream & Characters(const std::string &c)
			{
			std::string s(c);
			storedMessage.append(EscapeXMLMessage(s));
			return *this;
			}
			
		NxsXMLSocketOutputStream & EndElement(const std::string & uri, const std::string &  name, const std::string &  qName)
			{
			return (uri.empty() ? EndElementQName(qName) : EndElementQName(uri + name));
			}
		NxsXMLSocketOutputStream & StartElement (const std::string & uri, const std::string &  name, const std::string &  qName, const VecNxsSAXAttribute & atts)
			{
			return (uri.empty() ? StartElementQName(qName, atts) : StartElementQName(uri + name, atts));
			}
		NxsXMLSocketOutputStream & EndElementQName(const std::string &  qName)
			{
			return CloseElement(qName);
			}
		NxsXMLSocketOutputStream & StartElementQName(const std::string &  qName, const VecNxsSAXAttribute & atts);

	protected :
		const std::string & GetPrePrompt() const
			{
			return prePromptString;
			}
		void				ConvertMessageToXML();
		const std::string   GetOpenTag() const;
		const std::string GetEscapedMessage() const
			{
			if (!HasContent())
				return std::string();
			std::string content(storedMessage.begin() + static_cast<std::string::difference_type>(beginContentIndex), storedMessage.end());
			return EscapeXMLMessage(content);
			}
		const std::string & GetCloseTag() const
			{
			return closeTag;
			}
		bool HasContent() const
			{
			return (storedMessage.length() > beginContentIndex);
			}
		void				ResetMessage();
		void				PrintLine(const std::string &line) const;
		void				RawPrint();
		const std::string   startOpening;
		std::string			currOpenTag;
		const std::string   prePromptString;
		const std::string 	closeTag;
		std::size_t			width;
		NxsTable			table;		/* for use in formatting tabular output */
		 // the following members not persistent - they revert when PrintMessage is called
		bool				_isXML; /* true if the current message is already xml (and does not need wrapping) - not persistent*/
		bool				needsCloseTag;
		unsigned			beginContentIndex;
		friend void SetStreamAttributes(NxsXMLSocketOutputStream & outStream, const VecNxsSAXAttribute & attributes);
	};
	
typedef NxsXMLSocketOutputStream & (*NxsXMLSocketOutputStreamManipPtr)(NxsXMLSocketOutputStream &); // typedef the outputstream manipulator used for calls such as << endl

inline void  NxsXMLSocketOutputStream::ResetMessage()
	{
	storedMessage = GetOpenTag();
	beginContentIndex = (UInt) storedMessage.length();
	}
	
inline void	NxsXMLSocketOutputStream::SetWidth(
  unsigned outputWidth)
	{
	width = outputWidth;
	}
	
inline unsigned	NxsXMLSocketOutputStream::GetWidth() const
	{
	return (unsigned)width;
	}
	
NxsXMLSocketOutputStream & operator<<(NxsXMLSocketOutputStream & , NxsXMLSocketOutputStreamManipPtr);
//NxsXMLSocketOutputStream & endl(NxsXMLSocketOutputStream& om);
inline NxsXMLSocketOutputStream  &operator<<(NxsXMLSocketOutputStream & out, NxsXMLSocketOutputStreamManipPtr funcPtr)
	{
	return (*funcPtr)(out);
	}

inline  NxsXMLSocketOutputStream &operator<<(NxsXMLSocketOutputStream & out, int i) 
	{
	out.storedMessage << i; 
	return out;
	}
	
inline  NxsXMLSocketOutputStream &operator<<(NxsXMLSocketOutputStream & out, long l) 
	{
	out.storedMessage << l; 
	return out;
	}
	
inline  NxsXMLSocketOutputStream &operator<<(NxsXMLSocketOutputStream & out, unsigned i) 
	{
	out.storedMessage << i; 
	return out;
	}
	
inline  NxsXMLSocketOutputStream &operator<<(NxsXMLSocketOutputStream & out, unsigned long l) 
	{
	out.storedMessage << l; 
	return out;
	}
	
inline  NxsXMLSocketOutputStream &operator<<(NxsXMLSocketOutputStream & out, double d) 
	{
	out.dblFormat.FormatDouble(out.storedMessage, d);
	return out;
	}
	
inline  NxsXMLSocketOutputStream &operator<<(NxsXMLSocketOutputStream & out, const char *c) 
	{
	out.storedMessage << c; 
	return out;
	}
	
inline  NxsXMLSocketOutputStream &operator<<(NxsXMLSocketOutputStream & out, char c) 
	{
	out.storedMessage << c; 
	return out;
	}
	
inline  NxsXMLSocketOutputStream &operator<<(NxsXMLSocketOutputStream & out, const std::string & s) 
	{
	out.storedMessage.append(s);
	return out;
	}
namespace ncl
{
inline  NxsXMLSocketOutputStream & flush(NxsXMLSocketOutputStream& om)
	{
	om.PrintMessage();
	return om;
	}
	
inline  NxsXMLSocketOutputStream & endl(NxsXMLSocketOutputStream& om)
	{
	om << '\n';
	return flush(om);
	}

} // namespace ncl

inline  void NxsXMLSocketOutputStream::SetNDigitsAfterDecimal(unsigned i)
	{
	dblFormat.digitsAfterDecimal = i;
	}

inline NxsXMLSocketOutputStream & NxsXMLSocketOutputStream::OpenElement(const std::string & inOpenTag)
	{
	SetIsXML();
	return *this << '<' << inOpenTag << '>';
	}
	
inline NxsXMLSocketOutputStream & NxsXMLSocketOutputStream::CloseElement(const std::string & closeTag)
	{
	return *this << "</" << closeTag << '>';
	}
	
#if defined(NCL_SOCKET_IO)
	class NxsPlotStream: public NxsXMLSocketOutputStream
		{
		public:
			 NxsPlotStream(const std::string &xmlTagName, const std::string prePrompterOutput = std::string())
				:NxsXMLSocketOutputStream(xmlTagName, prePrompterOutput)
				{}
		};
#endif // defined(NCL_SOCKET_IO)
#endif
