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

#if !defined(NCL_OUTPUT_OUTPUTSTREAM_H)
#define NCL_OUTPUT_OUTPUTSTREAM_H

#include "phycas/src/ncl/nxs_defs.hpp"
#include "boost/shared_ptr.hpp"

typedef std::pair<std::string, std::string> NxsSAXAttribute;
typedef std::vector<NxsSAXAttribute> VecNxsSAXAttribute;
template<class OUTSTREAM>
const std::string & EscapeMessageForOutput(std::string & message);
template<class OUTSTREAM>
const std::string GetEscapedMessageForOutput(const std::string  & message);

class NxsXMLSocketOutputStream;
const std::string & EscapeXMLMessage(std::string & message);

template<>
inline const std::string & EscapeMessageForOutput<NxsXMLSocketOutputStream>(std::string & message)
	{
	return EscapeXMLMessage(message);
	}

template<class OUTSTREAM>
inline const std::string & EscapeMessageForOutput(std::string & message)
	{
	return message;
	}

template<class OUTSTREAM>
const std::string GetEscapedMessageForOutput(const std::string  & message)
	{
	std::string escapedMessage(message);
	return EscapeMessageForOutput<OUTSTREAM>(escapedMessage);
	}

#if defined(NCL_USE_STD_OUTPUT)
#	include <iostream>
	namespace ncl
	{
	inline NxsWritableStream & endl(NxsWritableStream & outStream)
		{
		return std::endl(outStream);
		}

	inline NxsWritableStream & flush(NxsWritableStream & outStream)
		{
		return std::flush(outStream);
		}

	inline void Prompt(NxsWritableStream & outStream, const char * w)
		{
		outStream << '\n' << w << std::flush;
		}

	} // namespace ncl

#else //defined(NCL_USE_STD_OUTPUT)
#	if defined(NCL_USE_NXS_CONSOLE_OUTPUT)
#		include "ncl/output/nxs_console_out_stream_impl.hpp"
#	else // defined(NCL_USE_NXS_CONSOLE_OUTPUT)
#		include NCL_USER_OUTPUT_HEADER
#	endif // defined(NCL_USE_NXS_CONSOLE_OUTPUT)
#endif //defined(NCL_USE_STD_OUTPUT)

#include "phycas/src/ncl/output/nxs_output_stream_wrapper.hpp"


#if defined(NCL_USE_NXS_CONSOLE_OUTPUT) || defined(NCL_USE_STD_OUTPUT)
	class NxsPlotStream: public NxsWritableStream
		{
		};
#endif

typedef boost::shared_ptr<NxsWritableStream> NxsWritableStreamShPtr;
// nxs_sax_output_wrapper will be useful whenever, NxsXMLSocketOutputStream is used (so we don't want to have
// to explicitly include it every time, but it requires the full nxs_output_stream.hpp header, so we put this include at the bottom
class NxsOutFilePath;
class NxsOutputDestinationDescription;
typedef std::vector<NxsWritableStream *> VecStreamPtrs;
		
class NxsOutputStreamWrapper: public boost::noncopyable
	{
	public:
		NxsOutputStreamWrapper(const NxsOutputDestinationDescription & nodd, NxsOutputManager & );
		~NxsOutputStreamWrapper();
		std::ofstream * const GetFilePtr() const
			{
			return filePtr;
			}
		const UInt GetNumNCLStreamPtr() const
			{
			return (UInt) streamPtrs.size();
			}
		NxsWritableStream * const GetNCLStreamPtr(unsigned i) const
			{
			return (i >= streamPtrs.size() ? NULL : streamPtrs[i]);
			}
		NxsPlotStream * GetPlotStreamPtr() const
			{
			return plotPtr;
			}
		void SetPlotStreamPtr(NxsPlotStream *p)
			{
			plotPtr = p;
			}
		operator bool () const // implicit cast to bool
			{
			return (streamPtrs.empty()) || (filePtr != NULL) || (plotPtr != NULL);
			}
		VecStreamPtrs		streamPtrs;
		
	private:
		void	ResetFile(NxsOutFilePath & nofp);
		NxsPlotStream     * plotPtr;
		std::ofstream     * filePtr;
		friend NxsOutputStreamWrapper & ncl::endl(NxsOutputStreamWrapper & outStream);
		friend NxsOutputStreamWrapper & ncl::flush(NxsOutputStreamWrapper & outStream);
	};
namespace ncl
{
inline NxsOutputStreamWrapper & endl(NxsOutputStreamWrapper & outStream)
	{
	for (VecStreamPtrs::iterator sIt = outStream.streamPtrs.begin(); sIt != outStream.streamPtrs.end(); ++sIt)
		ncl::endl(**sIt);
	if (outStream.GetFilePtr() != NULL)
		std::endl(*outStream.GetFilePtr());
#   if !defined(NCL_USE_NXS_CONSOLE_OUTPUT) && ! defined(NCL_USE_STD_OUTPUT)
		if (outStream.GetPlotStreamPtr())
			ncl::endl(*outStream.GetPlotStreamPtr());
#   endif
	return outStream;
	}

inline NxsOutputStreamWrapper & flush(NxsOutputStreamWrapper & outStream)
	{
	for (VecStreamPtrs::iterator sIt = outStream.streamPtrs.begin(); sIt != outStream.streamPtrs.end(); ++sIt)
		ncl::flush(**sIt);
	if (outStream.GetFilePtr() != NULL)
		std::flush(*outStream.GetFilePtr());
#   if !defined(NCL_USE_NXS_CONSOLE_OUTPUT) && ! defined(NCL_USE_STD_OUTPUT)
		if (outStream.GetPlotStreamPtr())
			ncl::flush(*outStream.GetPlotStreamPtr());
#   endif
	return outStream;
	}

} // namespace ncl

inline NxsOutputStreamWrapper & operator<<(NxsOutputStreamWrapper & outStream, NxsOutputStreamWrapperManipPtr funcPtr)
	{
	return (*funcPtr)(outStream);
	}

inline  NxsOutputStreamWrapper & operator<<(NxsOutputStreamWrapper & outStream, int i) 
	{
	for (VecStreamPtrs::iterator sIt = outStream.streamPtrs.begin(); sIt != outStream.streamPtrs.end(); ++sIt)
		**sIt << i;
	if (outStream.GetFilePtr() != NULL)
		*outStream.GetFilePtr() << i;
	return outStream;
	}
	
inline  NxsOutputStreamWrapper & operator<<(NxsOutputStreamWrapper & outStream, long l) 
	{
	for (VecStreamPtrs::iterator sIt = outStream.streamPtrs.begin(); sIt != outStream.streamPtrs.end(); ++sIt)
		**sIt << l;
	if (outStream.GetFilePtr() != NULL)
		*outStream.GetFilePtr() << l;
	return outStream;
	}
	
inline  NxsOutputStreamWrapper & operator<<(NxsOutputStreamWrapper & outStream, unsigned i) 
	{
	for (VecStreamPtrs::iterator sIt = outStream.streamPtrs.begin(); sIt != outStream.streamPtrs.end(); ++sIt)
		**sIt << i;
	if (outStream.GetFilePtr() != NULL)
		*outStream.GetFilePtr() << i;
	return outStream;
	}
		
inline  NxsOutputStreamWrapper & operator<<(NxsOutputStreamWrapper & outStream, double d) 
	{
	for (VecStreamPtrs::iterator sIt = outStream.streamPtrs.begin(); sIt != outStream.streamPtrs.end(); ++sIt)
		**sIt << d;
	if (outStream.GetFilePtr() != NULL)
		*outStream.GetFilePtr() << d;
	return outStream;
	}
	
inline  NxsOutputStreamWrapper & operator<<(NxsOutputStreamWrapper & outStream, const char *c) 
	{
	for (VecStreamPtrs::iterator sIt = outStream.streamPtrs.begin(); sIt != outStream.streamPtrs.end(); ++sIt)
		**sIt << c;
	if (outStream.GetFilePtr() != NULL)
		*outStream.GetFilePtr() << c;
	return outStream;
	}
	
inline  NxsOutputStreamWrapper & operator<<(NxsOutputStreamWrapper & outStream, char c) 
	{
	for (VecStreamPtrs::iterator sIt = outStream.streamPtrs.begin(); sIt != outStream.streamPtrs.end(); ++sIt)
		**sIt << c;
	if (outStream.GetFilePtr() != NULL)
		*outStream.GetFilePtr() << c;
	return outStream;
	}
	
inline  NxsOutputStreamWrapper & operator<<(NxsOutputStreamWrapper & outStream, const std::string &s) 
	{
	for (VecStreamPtrs::iterator sIt = outStream.streamPtrs.begin(); sIt != outStream.streamPtrs.end(); ++sIt)
		**sIt << s;
	if (outStream.GetFilePtr() != NULL)
		*outStream.GetFilePtr() << s;
	return outStream;
	}

template<class OUTSTREAM>
void SetNewOutputContext(OUTSTREAM & outStream, const VecNxsSAXAttribute & attributes);
void SetStreamAttributes(NxsXMLSocketOutputStream & outStream, const VecNxsSAXAttribute & attributes);

template<class OUTSTREAM>
void SetNewOutputContext(OUTSTREAM & /*outStream*/, const VecNxsSAXAttribute & /*attributes*/)
	{
	// nothing to be done in most output stream styles (only used in rich output stream types, like xml socket
	}

template<>
inline void SetNewOutputContext<NxsXMLSocketOutputStream>(NxsXMLSocketOutputStream & outStream, const VecNxsSAXAttribute & attributes)
	{
	SetStreamAttributes(outStream, attributes);
	}

#endif // !defined(NCL_OUTPUT_OUTPUTSTREAM_H)

