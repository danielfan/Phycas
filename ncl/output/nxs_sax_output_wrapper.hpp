#ifndef NCL_NXS_SAX_OUTPUT_WRAPPER_HPP
#define NCL_NXS_SAX_OUTPUT_WRAPPER_HPP
#include <list>
#include "ncl/output/nxs_output_stream.hpp"
class NxsXMLSocketOutputStream;
template <class OUT_STREAM_TYPE>
class NxsGenericSaxOutputWrapper
	{
	public:
		explicit NxsGenericSaxOutputWrapper(OUT_STREAM_TYPE & outStream, const std::string & elName)
			:outputStream(outStream)
			{ // do nothing with markup for most elements if we are not using NxsXMLSocketOutputStream
			}
			
		explicit NxsGenericSaxOutputWrapper(OUT_STREAM_TYPE & outStream, const std::string & elName, const VecNxsSAXAttribute &atts)
			:outputStream(outStream)
			{// do nothing with markup for most elements if we are not using NxsXMLSocketOutputStream
			}
			
		~NxsGenericSaxOutputWrapper()
			{ // currently we are using NxsGenericSaxOutputWrapper in contexts where each element could get a newline.  This is a hack.
			outputStream << '\n';
			}
	
		OUT_STREAM_TYPE & outputStream;
	};

#include "ncl/output/nxs_xml_socket_output_stream.hpp"
template <>
class NxsGenericSaxOutputWrapper<NxsXMLSocketOutputStream>
	{
	public:
		explicit NxsGenericSaxOutputWrapper(NxsXMLSocketOutputStream & outStream, const std::string & elName)
			:outputStream(outStream),
			tagName(elName)
			{
			outStream.OpenElement(elName);
			}
			
		explicit NxsGenericSaxOutputWrapper(NxsXMLSocketOutputStream & outStream, const std::string & elName, const VecNxsSAXAttribute &atts)
			:outputStream(outStream),
			tagName(elName)
			{
			outStream.StartElementQName(elName, atts);
			}
			
		~NxsGenericSaxOutputWrapper()
			{
			outputStream.CloseElement(tagName);
			}
	
		NxsXMLSocketOutputStream & outputStream;
		const std::string tagName;
	};
	
template <>
class NxsGenericSaxOutputWrapper<NxsOutputStreamWrapper>
	{
	public:
		explicit NxsGenericSaxOutputWrapper(NxsOutputStreamWrapper & outStream, const std::string & elName)
			:outputStream(outStream)
			{
			for (VecStreamPtrs::const_iterator sIt = outStream.streamPtrs.begin(); sIt != outStream.streamPtrs.end(); ++sIt)
				streamWrappers.push_back(NxsSaxOutputWrapper(**sIt, elName));
			if (outStream.GetPlotStreamPtr() != NULL)
				plotWrappers.push_back(NxsPlotSaxOutputWrapper(*outStream.GetPlotStreamPtr(), elName));
			if (outStream.GetFilePtr() != NULL)
				fileWrappers.push_back(NxsFileSaxOutputWrapper(*outStream.GetFilePtr(), elName));
			}
			
		explicit NxsGenericSaxOutputWrapper(NxsOutputStreamWrapper & outStream, const std::string & elName, const VecNxsSAXAttribute &atts)
			:outputStream(outStream)
			{
			for (VecStreamPtrs::const_iterator sIt = outStream.streamPtrs.begin(); sIt != outStream.streamPtrs.end(); ++sIt)
				streamWrappers.push_back(NxsSaxOutputWrapper(**sIt, elName, atts));
			if (outStream.GetPlotStreamPtr() != NULL)
				plotWrappers.push_back(NxsPlotSaxOutputWrapper(*outStream.GetPlotStreamPtr(), elName, atts));
			if (outStream.GetFilePtr() != NULL)
				fileWrappers.push_back(NxsFileSaxOutputWrapper(*outStream.GetFilePtr(), elName, atts));
			}
			
		~NxsGenericSaxOutputWrapper()
			{
			}
	
		NxsOutputStreamWrapper & outputStream;
		// need to create objects whose lifetime is tied to this instance of a NxsGenericSaxOutputWrapper so
		// that all of the wrapped streams will get open and close messages
		typedef NxsGenericSaxOutputWrapper<NxsWritableStream> NxsSaxOutputWrapper;
		typedef NxsGenericSaxOutputWrapper<NxsPlotStream> NxsPlotSaxOutputWrapper;
		typedef NxsGenericSaxOutputWrapper<std::ofstream> NxsFileSaxOutputWrapper;
		std::list<NxsSaxOutputWrapper> streamWrappers;
		std::list<NxsPlotSaxOutputWrapper> plotWrappers;
		std::list<NxsFileSaxOutputWrapper> fileWrappers;
		
	};
	
typedef NxsGenericSaxOutputWrapper<NxsWritableStream> NxsSaxOutputWrapper;
	
NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & , int i);	
NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & , long l);	
NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & , unsigned l);	
NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & , unsigned long l);	
NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & , double d);
NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & , const char *c);
NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & , char c);
NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & , const std::string &s);

inline NxsSaxOutputWrapper & operator<<(NxsSaxOutputWrapper & ow, int i)
	{
	ow.outputStream << i;
	return ow;
	}

inline NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & ow, long l)
	{
	ow.outputStream << l;
	return ow;
	}

inline NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & ow, unsigned l)
	{
	ow.outputStream << l;
	return ow;
	}

inline NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & ow, unsigned long l)
	{
	ow.outputStream << l;
	return ow;
	}

inline NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & ow, double d)
	{
	ow.outputStream << d;
	return ow;
	}

inline NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & ow, const char * c)
	{
	ow.outputStream << GetEscapedMessageForOutput<NxsWritableStream>(std::string(c));
	return ow;
	}

inline NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & ow, char c)
	{
	ow.outputStream << GetEscapedMessageForOutput<NxsWritableStream>(std::string(1,c));
	return ow;
	}

inline NxsSaxOutputWrapper 	& operator<<(NxsSaxOutputWrapper & ow, const std::string & s)
	{
	ow.outputStream << GetEscapedMessageForOutput<NxsWritableStream>(s);
	return ow;
	}


#endif

