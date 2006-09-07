#ifndef NCL_NXSSTDOUTPUTSTREAM_H
#define NCL_NXSSTDOUTPUTSTREAM_H

#include "ncl/nxs_defs.hpp"
#include "ncl/output/nxs_typist.hpp"
#include "ncl/output/nxs_table_cell.hpp"
#include "ncl/output/nxs_table.hpp"
#include "ncl/misc/string_extensions.hpp"
#include <iostream>
#include <boost/noncopyable.hpp>

#if defined (SUPPORT_TREE_OUTPUT)
	class Tree;
#endif
class NxsIndexInfo;

/*----------------------------------------------------------------------------------------------------------------------
|	Provides "put-to" operator interface for output to a console and/or a file.
|	calling NxsStdOutputStream << endl causes the output stream to be flused.   
|	Formatting and wrapping are provided by the NxsTypist and NxsTable classes.  
|	The NxsTypist is hidden from the users of NxsStdOutputStream, but tabular output must go through GetTablePtr() and
|	PrintTable() functions.
*/
class NxsStdOutputStream NON_COPYABLE
	{
	public :
		STATIC_DATA const unsigned kDefaultPrintWidth = 80;
		
		NxsStdOutputStream();
		
		unsigned	GetWidth() const;
		unsigned	GetSuccessiveNewlines() const;
		
		void 		DirectOutputToScreen(bool displayOutput = true);
		void 		DirectOutputToFile(std::ostream *newLogFile);
		void		DisplaySet(const NxsIndexSet &s, bool useNumbers, const NxsIndexInfo &);
		NxsTable  * GetTablePtr();
		void 		PrintMessage();
		void 		PrintTable(NxsTable *t);
		void 		SetNDigitsAfterDecimal(unsigned i);
		void		SetFont(const std::string &f) {fontName = f;}	// only used in SIOUX_PHOREST NxsStdOutputStream &operator<<(const Tree &);
		void		SetSuccessiveNewlines(unsigned n);
		void		SetWidth(unsigned outputWidth);
		void 		StopLoggingToFile();
	
		//	message and dblFormat are public to avoid having to make all of the output operators friends.  These fields should not be manipulated
		//	 externally (unless you'r sure you know what you are doing)
		std::string		message; 			/* current stream (output that has not been printed) */
		DblFormatter	dblFormat;			/* class used to format doubles in a consistent manner */
	protected :
	
		void PrintLine(const std::string &line) const;
		
		bool		    outputToMonitor; 	/* if true, output is flushed to cout */
		std::ostream   *logFile;			/* pointer to log file (NULL if not logging to a file) */
		std::string		fontName;			/* name of the current font for outputting to the console */ 
	private:
		NxsTypist		typist;		/* used to wrap text and avoid excessive white space */
		NxsTable		table;		/* for use in formatting tabular output */
		
		friend class NxsTypist;
	};
// The NxsOutputStream Interface
//
NxsStdOutputStream 	& operator<<(NxsStdOutputStream & , int i);	
NxsStdOutputStream 	& operator<<(NxsStdOutputStream & , long l);	
NxsStdOutputStream 	& operator<<(NxsStdOutputStream & , unsigned l);	
NxsStdOutputStream 	& operator<<(NxsStdOutputStream & , unsigned long l);	
NxsStdOutputStream 	& operator<<(NxsStdOutputStream & , double d);
NxsStdOutputStream 	& operator<<(NxsStdOutputStream & , const char *c);
NxsStdOutputStream 	& operator<<(NxsStdOutputStream & , char c);
NxsStdOutputStream 	& operator<<(NxsStdOutputStream & , const std::string &s);
NxsStdOutputStream & operator<<(NxsStdOutputStream & , NxsWritableStreamManipPtr);
#if defined (SUPPORT_TREE_OUTPUT)
	NxsStdOutputStream & operator<<(NxsStdOutputStream & , const Tree &);
#endif

inline NxsTable  *NxsStdOutputStream::GetTablePtr()
	{
	table.Reset();
	return &table;
	}

inline void	NxsStdOutputStream::SetWidth(
  unsigned outputWidth)
	{
	typist.SetWidth(outputWidth);
	table.SetPrnWidth(outputWidth);
	}
	
inline unsigned	NxsStdOutputStream::GetWidth() const
	{
	return typist.GetWidth();
	}
	
inline void NxsStdOutputStream::PrintMessage()
  	{
	if (!message.empty()) //@pol typist should deal with message of size 0 itself
		{
  		typist.Run(message);
  		std::cout << std::flush;
		message.clear();
		}
	}

namespace ncl
{
inline void Prompt(NxsStdOutputStream & outStream, const char * w)
	{
	outStream.PrintMessage();
	std::cout << w << std::flush;
	}
} // namespace ncl		

inline void NxsStdOutputStream::DirectOutputToScreen(bool displayOutput)
	{
	outputToMonitor = displayOutput;
	}
		
inline void NxsStdOutputStream::DirectOutputToFile(std::ostream *newLogFile)
	{
	logFile = newLogFile;
	}

inline void NxsStdOutputStream::StopLoggingToFile()
	{
	DirectOutputToFile(NULL);
	}
	
inline NxsStdOutputStream & operator<<(NxsStdOutputStream & outStream, NxsWritableStreamManipPtr funcPtr)
	{
	return (*funcPtr)(outStream);
	}

/**
 * @method PrintLine [void:public]
 * @param line [std::string] line of text to send to output
 *
 * Called by Typist to output individual lines chomped from message.
 * This method is virtual to allow it to be overridden in derived
 * classes to handle output on different platforms.
 */
inline void NxsStdOutputStream::PrintLine(const std::string &line) const
	{
	if (outputToMonitor)
		std::cout << line << '\n';
	if (logFile != NULL)
		*logFile << line << '\n';
	}
	
inline  NxsStdOutputStream & operator<<(NxsStdOutputStream & outStr, int i) 
	{
	outStr.message << i; 
	return outStr;
	}
	
inline  NxsStdOutputStream & operator<<(NxsStdOutputStream & outStr, long l) 
	{
	outStr.message << l; 
	return outStr;
	}
	
inline  NxsStdOutputStream & operator<<(NxsStdOutputStream & outStr, unsigned i) 
	{
	outStr.message << i; 
	return outStr;
	}
	
inline  NxsStdOutputStream & operator<<(NxsStdOutputStream & outStr, unsigned long l) 
	{
	outStr.message << l; 
	return outStr;
	}
	
inline  NxsStdOutputStream & operator<<(NxsStdOutputStream & outStr, double d) 
	{
	outStr.dblFormat.FormatDouble(outStr.message,  d);
	return outStr;
	}
	
inline  NxsStdOutputStream & operator<<(NxsStdOutputStream & outStr, const char *c) 
	{
	outStr.message << c; 
	return outStr;
	}
	
inline  NxsStdOutputStream & operator<<(NxsStdOutputStream & outStr, char c) 
	{
	outStr.message.push_back(c);
	return outStr;
	}
	
inline  NxsStdOutputStream & operator<<(NxsStdOutputStream & outStr, const std::string &s) 
	{
	outStr.message.append(s);
	return outStr;
	}

namespace ncl
{
inline  NxsStdOutputStream &flush(NxsStdOutputStream& om)
	{
	om.PrintMessage();
	return om;
	}
	

inline  NxsStdOutputStream &endl(NxsStdOutputStream& om)
	{
	om << '\n';
	return flush(om);
	}
} // namespace ncl
inline  void NxsStdOutputStream::SetNDigitsAfterDecimal(unsigned i)
	{
	dblFormat.digitsAfterDecimal = i;
	}

inline void NxsStdOutputStream::PrintTable(NxsTable *p)
	{
	NXS_ASSERT(p == &table);	//you should be getting the Table pointer from NxsStdOutputStream::GetTablePtr so that it has the correct width etc
	try 
		{
		p->Show(*this);
		p->Reset();
		}
	catch(NxsTable::NxsX_InsufficientWidth &)
		{
		*this << "output width is too small\n" << ncl::flush;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `succ_newl', the number of successive newlines at the end of the the most recent output.
*/
inline unsigned	NxsStdOutputStream::GetSuccessiveNewlines() const
	{
	return typist.GetSuccessiveNewlines();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of `succ_newl', the number of successive newlines at the end of the the most recent output.
*/
inline void NxsStdOutputStream::SetSuccessiveNewlines(unsigned s)
	{
	typist.SetSuccessiveNewlines(s);
	}


	
#endif
