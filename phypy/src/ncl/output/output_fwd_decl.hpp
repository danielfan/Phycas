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

#if ! defined (OUTPUT_FWD_DECL_HPP)
#define OUTPUT_FWD_DECL_HPP
//	if NCL_USE_NXS_CONSOLE_OUTPUT is defined, the the file nxsstdoutputstream.h will be included by nxsoutput.h (and nxsstdoutputstream.cpp) 
//	must be included in the project that builds ncl.
//	if NCL_USE_NXS_CONSOLE_OUTPUT is undefined, then the NCL_USER_OUTPUT_HEADER should be defined to the name of the header file that 
//		includes the typedefs and classes that are required by nxsoutput.h   The file will be included in the nxsoutput.h header
#if defined(NCL_USE_STD_OUTPUT)
#	include <iostream>
	class StdOutputManager;
	typedef std::ostream			NxsWritableStream;
	typedef StdOutputManager		NxsOutputManager;
#	undef NCL_HAS_TABULAR_OUTPUT
#   undef NCL_USER_OUTPUT_HEADER
#   undef NCL_OUTPUT_XML
#else
#	if defined(NCL_USE_NXS_CONSOLE_OUTPUT)
		class NxsStdOutputStream;
		class NxsStdOutputManager;
		typedef NxsStdOutputStream  NxsWritableStream;
		typedef NxsStdOutputManager		NxsOutputManager;
#		define NCL_HAS_TABULAR_OUTPUT
#		define NCL_USE_SINGLE_STREAM_OUTPUT_MGR
#		undef NCL_USER_OUTPUT_HEADER
#		undef NCL_OUTPUT_XML
#	else
#		if ! defined(NCL_USER_SUPPLIED_OUTPUT)
#			error ("define NCL_USE_STD_OUTPUT or NCL_USE_NXS_CONSOLE_OUTPUT or (NCL_USER_SUPPLIED_OUTPUT AND NCL_USER_OUTPUT_HEADER)  ")
#		endif
#		include NCL_USER_OUTPUT_FWD_DECLARATIONS
#	endif
#	if !defined(NCL_USE_SINGLE_STREAM_OUTPUT_MGR) && !defined(NCL_USER_OUTPUT_MGR_HEADER)
#		error define NCL_USER_OUTPUT_MGR_HEADER "your_header_file_name.h"
#	endif
#endif

#include "ncl/singleton.hpp"
#if defined(USE_LOKI_SINGLETON)
	typedef Loki::CreateUsingNew<NxsOutputManager> NxsOutputManagerCreator;
	typedef Loki::SingletonHolder<NxsOutputManager, Loki::CreateUsingNew, Loki::SingletonWithLongevity> NxsOutputManagerSingletonHolder;
#else
	typedef ncl::SingletonHolder<NxsOutputManager> NxsOutputManagerCreator;
	typedef ncl::SingletonHolder<NxsOutputManager> NxsOutputManagerSingletonHolder;
#endif

//for greater expresssivity, we typedef NxsWritableStream to a name that conveys the context of the stream
typedef NxsWritableStream 		NxsErrorStream;
typedef NxsWritableStream  		NxsHiddenQueryStream;
typedef NxsWritableStream  		NxsOutputStream;
typedef NxsWritableStream  		NxsOutputCommentStream;
typedef NxsWritableStream  		NxsStatusStream;
typedef NxsWritableStream  		NxsReturnValStream;
typedef NxsWritableStream  		NxsWarningStream;
typedef NxsWritableStream & (*NxsWritableStreamManipPtr)(NxsWritableStream &); // typedef the outputstream manipulator used for calls such as << endl

class NxsPlotStream; // derived from NxsWritableStream.  Not a typedef so that we can specialize templates for NxsPlotStream

typedef unsigned				NxsOutputOperationStatusID;

class NxsOutputStreamWrapper;
typedef NxsOutputStreamWrapper & (*NxsOutputStreamWrapperManipPtr)(NxsOutputStreamWrapper &); // typedef the outputstream manipulator used for calls such as << endl

namespace ncl
	{
	NxsWritableStream & flush(NxsWritableStream & om);
	NxsWritableStream & endl(NxsWritableStream & om);
	void Prompt(NxsWritableStream & outStream, const char * w);
	NxsOutputStreamWrapper & flush(NxsOutputStreamWrapper & om);
	NxsOutputStreamWrapper & endl(NxsOutputStreamWrapper & om);
	}

enum NCLQueryOrPromptType
	{
	kGenericOutput,
	kPromptForFileChooser,
	kPromptForAnyString,
	kPromptChoices,
	kPromptAlert,
	kPromptCancelOK,
	kPromptNoYes
	};
enum NCLOutputFormatHint
		{
		kExpanded			= 0x00,
		kNoTabular			= 0x00,
		kCompact			= 0x01,
		kTabular			= 0x02
		};

#endif

