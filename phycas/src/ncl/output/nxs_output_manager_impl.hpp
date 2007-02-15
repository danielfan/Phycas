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

#if ! defined(NCL_OUTPUT_MGR_H) 
#define NCL_OUTPUT_MGR_H
#include <boost/shared_ptr.hpp>
#include "phycas/src/ncl/nxs_defs.hpp"
#include "phycas/src/ncl/output/nxs_input.hpp"
class NxsOutputDestinationDescription;
class NxsOutputStreamWrapper;
typedef boost::shared_ptr<NxsOutputStreamWrapper> NxsOutputStreamWrapperShPtr;
#if defined(NCL_USE_STD_OUTPUT)
#	include "phycas/src/ncl/output/nxs_user_query.hpp" //@@TESTING OUTPUT

	class StdOutputManager
		{
#			if defined(NXS_IO_SUPPORT_INTERACTIVE_MODE)
				NxsInput		inputStream;
#			endif
			NxsUserQuery 	query;
		public:
			STATIC_SINGLETON StdOutputManager & GetInstance();
			STATIC_SINGLETON const StdOutputManager  & GetConstInstance()
				{
				return  GetInstance();
				}
			NxsOutputStreamWrapperShPtr		GetGenericOutputStreamShPtr(const NxsOutputDestinationDescription &);
			NxsErrorStream 		   *GetErrorStreamPtr()
				{
#				if defined(USE_STD_OUT_FOR_STD_ERR) 
					return &std::cout;
#				else
					return &std::cerr;
#				endif
				}
			NxsHiddenQueryStream  * GetHiddenQueryStreamPtr()
				{
				return & std::cout;
				}
			NxsIO::NxsInput * GetInputStreamPtr()
				{
#				if defined(NXS_IO_SUPPORT_INTERACTIVE_MODE)
					return &inputStream;
#				else
					return NULL;
#				endif
				}
			NxsOutputCommentStream *GetOutputCommentStreamPtr()
				{
				return &std::cout;
				}
			NxsOutputStream 	   *GetOutputStreamPtr()
				{
				return &std::cout;
				}
			NxsPlotStream		  * GetPlotStreamPtr() const
				{
				return NULL;
				}
			NxsUserQuery		   *GetUserQueryPtr()
				{
				return &query;
				}
			NxsWarningStream 	   *GetWarningStreamPtr() 
				{
				return &std::cout;
				}
			NxsStatusStream 	   *GetStatusStreamPtr() 
				{
				return &std::cout;
				}
			NxsReturnValStream 	   *GetReturnValueStreamPtr() 
				{
				return &std::cout;
				}
			NxsOutputOperationStatusID  StartStatusDisplay(const std::string &m, const bool updatesWillBePosted)
				{
				std::cout << m << std::endl;
				return NxsOutputOperationStatusID(1); //@@
				}
	
			void	UpdateStatusDisplay( const float , const NxsOutputOperationStatusID, const std::string)
				{
				}
			
			void EndStatusDisplay(const NxsOutputOperationStatusID, const std::string msg)
				{
				std::cout << msg << std::endl;
				}
				
			void AlertStatusDisplay(const std::string &msg)
				{
				std::cout << msg << std::endl;
				}
			// NOTE The functions:
			//		unsigned 				GetOutputWidth() const;
			// 		void					SetOutputWidth(unsigned outputWidth);
			//	aren't supported for NCL_USE_STD_OUTPUT
						
		private:
			StdOutputManager()
				:
#				if defined(NXS_IO_SUPPORT_INTERACTIVE_MODE)
					inputStream(),
					userQuery(&std::cout, &inputStream)
#				else
					userQuery(&std::cout, NULL)
#				endif
				{}
			friend struct NxsOutputManagerCreator;
		};
		
			
#else	//defined (NCL_USE_STD_OUTPUT)
#	if defined(NCL_USE_SINGLE_STREAM_OUTPUT_MGR)
#		include "phycas/src/ncl/output/nxs_output.hpp" 
#		include "phycas/src/ncl/output/nxs_user_query.hpp" //@@TESTING OUTPUT
		/*----------------------------------------------------------------------------------------------------------------------
		|	Class that supplies pointers to the "fundamental" streams for all types of output and the user query.
		|	In this implementation, there is only one output stream and all more specific types alias it.
		|	GUI versions may require a different implementation.
		*/
		class NxsStdOutputManager
			{
			public:
				STATIC_SINGLETON NxsStdOutputManager & GetInstance();
				STATIC_SINGLETON const NxsStdOutputManager  & GetConstInstance()
					{
					return  GetInstance();
					}
				NxsErrorStream 		  * GetErrorStreamPtr() const;
				NxsHiddenQueryStream  * GetHiddenQueryStreamPtr() const;
				NxsIO::NxsInput * GetInputStreamPtr()
					{
#					if defined(NXS_IO_SUPPORT_INTERACTIVE_MODE)
						return &inputStream;
#					else
						return NULL;
#					endif
					}
				NxsOutputCommentStream* GetOutputCommentStreamPtr() const;
				NxsOutputStream 	  * GetOutputStreamPtr() const;
				NxsPlotStream		  * GetPlotStreamPtr() const
					{
					return NULL;
					}
				NxsReturnValStream 	  * GetReturnValueStreamPtr() const;
				NxsUserQuery		  * GetUserQueryPtr() const;
				NxsWarningStream 	  * GetWarningStreamPtr() const;
				NxsOutputStreamWrapperShPtr		GetGenericOutputStreamShPtr(const NxsOutputDestinationDescription &);

				unsigned 				GetOutputWidth() const;
				void					SetOutputWidth(unsigned outputWidth);
			
				NxsOutputOperationStatusID  StartStatusDisplay(const std::string &m, const bool )//updatesWillBePosted)
					{
					*statusStream << m << ncl::endl;
					return NxsOutputOperationStatusID(1); //@@
					}
		
				void	UpdateStatusDisplay( const float , const NxsOutputOperationStatusID, const std::string)
					{
					}
				
				void EndStatusDisplay(const NxsOutputOperationStatusID, const std::string msg)
					{
					*statusStream << msg << ncl::endl;
					}
					
				void AlertStatusDisplay(const std::string &msg)
					{
					*statusStream << msg << ncl::endl;
					}
			private:
				NxsStdOutputManager();
			protected:	
				mutable NxsOutputStream   outputStream;			/* keep outputStream first in the declartion */
				NxsOutputCommentStream  * outputCommentStream; 	/* alias to outputStream */
				NxsHiddenQueryStream	* hiddenQuery;			/* alias to outputStream */
				NxsWarningStream  		* warningStream;		/* alias to outputStream */
				NxsErrorStream 			* errorStream;			/* alias to outputStream*/
				NxsStatusStream 		* statusStream;			/* alias to outputStream*/
				NxsReturnValStream 	    * returnValStream;		/* alias to outputStream*/
#				if defined(NXS_IO_SUPPORT_INTERACTIVE_MODE)
					NxsIO::NxsInput			  inputStream;
#				endif
				mutable NxsUserQuery      userQuery;			/* user IO with an alias to outputStream */
#				if defined(USE_LOKI_SINGLETON)
					friend struct Loki::CreateUsingNew<NxsOutputManager>;
#				else
					friend struct ncl::SingletonHolder<NxsOutputManager>;
#				endif
			};

		/*----------------------------------------------------------------------------------------------------------------------
		|	Creates an output stream and makes aliases to it for the outpuComment, warning, error streams and user query.
		*/
		inline NxsStdOutputManager::NxsStdOutputManager()
			:outputStream(),
			outputCommentStream(NULL),
			warningStream(NULL),
			errorStream(NULL),
			statusStream(NULL),
#			if defined(NXS_IO_SUPPORT_INTERACTIVE_MODE)
				inputStream(),
				userQuery(&outputStream, &inputStream)
#			else
				userQuery(&outputStream, NULL)
#			endif
			{
			outputCommentStream = &outputStream; 
			warningStream = &outputStream;
			errorStream = &outputStream;
			hiddenQuery = &outputStream;
			statusStream = &outputStream;
			returnValStream = &outputStream;
			}
			
		/*----------------------------------------------------------------------------------------------------------------------
		|	returns a pointer to the stream that is used to display output comments found in the file
		*/
		inline NxsOutputCommentStream *NxsStdOutputManager::GetOutputCommentStreamPtr() const
			{
			return outputCommentStream;
			}

		/*----------------------------------------------------------------------------------------------------------------------
		|	returns a pointer to the stream that is used to display generic output to the user
		*/
		inline NxsOutputStream *NxsStdOutputManager::GetOutputStreamPtr() const
			{
			return &outputStream;
			}

		/*----------------------------------------------------------------------------------------------------------------------
		|	returns a pointer to the stream that is used to display warnings to the user
		*/
		inline NxsWarningStream *NxsStdOutputManager::GetWarningStreamPtr() const 
			{
			return warningStream;
			}

		/*----------------------------------------------------------------------------------------------------------------------
		|	returns a pointer to the stream that is used to display errors to the user
		*/
		inline NxsHiddenQueryStream *NxsStdOutputManager::GetErrorStreamPtr() const
			{
			return errorStream;
			}
		
		/*----------------------------------------------------------------------------------------------------------------------
		|	returns a pointer to the stream that is used to information to another (driving) application
		*/
		inline NxsErrorStream *NxsStdOutputManager::GetHiddenQueryStreamPtr() const
			{
			return hiddenQuery;
			}
		inline NxsReturnValStream * NxsStdOutputManager::GetReturnValueStreamPtr() const
			{
			return returnValStream;
			}
				
		/*----------------------------------------------------------------------------------------------------------------------
		|	returns a pointer to the NxsUserQuery object used to ask users "multiple-choice"-style questions
		*/
		inline NxsUserQuery *NxsStdOutputManager::GetUserQueryPtr() const
			{
			return &userQuery;
			}
				
		/*----------------------------------------------------------------------------------------------------------------------
		|	Changes the output streams width (# of characters before a carriage return will be added).
		*/
		inline void NxsStdOutputManager::SetOutputWidth(unsigned outputWidth)
			{
			outputStream.SetWidth(outputWidth);
			}

		/*----------------------------------------------------------------------------------------------------------------------
		|	returns the output streams width (# of characters before a carriage return will be added).
		*/
		inline unsigned NxsStdOutputManager::GetOutputWidth() const
			{
			return outputStream.GetWidth();
			}
#	elif defined(NCL_USE_SUPPRESS_OUTPUT_MGR)
#	else
#		include NCL_USER_OUTPUT_MGR_HEADER
#	endif //defined(NCL_USE_SINGLE_STREAM_OUTPUT_MGR)
#endif //defined (NCL_USE_STD_OUTPUT)

namespace NxsOutput
{

	// utility class verifies that the stream pointer is not null before passing the call to the output operator on
template <typename OUT_STREAM_TYPE>
class NullCheckingStreamWrapper
	{
	public :
		NullCheckingStreamWrapper(OUT_STREAM_TYPE * streamPtr)
			:stream(streamPtr)
			{}
	
		OUT_STREAM_TYPE * const stream;
	};
typedef NullCheckingStreamWrapper<NxsWritableStream>  NxsNullCheckingStream;
	// abbreviations for NxsOutputManager::GetInstance().Get_XXXXX_StreamPtr()
NxsErrorStream * GetErrorStreamPtr();
NxsHiddenQueryStream * GetHiddenQueryStreamPtr();
NxsOutputStream * GetOutputStreamPtr();
NxsOutputCommentStream * GetOutputCommentStreamPtr();
NxsPlotStream * GetPlotStreamPtr();
NxsWarningStream * GetWarningStreamPtr();
NxsReturnValStream * GetReturnValueStreamPtr();
NxsUserQuery * GetUserQueryPtr();
NxsOutputStreamWrapperShPtr GetGenericOutputStreamShPtr(const NxsOutputDestinationDescription &);
enum 	{
		kError,
		kHiddenQuery,
		kOutput,
		kOutputComment,
		kPlot,
		kReturnValue,
		kWarning
		};
template<int STREAM_ID>
NxsNullCheckingStream GetNullCheckingStream();
const NxsNullCheckingStream & operator<<(const NxsNullCheckingStream & out, const char * rOP);
const NxsNullCheckingStream & operator<<(const NxsNullCheckingStream & out, const std::string & rOP);
const NxsNullCheckingStream & operator<<(const NxsNullCheckingStream & out, NxsWritableStreamManipPtr);

#if defined(NCL_SUPPRESS_OUTPUT)
	inline NxsOutputStream * GetOutputStreamPtr()
		{
		return NULL;
		}

	inline NxsPlotStream * GetPlotStreamPtr()
		{
		return NULL;
		}

	inline NxsErrorStream * GetErrorStreamPtr()
		{
		return NULL;
		}

	inline NxsHiddenQueryStream * GetHiddenQueryStreamPtr()
		{
		return NULL;
		}
	
	inline NxsReturnValStream * GetReturnValueStreamPtr()
		{
		return NULL;
		}

	inline NxsOutputCommentStream * GetOutputCommentStreamPtr()
		{
		return NULL;
		}

	inline NxsWarningStream * GetWarningStreamPtr()
		{
		return NULL;
		}

	inline NxsUserQuery * GetUserQueryPtr()
		{
		return NULL;
		}

	inline NxsOutputStreamWrapperShPtr GetGenericOutputStreamShPtr(const NxsOutputDestinationDescription & nodd)
		{
		NxsOutputManager om;
		return om.GetGenericOutputStreamShPtr(nodd);
		}


	inline const NxsNullCheckingStream & operator<<(const NxsNullCheckingStream & out, const char * rOP)
		{
		if (out.stream != NULL)
			*out.stream << rOP;
		return out;
		}
		
	inline const NxsNullCheckingStream & operator<<(const NxsNullCheckingStream & out, const std::string & rOP)
		{
		if (out.stream != NULL)
			*out.stream << rOP;
		return out;
		}
		
	inline const NxsNullCheckingStream & operator<<(const NxsNullCheckingStream & out, NxsWritableStreamManipPtr rOP)
		{
		if (out.stream != NULL)
			*out.stream << rOP;
		return out;
		}
		
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kError>()
		{
		return NxsNullCheckingStream(NULL);
		}
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kHiddenQuery>()
		{
		return NxsNullCheckingStream(NULL);
		}
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kOutput>()
		{
		return NxsNullCheckingStream(NULL);
		}
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kOutputComment>()
		{
		return NxsNullCheckingStream(NULL);
		}
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kPlot>()
		{
		return NxsNullCheckingStream(NULL);
		}
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kReturnValue>()
		{
		return NxsNullCheckingStream(NULL);
		}
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kWarning>()
		{
		return NxsNullCheckingStream(NULL);
		}
#else // defined(NCL_USE_SUPPRESS_OUTPUT_MGR)
	inline NxsOutputStream * GetOutputStreamPtr()
		{
		return NxsOutputManager::GetInstance().GetOutputStreamPtr();
		}

	inline NxsPlotStream * GetPlotStreamPtr()
		{
		return NxsOutputManager::GetInstance().GetPlotStreamPtr();
		}

	inline NxsErrorStream * GetErrorStreamPtr()
		{
		return NxsOutputManager::GetInstance().GetErrorStreamPtr();
		}

	inline NxsHiddenQueryStream * GetHiddenQueryStreamPtr()
		{
		return NxsOutputManager::GetInstance().GetHiddenQueryStreamPtr();
		}
	inline NxsReturnValStream * GetReturnValueStreamPtr()
		{
		return NxsOutputManager::GetInstance().GetReturnValueStreamPtr();
		}

	inline NxsOutputCommentStream * GetOutputCommentStreamPtr()
		{
		return NxsOutputManager::GetInstance().GetOutputCommentStreamPtr();
		}

	inline NxsWarningStream * GetWarningStreamPtr()
		{
		return NxsOutputManager::GetInstance().GetWarningStreamPtr();
		}

	inline NxsUserQuery * GetUserQueryPtr()
		{
		return NxsOutputManager::GetInstance().GetUserQueryPtr();
		}

	inline NxsOutputStreamWrapperShPtr GetGenericOutputStreamShPtr(const NxsOutputDestinationDescription & nodd)
		{
		return  NxsOutputManager::GetInstance().GetGenericOutputStreamShPtr(nodd);
		}


	inline const NxsNullCheckingStream & operator<<(const NxsNullCheckingStream & out, const char * rOP)
		{
		if (out.stream != NULL)
			*out.stream << rOP;
		return out;
		}
		
	inline const NxsNullCheckingStream & operator<<(const NxsNullCheckingStream & out, const std::string & rOP)
		{
		if (out.stream != NULL)
			*out.stream << rOP;
		return out;
		}
		
	inline const NxsNullCheckingStream & operator<<(const NxsNullCheckingStream & out, NxsWritableStreamManipPtr rOP)
		{
		if (out.stream != NULL)
			*out.stream << rOP;
		return out;
		}
		
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kError>()
		{
		return NxsNullCheckingStream(GetErrorStreamPtr());
		}
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kHiddenQuery>()
		{
		return NxsNullCheckingStream(GetHiddenQueryStreamPtr());
		}
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kOutput>()
		{
		return NxsNullCheckingStream(GetOutputStreamPtr());
		}
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kOutputComment>()
		{
		return NxsNullCheckingStream(GetOutputCommentStreamPtr());
		}
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kPlot>()
		{
		return NxsNullCheckingStream(GetPlotStreamPtr());
		}
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kReturnValue>()
		{
		return NxsNullCheckingStream(GetReturnValueStreamPtr());
		}
	template<>
	inline NxsNullCheckingStream GetNullCheckingStream<kWarning>()
		{
		return NxsNullCheckingStream(GetWarningStreamPtr());
		}
#endif

} // namespace NxsOuput
#if ! defined(NCL_SUPPRESS_OUTPUT)
	inline NxsOutputManager & NxsOutputManager::GetInstance()
		{
		return NxsOutputManagerSingletonHolder::Instance();
		}
#endif

	

#endif
