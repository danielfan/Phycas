#ifndef NCL_NXS_XML_SOCKET_OUTPUT_MGR_H
#define NCL_NXS_XML_SOCKET_OUTPUT_MGR_H

#include "ncl/nxs_defs.hpp"

#include <deque>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include "ncl/output/nxs_user_query.hpp" //@POL-15feb2004 Mark, ok?
#include "ncl/output/nxs_input.hpp"

#if defined(NCL_USER_OUTPUT_HEADER)
#   include NCL_USER_OUTPUT_HEADER
#endif

class Socket;
class ServerSocket;
class NxsOutputDestinationDescription;
typedef boost::shared_ptr<ServerSocket> ServerSocketShPtr;
#   if defined(WIN_PHOREST)
#		define PHYCAS__STDCALL __stdcall
#   else
#		define PHYCAS__STDCALL
#   endif
typedef UInt (PHYCAS__STDCALL *AnswerSocketCallBack)(void *);

/*----------------------------------------------------------------------------------------------------------------------
|	Class that supplies pointers to the "fundamental" streams for all types of output and the user query.
|	In this implementation, there is only one output stream and all more specific types alias it.
|	GUI versions may require a different implementation.
*/
class NxsXMLSocketOutputManager
	{
	public:
		STATIC_SINGLETON NxsXMLSocketOutputManager & GetInstance();
		STATIC_SINGLETON const NxsXMLSocketOutputManager  & GetConstInstance()
			{
			return  GetInstance();
			}
		NxsErrorStream 		 	  * GetErrorStreamPtr();
		NxsHiddenQueryStream 	  * GetHiddenQueryStreamPtr() 
			{
			return & hiddenQueryStream;
			}
		NxsIO::NxsInput * GetInputStreamPtr()
			{
			return &inputStream;
			}
		NxsOutputStream 		  * GetOutputStreamPtr();
		NxsOutputCommentStream	  * GetOutputCommentStreamPtr();
		NxsPlotStream			  * GetPlotStreamPtr()
			{
			return & plotStream;
			}
		NxsReturnValStream 		  * GetReturnValueStreamPtr()
			{
			return returnValStream;
			}
		NxsStatusStream 		  * GetStatusStreamPtr();
		NxsUserQuery			  * GetUserQueryPtr();
		NxsWarningStream		  * GetWarningStreamPtr();
		NxsOutputStreamWrapperShPtr	GetGenericOutputStreamShPtr(const NxsOutputDestinationDescription &);
		void Prompt()
			{
			outputStream.Prompt();
			}

		unsigned 				GetOutputWidth() const;
		void					SetOutputWidth(unsigned outputWidth);
		void 					SetSpawnedSocket(Socket *s)
		 	{
			spawnedSocket = s;
		 	} 
		Socket				 * GetSpawnedSocket()
			{
			return spawnedSocket;
			}
			
			///start of a process that might merit a new progress bar 
		NxsOutputOperationStatusID  StartStatusDisplay(const std::string &m, const bool updatesWillBePosted);	
			/// update progress bar with specified ID, 
		void						UpdateStatusDisplay(const float proportionDone, const NxsOutputOperationStatusID, const std::string optionalMsg = std::string());
			/// process with specified ID is finished
		void						EndStatusDisplay(const NxsOutputOperationStatusID, const std::string m =  std::string()); 
			/// One-time alert of an event (no progress bar should be created in a GUI app)
		void						AlertStatusDisplay(const std::string &msg); 
		
		void						AcceptConnection(AnswerSocketCallBack);
		
		bool						SetPortNumber(UInt newPortN);
		bool IsOpen() const
			{
			return (spawnedSocket != NULL);
			}
		void SendLine(const std::string & r) const;
	private:
		NxsXMLSocketOutputManager();
		~NxsXMLSocketOutputManager();
		UInt						portN;
		ServerSocketShPtr			serverSocket; //master socket used to accept new connections
		Socket					  * spawnedSocket;
		
	protected:	
		NxsErrorStream 				errorStream;			
		NxsHiddenQueryStream		hiddenQueryStream;
		NxsOutputStream 			outputStream;			// keep outputStream before  aliases 
		NxsOutputCommentStream		outputCommentStream; 	
		NxsPlotStream				plotStream; 	
		NxsWarningStream  			warningStream;			
		NxsOutputStream				queryOut; 
		NxsStatusStream			  * statusStream;			/// temp alias to output
		NxsReturnValStream		  * returnValStream;		//alias to output
		NxsIO::NxsInput				inputStream;
		NxsUserQuery    			userQuery;				/* user IO with an alias to outputStream */
		friend struct NxsXMLSocketOutputManagerCreator;
	};

/*----------------------------------------------------------------------------------------------------------------------
|	returns a pointer to the stream that is used to display output comments found in the file
*/
inline NxsOutputCommentStream * NxsXMLSocketOutputManager::GetOutputCommentStreamPtr()
	{
	return & outputCommentStream;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns a pointer to the stream that is used to display generic output to the user
*/
inline NxsOutputStream * NxsXMLSocketOutputManager::GetOutputStreamPtr()
	{
	return & outputStream;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns a pointer to the stream that is used to display warnings to the user
*/
inline NxsWarningStream * NxsXMLSocketOutputManager::GetWarningStreamPtr()
	{
	return & warningStream;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns a pointer to the stream that is used to display errors to the user
*/
inline NxsErrorStream * NxsXMLSocketOutputManager::GetErrorStreamPtr()
	{
	return & errorStream;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	returns a pointer to the NxsUserQuery object used to ask users "multiple-choice"-style questions
*/
inline NxsUserQuery * NxsXMLSocketOutputManager::GetUserQueryPtr()
	{
	return &userQuery;
	}
		
/*----------------------------------------------------------------------------------------------------------------------
|	Changes the output streams width (# of characters before a carriage return will be added).
*/
inline void NxsXMLSocketOutputManager::SetOutputWidth(unsigned outputWidth)
	{
	outputStream.SetWidth(outputWidth);
	//userQuery.SetWidth(outputWidth);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns the output streams width (# of characters before a carriage return will be added).
*/
inline unsigned NxsXMLSocketOutputManager::GetOutputWidth() const
	{
	return outputStream.GetWidth();
	}


#endif //ifndef NCL_NXS_XML_SOCKET_OUTPUT_MGR_H
