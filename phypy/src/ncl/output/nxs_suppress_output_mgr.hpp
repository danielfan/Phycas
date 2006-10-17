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

#ifndef NCL_NXS_SUPPRESS_OUTPUT_MGR_H
#define NCL_NXS_SUPPRESS_OUTPUT_MGR_H

#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/output/nxs_output_stream.hpp"

class NxsOutputDestinationDescription;
/*----------------------------------------------------------------------------------------------------------------------
|	Class that supplies pointers to the "fundamental" streams for all types of output and the user query.
|	In this implementation, there is only one output stream and all more specific types alias it.
|	GUI versions may require a different implementation.
*/
class NxsSuppressOutputManager
	{
	public:
		STATIC_SINGLETON NxsSuppressOutputManager & GetInstance() 
			{
			STATIC_SINGLETON NxsSuppressOutputManager theOutputSuppressor;
			return theOutputSuppressor;
			}
		NxsErrorStream 		 	  * GetErrorStreamPtr() { return NULL;}
		NxsHiddenQueryStream 	  * GetHiddenQueryStreamPtr()  { return NULL;}
		NxsIO::NxsInput			  * GetInputStreamPtr() { return NULL;}
		NxsOutputStream 		  * GetOutputStreamPtr() { return NULL;}
		NxsOutputCommentStream	  * GetOutputCommentStreamPtr() { return NULL;}
		NxsPlotStream			  * GetPlotStreamPtr() { return NULL;}
		NxsReturnValStream 		  * GetReturnValueStreamPtr() { return NULL;}
		NxsStatusStream 		  * GetStatusStreamPtr() { return NULL;}
		NxsUserQuery			  * GetUserQueryPtr() { return NULL;}
		NxsWarningStream		  * GetWarningStreamPtr() { return NULL;}
		NxsOutputStreamWrapperShPtr	GetGenericOutputStreamShPtr(const NxsOutputDestinationDescription &);
		void						Prompt() {}		

		unsigned 				GetOutputWidth() const
			{
			return 80U; //bogus value here that is unlikely to break code expecting a reasonable number
			}

		void					SetOutputWidth(unsigned ) {}
			///start of a process that might merit a new progress bar 
		NxsOutputOperationStatusID  StartStatusDisplay(const std::string &, const bool )
			{
			STATIC_DATA NxsOutputOperationStatusID dummy(0);
			return dummy++;
			}
			/// update progress bar with specified ID, 
		void UpdateStatusDisplay(const float proportionDone, const NxsOutputOperationStatusID, const std::string optionalMsg = std::string())
			{
			}
			/// process with specified ID is finished
		void EndStatusDisplay(const NxsOutputOperationStatusID, const std::string m =  std::string()) 
			{
			}
			/// One-time alert of an event (no progress bar should be created in a GUI app)
		void AlertStatusDisplay(const std::string &msg)
			{
			}
			
		NxsSuppressOutputManager() {} //@POL 27-Oct-2005 Mark, hope it was ok to provide empty definitions for constructor and destructor
		~NxsSuppressOutputManager() {}
	};

#endif //ifndef NCL_NXS_XML_SOCKET_OUTPUT_MGR_H
