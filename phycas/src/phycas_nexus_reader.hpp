/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
|																			  |
|  This program is free software; you can redistribute it and/or modify		  |
|  it under the terms of the GNU General Public License as published by		  |
|  the Free Software Foundation; either version 2 of the License, or		  |
|  (at your option) any later version.										  |
|																			  |
|  This program is distributed in the hope that it will be useful,			  |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of			  |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the			  |
|  GNU General Public License for more details.								  |
|																			  |
|  You should have received a copy of the GNU General Public License along	  |
|  with this program; if not, write to the Free Software Foundation, Inc.,	  |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.				  |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <boost/shared_ptr.hpp>
#include "ncl/nxspublicblocks.h"
#include "phycas/src/nxs_file_path.hpp"
#include "phycas/src/char_super_matrix.hpp"

enum CmdResult
	{
	kCmdSucceeded,
	kCmdFailedSilent,
	kCmdFailedGenerateMessage
	};

class PhycasNexusReader: public PublicNexusReader
	{
	public: 
		PhycasNexusReader(const int blocksToRead, NxsReader::WarningHandlingMode mode=NxsReader::WARNINGS_TO_STDERR) 
			:PublicNexusReader(blocksToRead, mode),
			activeTaxaBlockIndex(UINT_MAX),
			activeCharactersBlockIndex(UINT_MAX)
			{
			this->SetWarningOutputLevel(SKIPPING_CONTENT_WARNING);
			}
			
		~PhycasNexusReader();
			
		void ReadFilePath(const std::string & filePath)
			{
			NxsInFilePath fp(filePath);
			ReadNxsFilePath(fp);
			}
			
		std::vector<std::string> GetTaxLabels() const;
		
		const std::vector<NxsFullTreeDescription> & GetTrees() const;
		
		unsigned GetNChar() const;
		
		std::string GetErrorMessage()
			{
			return errorMsg;
			}
			
		virtual void Clear();
			
		virtual void ClearUsedBlockList()
			{
			PublicNexusReader::ClearUsedBlockList();
			activeTaxaBlockIndex = UINT_MAX;
			activeCharactersBlockIndex = UINT_MAX;
			}
			
		NxsCharactersBlock * getActiveCharactersBlock() const;
		
		virtual void NexusError(const std::string &msg, file_pos pos, unsigned line, unsigned col, CmdResult , NxsBlock* );
		
	private:
	
		void								ReadNxsFilePath(NxsInFilePath & filePath);
		NxsTaxaBlock *						getActiveTaxaBlock() const;

		std::string							errorMsg;
		unsigned 							activeTaxaBlockIndex;
		std::vector<NxsFullTreeDescription>	trees;
		unsigned 							activeCharactersBlockIndex;
		file_pos							filePositionOfError;
		unsigned							lineOfError ;
		unsigned							columnOfError ;

	};


typedef boost::shared_ptr<PhycasNexusReader> PhycasNexusReaderShPtr;

