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

#ifndef NCL_NXS_TAXA_MANAGER_H
#define NCL_NXS_TAXA_MANAGER_H

#include "phycas/src/ncl/taxa/base_taxa_manager.hpp"
class NxsTaxaBlock;
class NxsCommandManager;
typedef boost::weak_ptr<NxsTaxaBlock> NxsTaxaBlockAlias;
class NxsIndexSetCmdOption;
class NxsTaxaManager;
NxsIndexSetCmdOption *InstantiateTaxSetOption(const std::string &n, NxsIndexSet *manip, const std::string &d, NxsTaxaManager *tmPtr, bool  persist, CmdPermissionLevel pLevel);

class NxsTaxaManager;

class NxsTaxaManager : 
  public BaseTaxaManager, 
  public NxsBlockManager
	{
	private :
	public:
		CmdResult 				ReplaceTaxa(const VecString &v);
		NxsTaxaBlockAlias 		GetTaxaBlockReader();
		CmdResult				NewBlockRead(NxsBlock *);
		NxsTaxaManager();
		virtual ~NxsTaxaManager(){}
		
		void					WarnBeforeClearing(bool w)
			{
			warnBeforeClearingTaxa = w;
			}
	protected:
	
		void			DisplayTaxSets(bool compact);
		void			DisplayTaxSet(const NxsIndexSet &s, bool compact);
		
		
		const VecString reservedTaxAndSetNames;
		
		typedef boost::shared_ptr<NxsTaxaBlock> SharedTaxaBlockPtr;
		typedef std::vector<SharedTaxaBlockPtr> VecOfTaxaBlockPtrs;
		VecOfTaxaBlockPtrs 	taxaBlocks;
		bool				warnBeforeClearingTaxa;
	};
#endif
