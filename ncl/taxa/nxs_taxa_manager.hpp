#ifndef NCL_NXS_TAXA_MANAGER_H
#define NCL_NXS_TAXA_MANAGER_H

#include "ncl/taxa/base_taxa_manager.hpp"
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
