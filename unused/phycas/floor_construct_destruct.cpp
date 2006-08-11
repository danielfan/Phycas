#include "phycas/force_include.h"
#include "phycas/floor.hpp"
#include "ncl/output/nxs_input.hpp"
#include "ncl/misc/utilities.hpp"
#include "ncl/misc/nxs_file_path.hpp"
#include "ncl/output/nxs_input.hpp"
#include "phycas/rand/distribution_manager.hpp"
#include "phycas/modules/mcmc/mcmc.hpp"
#include "phycas/modules/sampler/sampler.hpp"
#include "phycas/misc/pwd.hpp"

#if defined(SCM_MODULE) //POL-23Nov2004
#	include "phycas/modules/scm/scm.hpp"
#endif

#if defined(GG_MODULE) //POL-14April2005
#	include "phycas/modules/gg/gg.hpp"
#endif

#if defined(PHYLODIVERSITY_MODULE) //POL-3May04
#	include "phycas/modules/phylodiversity/pd.hpp"
#endif

#if defined(SRQ_MODULE) //POL 9-September-2004
#	include "phycas/modules/srq/srq.hpp"
#endif

#include "phycas/taxa/taxa_manager.hpp"
#include "phycas/trees/trees_manager.hpp"
#include "phycas/characters/characters_manager.hpp"
#include "phycas/trees/split_manager.hpp" 
#include "phycas/trees/topology_manager.hpp" 
#include "phycas/trees/tree_node.hpp"
#if defined (CORBA_PHYCAS)
#   include "phycas/corba/corba_wrapper.hpp"
#endif

//	Blocks are included so that we can cast from NxsXXXBlock * returned from the GetXXXBlockReader() to NxsBlock *
//
#include "ncl/taxa/nxs_taxa_block.hpp"
#include "ncl/trees/nxs_trees_block.hpp"
#include "ncl/characters/nxs_characters_block.hpp"
#include "ncl/nxs_token.hpp"
#undef DEBUGGING_SCM_AS_SERVICE
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the `id' data member to "PHYCAS" and calls the FactoryDefaults member function to perform the 
|	remaining initializations. The data member `'' is set to NULL so that memory will be allocated for it in
|	FactoryDefaults.
*/
PhoFloor::PhoFloor()
	:
#	if defined(READ_PHYCAS_BLOCK)
		NxsCommandManagerBlock("PHYCAS", true, true, true, true),
#	endif
	NxsCharListener(NULL),
	NxsTreeListener(NULL),
	NxsTaxaListener(NULL),
	quitNow(false),
	outputMgrRef(NxsOutputManager::GetInstance()),
	notifyOfUnknownBlocks(true)
	{
#   if defined (CORBA_PHYCAS)
		cipresCORBA::CorbaWrapper::GetInstance(); //trigger corba initialization so we crash fast, if we are going to crash
#   endif
	phoTaxaMgr = boost::shared_ptr<PhoTaxaManager>(new PhoTaxaManager());
	TreeNode::gTaxaMgr = phoTaxaMgr;
	phoTreesMgr = boost::shared_ptr<PhoTreesManager>(new PhoTreesManager(*phoTaxaMgr.get()));
	phoCharactersMgr = boost::shared_ptr<PhoCharactersManager>(new PhoCharactersManager(*phoTaxaMgr.get()));

#	if defined(SCM_MODULE) //POL-23Nov2004
		scm = boost::shared_ptr<SCM>(new SCM(outputMgrRef, *phoTreesMgr.get()));
#	endif

#	if defined(GG_MODULE) //POL-14April2005
		gg = boost::shared_ptr<GG>(new GG(outputMgrRef, *phoTreesMgr.get()));
#	endif

#	if !defined (CONSOLE_PHOREST) && ! defined (NCL_OUTPUT_XML)
#		error Need to initialize output streams
#	endif
#   if defined(NCL_PRINT_COMMAND_STATE_TO_HIDDEN) && NCL_PRINT_COMMAND_STATE_TO_HIDDEN
		SetPrintStateToHiddenStream(true);
		AddCommandEnvironment(NxsCommandManager::IndexManagerInfo("char_set_manager", phoCharactersMgr.get()));
		AddCommandEnvironment(NxsCommandManager::IndexManagerInfo("tax_set_manager", phoTaxaMgr.get()));
		AddCommandEnvironment(NxsCommandManager::IndexManagerInfo("tree_set_manager", phoTreesMgr.get()));
#   endif	

#	if defined(READ_PHYCAS_BLOCK)
		//	Add all of the phycas commands (this call results in a call to PhoFloor::AddRecognizedCommands()
		NxsCommandManagerBlock::InitializeRecognizedCommands();
		NxsReader::Add(this);	// PhoFloor is the NxsReader and is a NxsBlock (for Phorest blocks);
#	endif
	//	Add the other recognized blocks
	NxsReader::Add(phoTaxaMgr->GetTaxaBlockReader());
	NxsReader::Add(phoTreesMgr->GetTreesBlockReader());
	NxsReader::Add(phoCharactersMgr->GetCharsBlockReader());
#	if defined(DEBUGGING_SCM_AS_SERVICE)
		StrVec vec;
		vec.push_back("(7:0.089964,(8:0.050264,9:0.034803):0.009528,16:0.055375)");
		vec.push_back("(1:0.232837,(((2:0.073950,4:0.049896):0.006671,3:0.084989):0.018558,(7:0.062061,((8:0.042580,9:0.042487):0.014108,16:0.050796):0.026376):0.013668):0.019356,(5:0.098338,6:0.081727):0.021046)");
		vec.push_back("(7:0.092413,(((((8:0.032800,(19:0.039050,20:0.031417):0.004826):0.002089,(12:0.044456,13:0.033808):0.000695):0.000952,(17:0.019621,18:0.018780):0.038687):0.007391,11:0.043775):0.005896,9:0.039206):0.005125,16:0.052926)");
		vec.push_back("(7:0.091431,((8:0.053804,(9:0.030934,10:0.035387):0.003327):0.002629,(14:0.047441,15:0.025483):0.011044):0.004281,16:0.053908)");
		NewickSCM(vec);
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Closes `logf' if it is open and deletes memory allocated to `next_command'.
*/
PhoFloor::~PhoFloor()
	{
	}



