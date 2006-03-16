#include "phycas/force_include.h"
#include "phycas/modules/scm/scm.hpp"
#include "ncl/trees/nxs_trees_block.hpp"
#include "phycas/trees/trees_manager.hpp"
#include "phycas/trees/tree.hpp"
#include "phycas/trees/tree_node.hpp"
#include "phycas/taxa/taxa_manager.hpp"
#include "ncl/nxs_token.hpp"
using std::string;
using std::vector;

std::string NewickSCM(const StrVec & toMerge)
	{
	boost::shared_ptr<PhoTaxaManager> prevGlobTaxMgr = boost::shared_ptr<PhoTaxaManager>(new PhoTaxaManager()); //@ temp global should go away
	swap(prevGlobTaxMgr, TreeNode::gTaxaMgr);
	string ret;
	try {
		PhoTreesManager treesMgr(*TreeNode::gTaxaMgr);
		NxsTreesBlock treesBlock(*TreeNode::gTaxaMgr, treesMgr);
		treesBlock.AllowImplicitNames(true);
		for (StrVec::const_iterator trIt = toMerge.begin(); trIt != toMerge.end(); ++trIt)
			{
			NxsToken token(*trIt);
			treesBlock.ReadNewickRepresentationAndAddToTaxa(token);
			}
		treesBlock.EndEncountered();
		TreeShPtr tree = StrictConsensusMerger(treesMgr);
		tree->AppendNewickRepresentation(ret, true, false);
		}
	catch (...)
		{
		swap(prevGlobTaxMgr, TreeNode::gTaxaMgr);
		throw;
		}
	swap(prevGlobTaxMgr, TreeNode::gTaxaMgr);
	return ret;
	}