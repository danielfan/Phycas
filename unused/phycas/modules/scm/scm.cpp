//#define DEBUGGING_SCM
#define DATA_FROM_DCM3       1
#define DATA_FROM_OTHERPLACE 0
#define SCM_DATA_SOURCE DATA_FROM_DCM3

#include <vector>
#include <algorithm>
#include <boost/bind.hpp>

#include "ncl/nxs_defs.hpp"
#include "phycas/force_include.h"
#include "phycas/modules/scm/scm.hpp"
#include "phycas/trees/tree_node.hpp"
#include "phycas/trees/tree.hpp"
#include "phycas/trees/draw_context.hpp"
//#include "phycas/trees/split_manager.hpp"
#include "phycas/taxa/taxa_manager.hpp"
#include "phycas/trees/trees_manager.hpp"
#include "phycas/trees/tree_manip.hpp" 
#include "ncl/output/temp_basic_output_operators.hpp"
#include "ncl/misc/nxs_index_set_output.hpp"
#include "phycas/modules/scm/scm_summary_output.hpp"
using std::vector;
using ncl::calc_mean_dbl;
using ncl::sort_to_get_median;
using std::string;
using ncl::endl;
using std::sprintf; //@VKJ 12/6/04 Added so compiles on CodeWarrior 


SCM::SCM(
  NxsOutputManager &o, 
  const PhoTreesManager &trm) 
	:outputMgr(o), 
	treesMgr(trm)
	{
	}

TreeShPtr StrictConsensusMerger(const PhoTreesManager & treesMgr)
	{
	NxsOutputManager & outputMgrRef = NxsOutputManager::GetInstance();
	SCM scm(outputMgrRef, treesMgr);
	const SCMSettings s;
	return scm.InvokeMergeTrees(s);
	}

#define USING_PAULS_MODIFICATIONS
#if defined(USING_PAULS_MODIFICATIONS)
// beginning of Paul's modifications
typedef std::set<std::string> InputTaxaSet;
typedef std::vector<InputTaxaSet> InputTaxaSetVector;

// Potentially throws these exceptions: 
//  XBadTreeDef (while building input trees)
//	XNotATaxon (while looking up root_at in taxaMgr to get its index)
//	XBadRootTaxon (while rerooting tree at root_at)
//
TreeShPtr SCM::MergeTrees(const PhoTreesManager & treesMgr, NxsOutputStream * const outStream, NxsOutputStreamWrapper * const outFilePtr)
	{
	InputTaxaSet taxaset;			// a set of taxon names
	InputTaxaSetVector taxavec;		// a vector of InputTaxaSet objects

#if defined(DEBUGGING_SCM)
	std::ofstream debugf("mergetrees.txt");
	debugf << "\nDebugging SCM::MergeTrees in file scm.cpp" << std::endl;
	debugf.flush();
#endif

	// Loop through all input trees, build each one and make an InputTaxaSet for each tree
	// containing the names of its taxa. Add this set to taxavec. When done, taxavec will
	// be used to find: 1) the common taxa (these are taxa present in all input trees); and
	// 2) the total number of taxa, num_taxa. One of the common taxa will be used to root
	// the SCM supertree.
	//
	const unsigned n = treesMgr.GetNumTrees();
	for (unsigned i = 0; i < n; ++i)
		{
		Tree t;

		// Build tree i using stored tree description
		//
		const FullTreeDescription &desc = treesMgr.GetTree(i);
		t.BuildTreeFromDescription(desc); 

#if defined(DEBUGGING_SCM)
		debugf << "\nInput tree " << i << " (of " << n << "):\n-->" << desc.newick << std::endl;
		debugf.flush();
#endif

		// Insert the taxa names appearing in every tree into InputTaxaSet taxaset, 
		// and then push_back this taxaset into InputTaxaSetVector taxavec.
// *******************************************************************************
// NOTE: using the names is probably not the way to go, should translate the names
// into taxon indices using the taxa manager, then use the indices hereafter
// *******************************************************************************
		//
		for (TreeNode *nd = t.GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
			{
			if (nd->IsLeafOrRoot())
				{
				taxaset.insert(nd->GetName());
#if defined(DEBUGGING_SCM)
				debugf << "  inserting into set: " << nd->GetName() << std::endl;
				debugf.flush();
#endif
				}
			}
		taxavec.push_back(taxaset);
		taxaset.clear();
		}

	// Find the taxa common to all input trees and store in variable commonTaxa.
	//
	// Example: 4 taxa sets from 4 input trees
	// In this example, I am presenting the indices of taxa whereas in Phycas these would 
	// be taxon names, the translation being done when the trees were read in. This example
	// is from MTH_22jun2005.
	//   taxaset 1 -->                   30 32    37 38 39 40 41 42       45
	//   taxaset 2 -->                   30 32          39          43 44 45 46 47 48 49
	//   taxaset 3 -->  5 13             30 32          39                45
	//   taxaset 4 -->       26 27 28 29 30 32 33       39                45
	// Loop starts with:
	//   commonTaxa: 30 32 37 38 39 40 41 42 45
	//   temp: <empty>
	// First time through the loop:
	//   commonTaxa: 30 32 37 38 39 40 41 42       45
	//   *it:        30 32       39          43 44 45 46 47 48 49
	//   temp: 30 32 39 45
	// Second time through the loop:
	//   commonTaxa:      30 32 39 45
	//   *it:        5 13 30 32 39 45
	//   temp: 30 32 39 45
	// Third time through the loop:
	//   commonTaxa:      30 32    39 45
	//   *it: 26 27 28 29 30 32 33 39 45
	//   temp: 30 32 39 45
	// Afterwards:
	//   commonTaxa: 30 32 39 45
	//
	InputTaxaSet temp;
	InputTaxaSetVector::iterator it = taxavec.begin();
	InputTaxaSet commonTaxa = *it;
	for (++it; it != taxavec.end(); ++it)
		{
		std::set_intersection(commonTaxa.begin(), commonTaxa.end(), it->begin(), it->end(), inserter(temp, temp.begin())); 
		commonTaxa = temp;
		temp.clear();
		}


#if defined(DEBUGGING_SCM)
{
	unsigned num_common_taxa = (unsigned)commonTaxa.size();
	debugf << "\nNumber of common taxa: " << num_common_taxa << "\nHere is a list of just the common taxa:\n  ";
	InputTaxaSet::const_iterator debug_it = commonTaxa.begin();
	for (; debug_it != commonTaxa.end(); ++debug_it)
		{
		debugf << *debug_it << " ";
		}
	debugf << std::endl;
	debugf.flush();
}
#endif

	// Find set of taxa present in any of the input trees and store in variable allTaxa.
	//
	// Example: 4 taxa sets from 4 input trees (MTH_22jun2005 continued)
	//   taxaset 1 -->                   30 32    37 38 39 40 41 42       45
	//   taxaset 2 -->                   30 32          39          43 44 45 46 47 48 49
	//   taxaset 3 -->  5 13             30 32          39                45
	//   taxaset 4 -->       26 27 28 29 30 32 33       39                45
	// Loop starts with:
	//   allTaxa: 30 32 37 38 39 40 41 42 45
	//   temp: <empty>
	// First time through the loop:
	//   allTaxa: 30 32 37 38 39 40 41 42       45
	//   *it:     30 32       39          43 44 45 46 47 48 49
	//   temp:    30 32 37 38 39 40 41 42 43 44 45 46 47 48 49
	// Second time through the loop:
	//   allTaxa:   30 32 37 38 39 40 41 42 43 44 45 46 47 48 49
	//   *it:  5 13 30 32       39                45
	//   temp: 5 13 30 32 37 38 39 40 41 42 43 44 45 46 47 48 49
	// Third time through the loop:
	//   allTaxa: 5 13             30 32    37 38 39 40 41 42 43 44 45 46 47 48 49
	//   *it:          26 27 28 29 30 32 33       39                45
	//   temp:    5 13 26 27 28 29 30 32 33 37 38 39 40 41 42 43 44 45 46 47 48 49
	// Afterwards:
	//   allTaxa: 5 13 26 27 28 29 30 32 33 37 38 39 40 41 42 43 44 45 46 47 48 49
	//
	it = taxavec.begin();
	InputTaxaSet allTaxa = *it;
	temp.clear();
	for (++it; it != taxavec.end(); ++it)
		{
		std::set_union(allTaxa.begin(), allTaxa.end(), it->begin(), it->end(), inserter(temp, temp.begin()));
		allTaxa = temp;
		temp.clear();
		}

	unsigned num_taxa_in_input_trees = (unsigned)allTaxa.size();
	unsigned num_taxa = treesMgr.taxaMgr.GetNumActive(); 
	
#if defined(DEBUGGING_SCM)
{
	debugf << "\nTotal number of active taxa: " << num_taxa;
	debugf << "\nTotal number of taxa in all input trees: " << num_taxa_in_input_trees;
	debugf << "\nHere is a list of all taxa in at least one input tree:";
	debugf << "\n  ";
	InputTaxaSet::const_iterator debug_it = allTaxa.begin();
	for (; debug_it != allTaxa.end(); ++debug_it)
		{
		debugf << *debug_it << " ";
		}
	debugf << std::endl;
	debugf.flush();
}
#endif

	typedef std::map<Split, Split> TreeMap;
	typedef std::pair<Split, Split> TreeMapPair;
	typedef std::pair<TreeMap::iterator, bool> TreePairInsertionResult;

	typedef std::map<unsigned, Split> SplitMap;
	typedef std::pair<unsigned, Split> SplitMapPair;
	typedef std::pair<SplitMap::iterator, bool> SplitPairInsertionResult;

	// Each element of the mask vector is a split with bits set for every taxon in one of the input trees
	//
	std::vector<Split> mask;

	// The id vector is used to store the treeid of every input tree
	//
	std::vector<TreeID> id;

	// The polytomy vector is used to store splits representing polytomies 
	// found in any input tree
	//
	std::vector<Split> polytomy;

	// Will root all input trees at the first common taxon
	//
	std::string root_at = *(commonTaxa.begin());

#if defined(DEBUGGING_SCM)
	debugf << "\nInput trees will be rooted at this common taxon: " << root_at << std::endl;
	debugf.flush();
#endif

	// This loop builds the mask, id and polytomy vectors
	//
	for (unsigned i = 0; i < n; i++)
		{
		// Build each tree in PhoTaxaManager treesMgr, using the specified root taxon to reroot the tree.
		// 
		Tree t;
		const FullTreeDescription &desc = treesMgr.GetTree(i);
		try 
			{
			t.BuildTreeFromString(desc.newick.c_str(), false, false); 
			}
		catch (XBadTreeDef & x)
			{
			//@POL should let this exception fall through to a higher level
			char message[256];
			sprintf(message, "The following problem was encountered in trying to build a\ntree from the specified tree description:\n  %s", x.msg.c_str());
			}

		// Throw XNotATaxon exception if specified root taxon is not found in the tree.
		//
		unsigned root_at_index = treesMgr.taxaMgr.FindIndex(root_at); 
		if (root_at_index == UINT_MAX)
			{
			throw XNotATaxon();
			}

		// Throw XBadRootTaxon exception if specified root taxon turns out not to be common to all input trees.
		//
		bool ok = t.RerootAtTip(root_at_index); 
		if (!ok)
			{
			throw XBadRootTaxon();
			}

		// Get the treeid of every input tree, and push_back it into vector<TreeID> id
		//
		t.DebugCheckTreeStructure();

		TreeID idt = t.GetID();

		id.push_back(idt);

		// Check tree topology to find out these polytomy splits
		//
		for (const TreeNode *nd = t.GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
			{
			if (!(nd->IsLeafOrRoot()))
				{
				unsigned nchildren = nd->CountChildren();
				if (nchildren != 2)
					{
					const TreeNode *tempnode = nd->GetLeftChild();
					polytomy.push_back(tempnode->GetSplit());
					for (unsigned j = 1; j < nchildren; j++)
						{
						tempnode = tempnode->GetRightSib();
						polytomy.push_back(tempnode->GetSplit());
						}
					}
				}
			}

		// Create a Split that has bits set for every taxon in this input tree, then
		// add it to vector<Split> mask
		//
		Split tmask;
		for (const TreeNode *nd = t.GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
			{
			if (nd->IsLeafOrRoot())
				{
				unsigned k = nd->GetNodeNumber();
				tmask.SetBit(k);
				}
			}
		mask.push_back(tmask);
		}

#if defined(DEBUGGING_SCM)
{
	debugf << "\nHere are the " << (unsigned)mask.size() << " splits in the mask vector:\n";
	for (unsigned x = 0; x < Split::GetNTaxa(); ++x)
		{
		debugf << (char)('a' + (x % 26));
		}
	debugf << std::endl;
	for (std::vector<Split>::const_iterator debug_it = mask.begin(); debug_it != mask.end(); ++debug_it)
		{
		debugf << debug_it->CreatePatternRepresentation() << "\n";
		}
	debugf << std::endl;
	debugf.flush();

	debugf << "\nHere are the " << (unsigned)polytomy.size() << " splits in the polytomy vector:\n";
	for (unsigned x = 0; x < Split::GetNTaxa(); ++x)
		{
		debugf << (char)('a' + (x % 26));
		}
	debugf << std::endl;
	for (std::vector<Split>::const_iterator debug_it = polytomy.begin(); debug_it != polytomy.end(); ++debug_it)
		{
		debugf << debug_it->CreatePatternRepresentation() << "\n";
		}
	debugf << std::endl;
	debugf.flush();

	debugf << "\nHere are the splits in each of the " << (unsigned)id.size() << " input trees:\n";
	unsigned tree_number = 0;
	for (std::vector<TreeID>::const_iterator id_it = id.begin(); id_it != id.end(); ++id_it)
		{
		debugf << "\n  Splits in tree number " << tree_number << ":\n    ";
		for (unsigned x = 0; x < Split::GetNTaxa(); ++x)
			{
			debugf << (char)('a' + (x % 26));
			}
		debugf << std::endl;
		const SplitSet &split_set = id_it->GetSplitSet();
		for (std::set<Split>::const_iterator debug_it = split_set.begin(); debug_it != split_set.end(); ++debug_it)
			{
			debugf << "    " << debug_it->CreatePatternRepresentation() << "\n";
			}
		tree_number++;
		}
	debugf << std::endl;
	debugf.flush();
}
#endif

	// Create common_mask with bits set only for common taxa in all trees
	//
	std::vector<Split>::iterator mask_it = mask.begin();
	Split common_mask = *mask_it;
	for (mask_it = mask.begin(); mask_it != mask.end(); mask_it++)
		{
		common_mask.IntersectWith(*mask_it);
		}

#if defined(DEBUGGING_SCM)
{
	debugf << "\nHere is the common_mask:" << std::endl;
	for (unsigned x = 0; x < Split::GetNTaxa(); ++x)
		{
		debugf << (char)('a' + (x % 26));
		}
	debugf << std::endl;
	debugf << common_mask.CreatePatternRepresentation() << std::endl;
	debugf.flush();
}
#endif

	// Check polytomy splits, and find those orphan taxa that do not directly connect with common taxa
	//
	std::set<unsigned> polytomytaxa;
	for (std::vector<Split>::iterator it = polytomy.begin(); it != polytomy.end(); it++)
		{
		//std::string s = it->CreatePatternRepresentation();
		//std::cout << "split = " << s << std::endl;
		Split temps = (*it) & common_mask;
		if (temps.CountOnBits() == 0)
			{
			for (unsigned j = 0; j < num_taxa; j++)
				{
				if ((*it).IsBitSet(j))
					polytomytaxa.insert(j);
				}
			}
		}

	// Create SplitSet all_splits to store the splits of all input trees
	// Create vector<SplitSet> subtree to store the substructure of every input tree. Here substructure means that the subtree structure of 
	//		every input tree, which is some splits only composed by common taxa ( here we use the common_mask to represent them ). The reason 
	//  	we do this is to find the intersection of all elements of vector<SplitSet> subtree, that is, the common substructure of all input trees.
	//      we name it SplitSet common_subtree here. By using SplitSet common_subtree, we can check whether the start point (common subtree) is 
	//      inconsistent or not. We use the term "inconsistent" here to represent there is polytomy existed in the common subtree. Since we always
	//		use DCM3 first to decompose the initial dataset, it seems that we need not worry about inconsistency of common subtree.
	//
	std::vector<SplitSet> subtree;
	SplitSet all_splits;
	for (std::vector<TreeID>::iterator id_it = id.begin(); id_it != id.end(); id_it++)
		{
		SplitSet tempss = id_it->GetSplitSet();
		//Split alltaxa;
		//for (SplitSet::iterator ssit = tempss.begin(); ssit != tempss.end(); ssit++)
		//	{
		//	std::cout << (*ssit).CreatePatternRepresentation() << std::endl;
		//	alltaxa |= (*ssit);
		//	}
		//	tempss.insert(alltaxa);
		//std::cout << tempss.size() << std::endl;
		SplitSet tempsub;
		for (SplitSet::iterator tempss_it = tempss.begin(); tempss_it != tempss.end(); tempss_it++)
			{
			// Get SplitSet all_splits here
			//
			all_splits.insert(*tempss_it);
			Split temps = (*tempss_it) & common_mask;
			if ((temps.CountOnBits() != (unsigned)0) && (temps.CountOnBits() != (unsigned)1))
				{
				tempsub.insert(temps);
				}
			}
		// Get vector<SplitSet> subtree here
		//
		subtree.push_back(tempsub);
		}

	// Get SplitSet common_subtree here. We first try to find the element of vector<SplitSet> subtree that has the smallest size,
	// and then SplitSet common_subtree
	//
	SplitSet common_subtree;

	// Define the possible maximum number of subtree_size
	//
	unsigned subtree_size = num_taxa_in_input_trees - 2;	
	std::vector<SplitSet>::iterator temp_subtree_it, subtree_it;
	for (subtree_it = subtree.begin(); subtree_it != subtree.end(); subtree_it++)
		{
		unsigned size = (unsigned)subtree_it->size();
		// Find the element of vector<SplitSet> subtree that has the smallest size
		//
		if (size <= subtree_size) //POL Daniel, we found a case where max subtree size equals num_taxa_in_input_trees
			{
			subtree_size = size;
			temp_subtree_it = subtree_it;
			}
		}

	// Use SplitSet common_subtree_candidate that has the smallest size to find splits belonging to SplitSet common_subtree
	//
	SplitSet common_subtree_candidate = *temp_subtree_it;
	for (SplitSet::iterator csc_it = common_subtree_candidate.begin(); csc_it != common_subtree_candidate.end(); csc_it++)
		{
		bool found = true;
		for (std::vector<SplitSet>::iterator subtree_it = subtree.begin(); subtree_it != subtree.end(); subtree_it++)
			{
			if (subtree_it->count(*csc_it) == 0)
				found = false;
			}
		if (found)
			{
			common_subtree.insert(*csc_it);
			//std::cout << (*csc_it).CreatePatternRepresentation() << std::endl;
			//std::cout << common_mask.CreatePatternRepresentation() << std::endl;
			}
		}

	// Check whether the common subtree is consistent or not
	//
	bool CommonSubtreeConsistent = true;
	//std::cout << common_subtree.size();
    //std::cout << ((common_subtree.begin())+1)->CreatePatternRepresentation() << std::endl;
	if ((common_mask.CountOnBits() - 2) != (common_subtree.size()))
		{
		CommonSubtreeConsistent = false;
		std::cout << "Oops! The CommonSubtree of trees in the input file is inconsistent. Be careful!" << std::endl;
		}

	// For every Orphan Taxon, calculates the count of every splits stored in the SplitSet all_splits once. Iteratively do that for all orphan taxa.
	//  
	//	(1)	if find the count number of splits is equal to 0 or 1, then insert the corresponding Split into SplitSet zero_splitset or one_splitset
	//      respectively. Meanwhile, insert the split whose count number is 1 into TreeMap tmapCandidate <key: Split & common_mask, value: Split> 
	//
	//	(2)	if find the lowest count number of all splits in SplitSet all_splits is more than 1, then insert the corresponding Split into 
	//      TreeMap tmapCandidate and insert the corresponding Split into SplitMap smapCandidate <key: Orphan Taxon number, value: Split>
	//
	//  (3) in the above two operations, if find any two Splits have the same key when inserting into tmapCandidate, then uses CombineWith
	//      method on the two Splits
	//
	//
	TreeMap tmapCandidate;
	SplitMap smapCandidate;
	SplitSet zero_splitset, one_splitset;

	Split s_for_count,orphantax_mask;
	
	unsigned common_size = common_mask.CountOnBits();
	std::vector<TreeMap> treemap(common_size - 2); // tmap_two_count, tmap_three_count;
	for (unsigned k = 0; k < num_taxa; k++)
		{
		if (!common_mask.IsBitSet(k))	// if TRUE means this position or bit corresponds to an orphan taxa
			{
			orphantax_mask = common_mask;
			orphantax_mask.SetBit(k);	// in orphantax_mask the bits are set for this orphan taxa and all common taxa
			SplitSet::iterator temp_ssit;
			std::vector< vector<Split> > count_vector(common_size); // zero_count, one_count, two_count, three_count;
			for (SplitSet::iterator ssit = all_splits.begin(); ssit != all_splits.end(); ssit++)
				{
				if ((*ssit & orphantax_mask).IsBitSet(k))	// if TRUE means this split have this orphan taxa set, and should be used for counting the count number
					{										 
					s_for_count = (*ssit & common_mask);	// when counting the count number, we just consider those common taxa bits 
					unsigned count = s_for_count.CountOnBits();
					(count_vector[count]).push_back(*ssit);
					}
				}
			for (unsigned i = 0; i < common_size; i++)
				{
				std::sort((count_vector[i]).begin(),(count_vector[i]).end());
				}

			if (!(count_vector[0]).empty())
				{
				for (std::vector<Split>::iterator vs_it = (count_vector[0]).begin(); vs_it != (count_vector[0]).end(); vs_it++)
					{
					zero_splitset.insert(*vs_it);
					}
				}
			if (!(count_vector[1]).empty())
				{
				for (std::vector<Split>::iterator vs_it = (count_vector[1]).begin(); vs_it != (count_vector[1]).end(); vs_it++)
					{
					one_splitset.insert(*vs_it);
					TreePairInsertionResult result = tmapCandidate.insert(TreeMapPair((*vs_it & common_mask), *vs_it));
					if (result.second == false)
						{
						((result.first)->second).CombineWith(*vs_it);
						}
					}
				}
			if ((count_vector[1]).empty())// && (!two_count.empty() || !three_count.empty()))
				{
				bool found = false;
				unsigned temp;
				for (unsigned i = 2; i < common_size; i++)
					{
					if (!(count_vector[i]).empty())
						{
						found = true;
						temp = i;
						break;
						}
					}
				if (found)
					{
					TreePairInsertionResult result = tmapCandidate.insert(TreeMapPair(((count_vector[temp])[0] & common_mask), (count_vector[temp])[0]));
					if (result.second == false)
						{
						((result.first)->second).CombineWith((count_vector[temp])[0]);
						}
					if (polytomytaxa.count(k) == 0)
						{
						Split tempsplit;
						tempsplit.SetBit(k);
						TreePairInsertionResult result2 = (treemap[temp-2]).insert(TreeMapPair((count_vector[temp])[0], tempsplit));
						if (result2.second == false)
							{
							((result2.first)->second).CombineWith(tempsplit);
							}
						}
					}
				}
			}
		}

#if defined(DEBUGGING_SCM)
{
	debugf << "\nSplits composing zero_splitset:\n";
	for (unsigned x = 0; x < Split::GetNTaxa(); ++x)
		{
		debugf << (char)('a' + (x % 26));
		}
	debugf << std::endl;
	for (SplitSet::iterator debug_it = zero_splitset.begin(); debug_it != zero_splitset.end(); ++debug_it)
		{
		debugf << debug_it->CreatePatternRepresentation() << "\n";
		}
	debugf << std::endl;
	debugf.flush();

	debugf << "\nSplits composing one_splitset:\n";
	for (unsigned x = 0; x < Split::GetNTaxa(); ++x)
		{
		debugf << (char)('a' + (x % 26));
		}
	debugf << std::endl;
	for (SplitSet::iterator debug_it = one_splitset.begin(); debug_it != one_splitset.end(); ++debug_it)
		{
		debugf << debug_it->CreatePatternRepresentation() << "\n";
		}
	debugf << std::endl;
	debugf.flush();
}
#endif

		// Create the SplitSet ssSuper to store the merged tree's Splits
		// Adds Splits to ssSuper:
		//  (1a) Deal with zero_splitset
		//  (1b) Deal with one_splitset
		//  (2) Splits' key in tmapCandidate is equal to 1
		//  (3) OR Splits whose (key == (key & Split in common_subtree)), then adds them to ssSuper. If no found, adds this common_subtree split to ssSuper.
		//  (4) Make copies of Splits in ssSuper which have the same key(after masking with common_mask) as Splits in smapCandidate, 
		//  	unset k bit, then adds them to ssSuper
		//  (5) OR all of the splits, which now are contained in the SplitSet ssSuper 
		//
		SplitSet ssSuper;

		// (1a)
		//
		for (SplitSet::iterator ssit = zero_splitset.begin(); ssit != zero_splitset.end(); ssit++)
			{
			//std::cout << (*ssit).CreatePatternRepresentation() << std::endl;
			ssSuper.insert(*ssit);
			}

#if defined(DEBUGGING_SCM)
{
	debugf << "\nssSuper after inserting splits in zero_splitset (1a)" << std::endl;
	for (unsigned x = 0; x < Split::GetNTaxa(); ++x)
		{
		debugf << (char)('a' + (x % 26));
		}
	debugf << std::endl;
	for (SplitSet::iterator ssit = ssSuper.begin(); ssit != ssSuper.end(); ++ssit)
		{
		debugf << ssit->CreatePatternRepresentation() << "\n";
		}
	debugf << std::endl;
	debugf.flush();
}
#endif

		// (1b)
		//  The operation here is similar to (1a)
		//
		while (one_splitset.size() != 0)	// the one_splitset will be empty at last
			{
			SplitSet tempset;
			std::vector<Split> vs;
			Split key = (*(one_splitset.begin()) & common_mask);

			// Find out those splits that share the same KEY ( split AND common_mask ),
			// and then insert them into SplitSet tempset and vector<Split> vs
			//
			for (SplitSet::iterator ssit = one_splitset.begin(); ssit != one_splitset.end(); ssit++)
				{
				if (key == (*ssit & common_mask))
					{
					tempset.insert(*ssit);
					vs.push_back(*ssit);
					}
				}

			// According to SplitSet tempset, delete those splits from SplitSet one_splitset
			//
			for (SplitSet::iterator ssit = tempset.begin(); ssit != tempset.end(); ssit++)
				{
				SplitSet::iterator found = one_splitset.find(*ssit);
				if (found != one_splitset.end())
					one_splitset.erase(found);
				}

			// Sort vector<Split> vs, let it from the biggest one to the smallest one
			//
			std::sort(vs.begin(),vs.end());

			// Check whether the former one in vector<Split> vs could contain the latter one iteratively
			// if TRUE means those splits could compose a whole clade, and then insert those splits into SplitSet ssSuper
			//
			bool whole_clade = true;
			for (std::vector<Split>::iterator it = vs.begin(); it != (vs.end() -1); it++)
				{
				if (*it != (*it & *(it + 1)))
					whole_clade = false;
				}
			if (whole_clade)
				{
				for (std::vector<Split>::const_iterator it = vs.begin(); it != vs.end(); it++)
					{
					ssSuper.insert(*it);
					}
				}
			}

#if defined(DEBUGGING_SCM)
{
	debugf << "\nssSuper after inserting splits in one_splitset (1b)\n";
	for (unsigned x = 0; x < Split::GetNTaxa(); ++x)
		{
		debugf << (char)('a' + (x % 26));
		}
	debugf << std::endl;
	for (SplitSet::iterator ssit = ssSuper.begin(); ssit != ssSuper.end(); ++ssit)
		{
		debugf << ssit->CreatePatternRepresentation() << "\n";
		}
	debugf << std::endl;
	debugf.flush();
}
#endif

		// (2)
		// Find out those splits from TreeMap tmapCandidate whose KEYs' CountOnBits are equal to 1, 
		// and then insert them into SplitSet ssSuper
		//		
		for (TreeMap::iterator tmap_it = tmapCandidate.begin(); tmap_it != tmapCandidate.end(); tmap_it++)
			{
			unsigned temp = ((tmap_it->first).CountOnBits());
			if (temp == 1)
				ssSuper.insert(tmap_it->second);
			}

		// (3)
		// Use the splits of common_subtree iteraively to find out splits from TreeMap tmapCandidate. 
		// The criterion is the KEY of split of TreeMap tmapCandidate can be contained by the split of common_subtree.
		// Next, OR those splits together (including the split of common_subtree), then adds it to ssSuper. 
		// If not find such kind of splits, add this common_subtree split to ssSuper.
		//
		for (SplitSet::iterator ssit = common_subtree.begin(); ssit != common_subtree.end(); ssit++)
			{
			Split temp = *ssit;
			for (TreeMap::iterator tmap_it = tmapCandidate.begin(); tmap_it != tmapCandidate.end(); tmap_it++)
				{
				Split mask = (tmap_it->first);
				if (mask == (mask & *ssit))
					{
					//std::cout << (tmap_it->second).CreatePatternRepresentation() << std::endl;
					temp = (temp | (tmap_it->second));
					}
				}
			ssSuper.insert(temp);
			}

		// (4)
		std::vector< vector<Split> > ssSuperSplit(common_size - 2); //, ssSuperSplit_commonTaxa3;
		for (SplitSet::iterator ssit = ssSuper.begin(); ssit != ssSuper.end(); ssit++)
			{
			Split mask = (*ssit & common_mask);
			if (mask.CountOnBits() > 1)
				{
				(ssSuperSplit[mask.CountOnBits() - 2]).push_back(*ssit);
				}
			}
		
		for (SplitSet::iterator ssit = common_subtree.begin(); ssit != common_subtree.end(); ssit++)
			{
			unsigned countbits = (*ssit).CountOnBits();
			if(!(ssSuperSplit[countbits - 2]).empty())
				{
				std::vector<Split> SuperSplitVec;
				for (std::vector<Split>::iterator vecit = (ssSuperSplit[countbits - 2]).begin(); vecit != (ssSuperSplit[countbits - 2]).end(); vecit++)
					{
					if (((*vecit) & common_mask) == (*ssit))
						SuperSplitVec.push_back(*vecit);
					}
				std::sort(SuperSplitVec.begin(),SuperSplitVec.end());
				if (!(treemap[countbits - 2]).empty())
					{
					std::vector<Split> tmapSplitVec;
					for (TreeMap::iterator tmap_it = (treemap[countbits - 2]).begin(); tmap_it != (treemap[countbits - 2]).end(); tmap_it++)
						{
						if (((tmap_it->first)&common_mask) == (*ssit))
                            tmapSplitVec.push_back(tmap_it->first);
						}
					std::sort(tmapSplitVec.begin(),tmapSplitVec.end());
					bool compatible = true;
						//mth changed July7 to fix crash in CIPRES
#					if 0 
						
						for (std::vector<Split>::iterator it = tmapSplitVec.begin(); it != (tmapSplitVec.end() - 1); it++)
							{
							if (*it != (*it & *(it+1)))
								{
								compatible = false;
								//std::cout << "splits in tmap_three_count are uncompatible!" << std::endl;
								}
							}
#					else
						if (!tmapSplitVec.empty())
							{
							std::vector<Split>::const_iterator it = tmapSplitVec.begin();
							std::vector<Split>::const_iterator prevIt = it;
							for (++it; it != tmapSplitVec.end(); ++it, ++prevIt)
								{
								if (*prevIt != ((*prevIt) & (*it)))
									{
									compatible = false;
									break;
									}
								}
							}
#					endif
					if (compatible)
						{
						Split split_from_ssSuper = SuperSplitVec[SuperSplitVec.size() - 1]; 
						for (std::vector<Split>::reverse_iterator r_it = tmapSplitVec.rbegin(); r_it != tmapSplitVec.rend(); r_it++)
							{
							TreeMap::iterator tmap_it2 = (treemap[countbits - 2]).find(*r_it);
							split_from_ssSuper ^= (tmap_it2->second);
							ssSuper.insert(split_from_ssSuper);
							}
						}
					else
						{
						Split split_from_ssSuper = SuperSplitVec[SuperSplitVec.size() - 1];
						Split temp;
						for (std::vector<Split>::iterator it = tmapSplitVec.begin(); it != tmapSplitVec.end(); it++)
							{
							TreeMap::iterator tmap_it2 = (treemap[countbits - 2]).find(*it);
							temp |= (tmap_it2->second);
							}
						split_from_ssSuper ^= temp;
						ssSuper.insert(split_from_ssSuper);
						}
					}
				}
			}

		// (5)
		// OR together all of the splits which now are contained in the SplitSet ssSuper
		//
		Split the_last_one;
		for (SplitSet::iterator ssit = ssSuper.begin(); ssit != ssSuper.end(); ssit++)
			{
			the_last_one = (the_last_one | (*ssit));
			}
		ssSuper.insert(the_last_one);

		// Bulids and displays the merged tree
		//
		TreeID idSCM;
		for (SplitSet::iterator ssit = ssSuper.begin(); ssit != ssSuper.end(); ssit++)
			{
			idSCM.AddSplit(*ssit);
			}

#if defined(DEBUGGING_SCM)
{
/* For MTH_22jun2005 example:
	Splits in ssSuper composing the supertree:
	abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvw
	----*-------*------------------------------------  1
	-------------------------*--*-------------------- 17
	----*-------*------------------*-----------------  2
	-------------------------*--*---*---------------- 18 
	------------------------------------**-----------  3
	------------------------------------**-*---------  4
	------------------------------------****---------  5
	------------------------------------****-*-------  6
	------------------------------------******-------  7
	------------------------------------------*-*----  8
	------------------------------------------***----  9
	----------------------------------------------**- 10
	---------------------------------------------***- 11
	---------------------------------------------**** 12
	------------------------------------------******* 13
	------------------------------------************* 14
	----*-------*------------------*----************* 15
	----*-------*-------------*----*----************* 16
	----*-------*------------**-*--**---************* 19
	----*-------*------------****--**---************* 20
*/
	debugf << "\nSplits in ssSuper composing the supertree:" << std::endl;
	for (unsigned x = 0; x < Split::GetNTaxa(); ++x)
		{
		debugf << (char)('a' + (x % 26));
		}
	debugf << std::endl;
	for (SplitSet::iterator ssit = ssSuper.begin(); ssit != ssSuper.end(); ssit++)
		{
		debugf << ssit->CreatePatternRepresentation() << "\n";
		}
	debugf << std::endl;
	debugf.flush();
	debugf.close();
}
#endif

		unsigned root_at_index = treesMgr.taxaMgr.FindIndex(root_at); 
		if (root_at_index == UINT_MAX)
			{
			throw XNotATaxon();
			}

	TreeShPtr tSCM = TreeShPtr(new Tree());
	TreeManip tm(tSCM.get());
	tm.SimpleBuildTreeFromID(num_taxa,idSCM,root_at_index);	//POL 22jun2005
	tSCM->RefreshID();
	tSCM->DebugCheckTreeStructure();
	return tSCM;
	}
// end of Paul's modifications
#else  // #if defined(USING_PAULS_MODIFICATIONS)
// beginning of Daniel's original code
TreeShPtr SCM::MergeTrees(const PhoTreesManager & treesMgr, NxsOutputStream * const outStream, NxsOutputStreamWrapper * const outFilePtr)
	{
	const unsigned n = treesMgr.GetNumTrees();

	// Find commonTaxa, use it to root the tree 
	// Find allTaxa, use it to define num_taxa
	//
	typedef std::set<std::string> InputTaxaSet;
	typedef std::vector<InputTaxaSet> InputTaxaSetVector;
	InputTaxaSet taxaset;
	InputTaxaSetVector taxavec;

	// Here the n equals to the number of trees in the input file
	//
	for (unsigned i = 0; i < n; ++i)
		{
		Tree t;

		// Build tree i using stored tree description
		//
		const FullTreeDescription &desc = treesMgr.GetTree(i);
		t.BuildTreeFromDescription(desc); 

		// Insert the taxa names appearing in every tree into InputTaxaSet taxaset, 
		// and then push_back this taxaset into InputTaxaSetVector taxavec
		//
		for (TreeNode *nd = t.GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
			{
			if (nd->IsLeafOrRoot())
				{
				taxaset.insert(nd->GetName());
				}
			}
		taxavec.push_back(taxaset);
		taxaset.clear();
		}

	// By dealing with every element ( here is InputTaxaSet taxaset ) in InputTaxaSetVector taxavec,
	// and using set_intersection operation iteratively to find the common taxa of all trees
	//
	InputTaxaSet commonTaxa = *(taxavec.begin());
	for (InputTaxaSetVector::iterator it = taxavec.begin(); it != (taxavec.end() - 1); it++)
		{
		InputTaxaSet temp = commonTaxa;
		commonTaxa.clear();
		std::set_intersection(temp.begin(), temp.end(), (it+1)->begin(), (it+1)->end(), inserter(commonTaxa, commonTaxa.begin())); 
		}

	// By dealing with every element ( here is InputTaxaSet taxaset ) in InputTaxaSetVector taxavec,
	// and using set_union operation iteratively to find all taxa of all trees
	// Create a set ( InputTaxaSet allTaxa ) to contain all taxa
	//
	InputTaxaSet allTaxa;
	for (InputTaxaSetVector::iterator it = taxavec.begin(); it != taxavec.end(); it++)
		{
		InputTaxaSet temp = allTaxa;
		allTaxa.clear();
		std::set_union(temp.begin(), temp.end(), it->begin(), it->end(), inserter(allTaxa, allTaxa.begin()));
		}

	// Now get the number of all taxa
	//
	unsigned num_taxa = (unsigned)allTaxa.size();
	
	typedef std::map<Split, Split> TreeMap;
	typedef std::pair<Split, Split> TreeMapPair;
	typedef std::pair<TreeMap::iterator, bool> TreePairInsertionResult;

	typedef std::map<unsigned, Split> SplitMap;
	typedef std::pair<unsigned, Split> SplitMapPair;
	typedef std::pair<SplitMap::iterator, bool> SplitPairInsertionResult;

	// Create vector<TreeID> id, and use it to store the treeid of every tree
	// Create vector<Split> mask, and use it to store the Split of every tree
	// 
	// I guess I need to make some explanations here. Every element of vector<Split> mask is a Split indeed, 
	// which has bits set for every taxon in the input tree, and then push_back it into vector<Split> mask
	//
	std::vector<Split> mask;
	std::vector<Split> polytomy;
	std::vector<TreeID> id;
//	SplitManager splitmgr((NxsTaxaManager*)&taxaMgr);
	// Will root all input trees at the first common taxon
	//
	std::string root_at = *(commonTaxa.begin());

	// Here the n equals to the number of trees in the input file
	//
	for (unsigned i = 0; i < n; i++)
		{
		// Build every tree existed in PhoTaxaManager treesMgr iteratively, and use the specified root taxon to reroot the tree.
		// 
		Tree t;
		const FullTreeDescription & desc = treesMgr.GetTree(i);
		try 
			{
			t.BuildTreeFromString(desc.newick.c_str(), false, false); 
			}
		catch (XBadTreeDef & x)
			{
			char message[256];
			sprintf(message, "The following problem was encountered in trying to build a\ntree from the specified tree description:\n  %s", x.msg.c_str());
			}

		// Throw XNotATaxon exception if specified root taxon is not found in the tree.
		//
		unsigned root_at_index = treesMgr.taxaMgr.FindIndex(root_at); 
		if (root_at_index == UINT_MAX)
			{
			throw XNotATaxon();
			}

		// Throw XBadRootTaxon exception if specified root taxon turns out not to be common to all input trees.
		//
		bool ok = t.RerootAtTip(root_at_index); 
		if (!ok)
			{
			throw XBadRootTaxon();
			}

		// Get the treeid of every input tree, and push_back it into vector<TreeID> id
		//
		t.DebugCheckTreeStructure();

		TreeID idt = t.GetID();

		id.push_back(idt);

		// Check tree topology to find out these polytomy splits
		//
		for (const TreeNode *nd = t.GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
			{
			if (!(nd->IsLeafOrRoot()))
				{
				unsigned nchildren = nd->CountChildren();
				if (nchildren != 2)
					{
					const TreeNode *tempnode = nd->GetLeftChild();
					polytomy.push_back(tempnode->GetSplit());
					for (unsigned j = 1; j < nchildren; j++)
						{
						tempnode = tempnode->GetRightSib();
						polytomy.push_back(tempnode->GetSplit());
						}
					}
				}
			}

		// Create a Split that has bits set for every taxon in this input tree, then
		// add it to vector<Split> mask
		//
		Split tmask;
		for (const TreeNode *nd = t.GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
			{
			if (nd->IsLeafOrRoot())
				{
				unsigned k = nd->GetNodeNumber();
				tmask.SetBit(k);
				}
			}
		mask.push_back(tmask);
		}

	// Create common_mask with bits set only for common taxa in all trees
	//
	std::vector<Split>::iterator mask_it = mask.begin();
	Split common_mask = *mask_it;
	for (mask_it = mask.begin(); mask_it != mask.end(); mask_it++)
		{
		common_mask.IntersectWith(*mask_it);
		}


	// Check polytomy splits, and find out these orphan taxa that do not directly connect with common taxa
	//
	std::set<unsigned> polytomytaxa;
	for (std::vector<Split>::iterator it = polytomy.begin(); it != polytomy.end(); it++)
		{
		//std::string s = it->CreatePatternRepresentation();
		//std::cout << "split = " << s << std::endl;
		Split temps = (*it) & common_mask;
		if (temps.CountOnBits() == 0)
			{
			for (unsigned j = 0; j < num_taxa; j++)
				{
				if ((*it).IsBitSet(j))
					polytomytaxa.insert(j);
				}
			}
		}

	// Create SplitSet all_splits to store the splits of all input trees
	// Create vector<SplitSet> subtree to store the substructure of every input tree. Here substructure means that the subtree structure of 
	//		every input tree, which is some splits only composed by common taxa ( here we use the common_mask to represent them ). The reason 
	//  	we do this is to find the intersection of all elements of vector<SplitSet> subtree, that is, the common substructure of all input trees.
	//      we name it SplitSet common_subtree here. By using SplitSet common_subtree, we can check whether the start point (common subtree) is 
	//      inconsistent or not. We use the term "inconsistent" here to represent there is polytomy existed in the common subtree. Since we always
	//		use DCM3 first to decompose the initial dataset, it seems that we need not worry about inconsistency of common subtree.
	//
	std::vector<SplitSet> subtree;
	SplitSet all_splits;
	for (std::vector<TreeID>::iterator id_it = id.begin(); id_it != id.end(); id_it++)
		{
		SplitSet tempss = id_it->GetSplitSet();
		//Split alltaxa;
		//for (SplitSet::iterator ssit = tempss.begin(); ssit != tempss.end(); ssit++)
		//	{
		//	std::cout << (*ssit).CreatePatternRepresentation() << std::endl;
		//	alltaxa |= (*ssit);
		//	}
		//	tempss.insert(alltaxa);
		//std::cout << tempss.size() << std::endl;
		SplitSet tempsub;
		for (SplitSet::iterator tempss_it = tempss.begin(); tempss_it != tempss.end(); tempss_it++)
			{
			// Get SplitSet all_splits here
			//
			all_splits.insert(*tempss_it);
			Split temps = (*tempss_it) & common_mask;
			if ((temps.CountOnBits() != (unsigned)0) && (temps.CountOnBits() != (unsigned)1))
				{
				tempsub.insert(temps);
				}
			}
		// Get vector<SplitSet> subtree here
		//
		subtree.push_back(tempsub);
		}

	// Get SplitSet common_subtree here. We first try to find the element of vector<SplitSet> subtree that has the smallest size,
	// and then SplitSet common_subtree
	//
	SplitSet common_subtree;
	// Define the possible maximum number of subtree_size
	//
	unsigned subtree_size = num_taxa - 2;	
	std::vector<SplitSet>::iterator temp_subtree_it, subtree_it;
	for (subtree_it = subtree.begin(); subtree_it != subtree.end(); subtree_it++)
		{
		unsigned size = (unsigned)subtree_it->size();
		// Find the element of vector<SplitSet> subtree that has the smallest size
		//
		if (size <= subtree_size) //POL Daniel, we found a case where max subtree size equals num_taxa
			{
			subtree_size = size;
			temp_subtree_it = subtree_it;
			}
		}

	// Use SplitSet common_subtree_candidate that has the smallest size to find splits belonging to SplitSet common_subtree
	//
	SplitSet common_subtree_candidate = *temp_subtree_it;
	for (SplitSet::iterator csc_it = common_subtree_candidate.begin(); csc_it != common_subtree_candidate.end(); csc_it++)
		{
		bool found = true;
		for (std::vector<SplitSet>::iterator subtree_it = subtree.begin(); subtree_it != subtree.end(); subtree_it++)
			{
			if (subtree_it->count(*csc_it) == 0)
				found = false;
			}
		if (found)
			{
			common_subtree.insert(*csc_it);
			//std::cout << (*csc_it).CreatePatternRepresentation() << std::endl;
			//std::cout << common_mask.CreatePatternRepresentation() << std::endl;
			}
		}

	// Check whether the common subtree is consistent or not
	//
	bool CommonSubtreeConsistent = true;
	//std::cout << common_subtree.size();
    //std::cout << ((common_subtree.begin())+1)->CreatePatternRepresentation() << std::endl;
	if ((common_mask.CountOnBits() - 2) != (common_subtree.size()))
		{
		CommonSubtreeConsistent = false;
		std::cout << "Oops! The CommonSubtree of trees in the input file is inconsistent. Be careful!" << std::endl;
		}

	// For every Orphan Taxon, calculates the count of every splits stored in the SplitSet all_splits once. Iteratively do that for all orphan taxa.
	//  
	//	(1)	if find the count number of splits is equal to 0 or 1, then insert the corresponding Split into SplitSet zero_splitset or one_splitset
	//      respectively. Meanwhile, insert the split whose count number is 1 into TreeMap tmapCandidate <key: Split & common_mask, value: Split> 
	//
	//	(2)	if find the lowest count number of all splits in SplitSet all_splits is more than 1, then insert the corresponding Split into 
	//      TreeMap tmapCandidate and insert the corresponding Split into SplitMap smapCandidate <key: Orphan Taxon number, value: Split>
	//
	//  (3) in the above two operations, if find any two Splits have the same key when inserting into tmapCandidate, then uses CombineWith
	//      method on the two Splits
	//
	//
	TreeMap tmapCandidate;
	SplitMap smapCandidate;
	SplitSet zero_splitset, one_splitset;

	Split s_for_count,orphantax_mask;
	
	unsigned common_size = common_mask.CountOnBits();
	std::vector<TreeMap> treemap(common_size - 2); // tmap_two_count, tmap_three_count;
	for (unsigned k = 0; k < num_taxa; k++)
		{
		if (!common_mask.IsBitSet(k))	// if TRUE means this position or bit corresponds to an orphan taxa
			{
			orphantax_mask = common_mask;
			orphantax_mask.SetBit(k);	// in orphantax_mask the bits are set for this orphan taxa and all common taxa
			SplitSet::iterator temp_ssit;
			std::vector< vector<Split> > count_vector(common_size); // zero_count, one_count, two_count, three_count;
			for (SplitSet::iterator ssit = all_splits.begin(); ssit != all_splits.end(); ssit++)
				{
				if ((*ssit & orphantax_mask).IsBitSet(k))	// if TRUE means this split have this orphan taxa set, and should be used for counting the count number
					{										 
					s_for_count = (*ssit & common_mask);	// when counting the count number, we just consider those common taxa bits 
					unsigned count = s_for_count.CountOnBits();
					(count_vector[count]).push_back(*ssit);
					}
				}
			for (unsigned i = 0; i < common_size; i++)
				{
				std::sort((count_vector[i]).begin(),(count_vector[i]).end());
				}

			if (!(count_vector[0]).empty())
				{
				for (std::vector<Split>::iterator vs_it = (count_vector[0]).begin(); vs_it != (count_vector[0]).end(); vs_it++)
					{
					zero_splitset.insert(*vs_it);
					}
				}
			if (!(count_vector[1]).empty())
				{
				for (std::vector<Split>::iterator vs_it = (count_vector[1]).begin(); vs_it != (count_vector[1]).end(); vs_it++)
					{
					one_splitset.insert(*vs_it);
					TreePairInsertionResult result = tmapCandidate.insert(TreeMapPair((*vs_it & common_mask), *vs_it));
					if (result.second == false)
						{
						((result.first)->second).CombineWith(*vs_it);
						}
					}
				}
			if ((count_vector[1]).empty())// && (!two_count.empty() || !three_count.empty()))
				{
				bool found = false;
				unsigned temp;
				for (unsigned i = 2; i < common_size; i++)
					{
					if (!(count_vector[i]).empty())
						{
						found = true;
						temp = i;
						break;
						}
					}
				if (found)
					{
					TreePairInsertionResult result = tmapCandidate.insert(TreeMapPair(((count_vector[temp])[0] & common_mask), (count_vector[temp])[0]));
					if (result.second == false)
						{
						((result.first)->second).CombineWith((count_vector[temp])[0]);
						}
					if (polytomytaxa.count(k) == 0)
						{
						Split tempsplit;
						tempsplit.SetBit(k);
						TreePairInsertionResult result2 = (treemap[temp-2]).insert(TreeMapPair((count_vector[temp])[0], tempsplit));
						if (result2.second == false)
							{
							((result2.first)->second).CombineWith(tempsplit);
							}
						}
					}
				}
			}
		}

	
		// Create the SplitSet ssSuper to store the merged tree's Splits
		// Adds Splits to ssSuper:
		//  (1a) Deal with zero_splitset
		//  (1b) Deal with one_splitset
		//  (2) Splits' key in tmapCandidate is equal to 1
		//  (3) OR Splits whose (key == (key & Split in common_subtree)), then adds them to ssSuper. If no found, adds this common_subtree split to ssSuper.
		//  (4) Make copies of Splits in ssSuper which have the same key(after masking with common_mask) as Splits in smapCandidate, 
		//  	unset k bit, then adds them to ssSuper
		//  (5) OR all of the splits, which now are contained in the SplitSet ssSuper 
		//
		SplitSet ssSuper;

		// (1a)
		//
		for (SplitSet::iterator ssit = zero_splitset.begin(); ssit != zero_splitset.end(); ssit++)
			{
			//std::cout << (*ssit).CreatePatternRepresentation() << std::endl;
			ssSuper.insert(*ssit);
			}

		// (1b)
		//  The operation here is similar to (1a)
		//
		while (one_splitset.size() != 0)	// the one_splitset will be empty at last
			{
			SplitSet tempset;
			std::vector<Split> vs;
			Split key = (*(one_splitset.begin()) & common_mask);

			// Find out those splits that share the same KEY ( split AND common_mask ),
			// and then insert them into SplitSet tempset and vector<Split> vs
			//
			for (SplitSet::iterator ssit = one_splitset.begin(); ssit != one_splitset.end(); ssit++)
				{
				if (key == (*ssit & common_mask))
					{
					tempset.insert(*ssit);
					vs.push_back(*ssit);
					}
				}

			// According to SplitSet tempset, delete those splits from SplitSet one_splitset
			//
			for (SplitSet::iterator ssit = tempset.begin(); ssit != tempset.end(); ssit++)
				{
				SplitSet::iterator found = one_splitset.find(*ssit);
				if (found != one_splitset.end())
					one_splitset.erase(found);
				}

			// Sort vector<Split> vs, let it from the biggest one to the smallest one
			//
			std::sort(vs.begin(),vs.end());

			// Check whether the former one in vector<Split> vs could contain the latter one iteratively
			// if TRUE means those splits could compose a whole clade, and then insert those splits into SplitSet ssSuper
			//
			bool whole_clade = true;
			for (std::vector<Split>::iterator it = vs.begin(); it != (vs.end() -1); it++)
				{
				if (*it != (*it & *(it + 1)))
					whole_clade = false;
				}
			if (whole_clade)
				{
				for (std::vector<Split>::const_iterator it = vs.begin(); it != vs.end(); it++)
					{
					ssSuper.insert(*it);
					}
				}
			}


		// (2)
		// Find out those splits from TreeMap tmapCandidate whose KEYs' CountOnBits are equal to 1, 
		// and then insert them into SplitSet ssSuper
		//		
		for (TreeMap::iterator tmap_it = tmapCandidate.begin(); tmap_it != tmapCandidate.end(); tmap_it++)
			{
			unsigned temp = ((tmap_it->first).CountOnBits());
			if (temp == 1)
				ssSuper.insert(tmap_it->second);
			}

		// (3)
		// Use the splits of common_subtree iteraively to find out splits from TreeMap tmapCandidate. 
		// The criterion is the KEY of split of TreeMap tmapCandidate can be contained by the split of common_subtree.
		// Next, OR those splits together (including the split of common_subtree), then adds it to ssSuper. 
		// If not find such kind of splits, add this common_subtree split to ssSuper.
		//
		for (SplitSet::iterator ssit = common_subtree.begin(); ssit != common_subtree.end(); ssit++)
			{
			Split temp = *ssit;
			for (TreeMap::iterator tmap_it = tmapCandidate.begin(); tmap_it != tmapCandidate.end(); tmap_it++)
				{
				Split mask = (tmap_it->first);
				if (mask == (mask & *ssit))
					{
					//std::cout << (tmap_it->second).CreatePatternRepresentation() << std::endl;
					temp = (temp | (tmap_it->second));
					}
				}
			ssSuper.insert(temp);
			}

		// (4)
		std::vector< vector<Split> > ssSuperSplit(common_size - 2); //, ssSuperSplit_commonTaxa3;
		for (SplitSet::iterator ssit = ssSuper.begin(); ssit != ssSuper.end(); ssit++)
			{
			Split mask = (*ssit & common_mask);
			if (mask.CountOnBits() > 1)
				{
				(ssSuperSplit[mask.CountOnBits() - 2]).push_back(*ssit);
				}
			}
		
		for (SplitSet::iterator ssit = common_subtree.begin(); ssit != common_subtree.end(); ssit++)
			{
			unsigned countbits = (*ssit).CountOnBits();
			if(!(ssSuperSplit[countbits - 2]).empty())
				{
				std::vector<Split> SuperSplitVec;
				for (std::vector<Split>::iterator vecit = (ssSuperSplit[countbits - 2]).begin(); vecit != (ssSuperSplit[countbits - 2]).end(); vecit++)
					{
					if (((*vecit) & common_mask) == (*ssit))
						SuperSplitVec.push_back(*vecit);
					}
				std::sort(SuperSplitVec.begin(),SuperSplitVec.end());
				if (!(treemap[countbits - 2]).empty())
					{
					std::vector<Split> tmapSplitVec;
					for (TreeMap::iterator tmap_it = (treemap[countbits - 2]).begin(); tmap_it != (treemap[countbits - 2]).end(); tmap_it++)
						{
						if (((tmap_it->first)&common_mask) == (*ssit))
                            tmapSplitVec.push_back(tmap_it->first);
						}
					std::sort(tmapSplitVec.begin(),tmapSplitVec.end());
					bool compatible = true;
					for (std::vector<Split>::iterator it = tmapSplitVec.begin(); it != (tmapSplitVec.end() - 1); it++)
						{
						if (*it != (*it & *(it+1)))
							{
							compatible = false;
							//std::cout << "splits in tmap_three_count are uncompatible!" << std::endl;
							}
						}
					if (compatible)
						{
						Split split_from_ssSuper = SuperSplitVec[SuperSplitVec.size() - 1]; 
						for (std::vector<Split>::reverse_iterator r_it = tmapSplitVec.rbegin(); r_it != tmapSplitVec.rend(); r_it++)
							{
							TreeMap::iterator tmap_it2 = (treemap[countbits - 2]).find(*r_it);
							split_from_ssSuper ^= (tmap_it2->second);
							ssSuper.insert(split_from_ssSuper);
							}
						}
					else
						{
						Split split_from_ssSuper = SuperSplitVec[SuperSplitVec.size() - 1];
						Split temp;
						for (std::vector<Split>::iterator it = tmapSplitVec.begin(); it != tmapSplitVec.end(); it++)
							{
							TreeMap::iterator tmap_it2 = (treemap[countbits - 2]).find(*it);
							temp |= (tmap_it2->second);
							}
						split_from_ssSuper ^= temp;
						ssSuper.insert(split_from_ssSuper);
						}
					}
				}
			}

		// (5)
		// OR together all of the splits which now are contained in the SplitSet ssSuper
		//
		Split the_last_one;
		for (SplitSet::iterator ssit = ssSuper.begin(); ssit != ssSuper.end(); ssit++)
			{
			the_last_one = (the_last_one | (*ssit));
			}
		ssSuper.insert(the_last_one);

		// Bulids and displays the merged tree
		//
		TreeID idSCM;
		for (SplitSet::iterator ssit = ssSuper.begin(); ssit != ssSuper.end(); ssit++)
			{
			idSCM.AddSplit(*ssit);
			}

		unsigned root_at_index = treesMgr.taxaMgr.FindIndex(root_at); 
		if (root_at_index == UINT_MAX)
			{
			throw XNotATaxon();
			}

	TreeShPtr tSCM = TreeShPtr(new Tree());
	TreeManip tm(tSCM.get());
	//tm.SimpleBuildTreeFromID(num_taxa,idSCM,root_at_index);	//POL 22jun2005
	tm.BuildTreeFromID(&treesMgr.taxaMgr, idSCM, root_at_index);	//POL 22jun2005
	tSCM->RefreshID();
	tSCM->DebugCheckTreeStructure();
	return tSCM;
	}
// end of Daniel's original code
#endif

CmdResult SCM::HandleSCM(SCMSettings * s)
	{
	try
		{
		if (s != NULL)
			InvokeMergeTrees(*s);
		}
	catch (XBadTreeInSCM &)
		{
		NxsErrorStream * errStream  = outputMgr.GetErrorStreamPtr();
		if (errStream != NULL)
			*errStream << "Illegal tree found in SCM operation (Trees should be cleared from memory)."<< ncl::endl;
		return kCmdFailedSilent;
		}
	return kCmdSucceeded;
	}

TreeShPtr	SCM::InvokeMergeTrees(const SCMSettings & s)
	{
	NxsOutputStream * outStream = outputMgr.GetOutputStreamPtr();
	// Output log file column header label strings
	NxsOutputStreamWrapperShPtr outFileShPtr = outputMgr.GetGenericOutputStreamShPtr(s.scmfile);
	NxsOutputStreamWrapper * outFilePtr = outFileShPtr.get();
	//const PhoTreesManager & treesMgr = PhoTreesManager::GetInstance();

	TreeShPtr tree = MergeTrees(treesMgr, outStream, outFilePtr);

	const unsigned n = treesMgr.GetNumTrees();
	SCMSummary summary;
	if (outStream != NULL || outFilePtr != NULL)
		{
		tree->AppendNewickRepresentation(summary.supertreeDescription, true, false);
		if  (outStream != NULL)
			{
			*outStream << "\nStrict consensus merger tree based on " << n << " input trees:\n" << endl;
			ASCIIDrawContext ascii(*outStream);
			tree->Draw(ascii);
			}
		if (outFilePtr != NULL)
			Emit<kTabSeparatedTabularLabelsOutStyle, SCMSummary *>(*outFilePtr, &summary);
		}
	return tree;
	}

