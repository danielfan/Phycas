#include "phycas/force_include.h"
#include <vector>
#include <boost/bind.hpp>

#include "ncl/nxs_defs.hpp"
#include "phycas/modules/phylodiversity/pd.hpp"
#include "phycas/trees/tree_node.hpp"
#include "phycas/trees/tree.hpp"
#include "phycas/trees/draw_context.hpp"
#include "phycas/taxa/taxa_manager.hpp"
#include "phycas/trees/trees_manager.hpp"
#include "ncl/output/temp_basic_output_operators.hpp"
#include "ncl/misc/nxs_index_set_output.hpp"
#include "phycas/modules/phylodiversity/phylodiversity_summary_output.hpp"
using std::vector;
using ncl::calc_mean_dbl;
using ncl::sort_to_get_median;
using std::string;
using ncl::endl;

void PhylodiversitySummary::calculateStatisticsFromCumulatives()
	{
	assert(cumI > 0.0);
	assert(cumTotal > 0.0);
#   if defined(NEW_PD_WAY)
		if (cumE > cumI)
			cumE = cumI;
		pEI = cumE/cumI;
		pIT = cumI/cumTotal;
		pET = cumE/cumTotal;
#   else
		pEI = (cumI - cumEmin)/cumI;
		pIT = (cumTotal - cumI)/cumTotal;
		pET = (cumTotal - cumEmax)/cumTotal;
#   endif
	}

CmdResult Phylodiversity::HandlePhylodiversity(PhylodiversitySettings *s)
	{
	if (!taxaMgr || ! treesMgr)
		return kCmdFailedSilent; // @ error message 
	const string focalTaxsetName = GetCapitalized(s->taxonset);
	const NxsIndexSet * const focalTaxa = taxaMgr->GetSet(focalTaxsetName.c_str());
	NxsOutputManager & outputMgrRef = NxsOutputManager::GetInstance();
	if (focalTaxa == NULL)
		{
		if (outputMgrRef.GetErrorStreamPtr() != NULL)
			*outputMgrRef.GetErrorStreamPtr() << "Error: taxon set named " << focalTaxsetName << " does not exist" << endl;
		return kCmdFailedSilent;
		}
	NxsOutputStream * outStream = outputMgrRef.GetOutputStreamPtr();
	NxsOutputStream * returnValStream = outputMgrRef.GetReturnValueStreamPtr();
	if (outStream != NULL && taxaMgr)
		Emit<kVerboseOutStyle>(*outStream, NxsIndexSetOutput(*focalTaxa, *taxaMgr, kTaxaSetOutput));
	
	const Split selSplit(*focalTaxa);
	
	const unsigned n = treesMgr->GetNumTrees();
	if (outStream != NULL)
		*outStream << "\nNumber of trees: " << n << endl;
	
	NxsOutputStreamWrapperShPtr pdFileShPtr = outputMgrRef.GetGenericOutputStreamShPtr(s->pdOut);
	NxsOutputStreamWrapper * outFilePtr = pdFileShPtr.get();
	if (outFilePtr)
		Emit<kTabSeparatedTabularLabelsOutStyle, PhylodiversitySummary *>(*outFilePtr, (PhylodiversitySummary *)NULL);
	Tree t; 
	for (unsigned i = 0; i < n; ++i)
		{
		t.BuildTreeFromDescription(treesMgr->GetTree(i)); 
		if (outStream != NULL && (i+1) % 100 == 0)
			*outStream << "\n  " << (i+1) << " trees processed..." << endl; //use endl because we are intercalating outStream and returnValStream
		PhylodiversitySummary pdSummary = ProcessTree(selSplit, t);
		pdSummary.treeNumber = i + 1;
		if (returnValStream != NULL)
			Emit<kVerboseOutStyle>(*returnValStream, pdSummary) << endl; //use endl because we are intercalating outStream and returnValStream
		if (outFilePtr != NULL)
			EmitGeneric<kTabSeparatedTabularOutStyle>(*outFilePtr, pdSummary) << endl;
		}
	return kCmdSucceeded;
	}		

/*----------------------------------------------------------------------------------------------------------------------
|	ProcessTree is called by Phylodiversity::ProcessTreesInMemory, which supplies a Tree built from a tree definition 
|	currently stored in memory (read in originally from a tree file by the gettrees command).
*/
PhylodiversitySummary Phylodiversity::ProcessTree(
  const Split & selSplit,
  const Tree  & tree)	
	{
	PhylodiversitySummary pdSummary;
	const TreeNode * r = tree.GetRoot();
	Split S = r->GetSplit();
	S.InvertSplit();
	if (S.SubsumedIn(selSplit))
		{	//@pol: currently, first taxon (i.e. the tip serving as the root) cannot be in the selected set
		pdSummary.isValid = false;
		return pdSummary;
		}
	VecDbl selectedTipLength, unselectedTipLength;
	for (const TreeNode * nd = tree.GetFirstPostorder(); !nd->IsRoot(); nd = nd->GetNextPostorder())
		{
		const bool is_pendant_edge = ! (nd->EdgeIsInternal());
		const double brlen = nd->GetFltEdgeLen();
		pdSummary.cumTotal += brlen;
			// Check for inclusion in inclusive phylodiversity (I)
			// ---------------------------------------------------
			// Want to add this node's edge length to cumI if at least one taxon on each side of
			// the split associated with this branch can be found within selSplit.
		S = nd->GetSplit();
		bool I_include = (nd->IsLeafOrRoot() && S.SubsumedIn(selSplit));
		if (!I_include && nd->EdgeIsInternal())
			{
			S.IntersectWith(selSplit);
			if (S.CountOnBits() > 0)
				{
				// If here, means at least one bit set in S was also set in selSplit
				// Now see if at least one bit in complement of S is set in selSplit
				//
				S = nd->GetSplit();
				S.InvertSplit();
				S.IntersectWith(selSplit);
				if (S.CountOnBits() > 0)
					I_include = true;
				}
			}
			/*  Comment from previous section of code that was #if 0'd out of compilation
			14 May 2004
				Going back to defining I as the traditional PD measure, but E is still
				the same as the old Emax. This is the only section of code that needs to
				be changed, and all that we need to do is igorne this section, so it has
				been conditionally compiled out, leaving the other NEW_PD_WAY sections
				intact.

			13 May 2004
				Now calculating I by including the subtending edge, so I is now max(E)
				E is the same as the old Emax, also including subtending edge
				Old Emin is now not calculated
				Old PEI = Emin/I, new version is PEI = E/I
				Old PET = Emax/T, new version is PET = E/T (no change in meaning)
				Old PIT = I/T, new version is still PIT = I/T (but I has changed meaning)
				To test whether current node's edge is the subtending edge in calculation of I,
				need to pass the following tests:
				1) no selected taxa are on one side of the branch (lower side)
				2) all selected taxa are on the other side (higher side)
				3) all children of current node must be selected (prevents all edges 
				back to the root from being included in the selected set of edges)

				Suppose taxa 2, 3 and 5 are in the selected set, nd is the current node,
				and nd's split has bits for taxa 1, 2, 3, 4, and 5 set:

						  9876543210
				selSplit  0000101100
					 ~nd  1111000001
						  ----------
				bitwise & 0000000000 <- if this equals 0, satisfies condition 1

						  9876543210
				selSplit  0000101100
					  nd  0000111110
						  ----------
				bitwise & 0000101100 <- if this equals selSplit, satisfies condition 2
			#endif
			*/

			// Create plotData structure for this node if it doesn't exist
			// and set plotData->selected
		const_cast<TreeNode *>(nd)->RefreshPlotData(); //@@ perhaps plot data should be mutable? (or we should not declare that the Tree arg to ProcessTree is const)
		const_cast<TreeNode *>(nd)->SetSelectionStatus(I_include);
		if (I_include)	
			pdSummary.cumI += brlen;
			// Check for inclusion in exclusive phylodiversity (Emin and/or Emax)
			// ------------------------------------------------------------------
			// Want to add this node's edge length to cumEmax if this split is subsumed 
			// within selSplit
		S = nd->GetSplit(); // refresh split rep from node
		const bool isInExclusivePD = S.SubsumedIn(selSplit);
		if (isInExclusivePD)
			{
#			if defined(NEW_PD_WAY)
				pdSummary.cumE += brlen;
#			else
				pdSummary.cumEmax += brlen;
					// If computing Emin, do not include internal edge if parent  edge is not also subsumed
				if (is_pendant_edge || nd->GetParent()->GetSplit().SubsumedIn(selSplit)) 						
					pdSummary.cumEmin += brlen;
#			endif //defined(NEW_PD_WAY)
			}
		if (is_pendant_edge)
			{
			if (isInExclusivePD)
				selectedTipLength.push_back(brlen);
			else
				unselectedTipLength.push_back(brlen);
			}
		}
	assert(selectedTipLength.size() > 0);
	assert(unselectedTipLength.size() > 0);
	pdSummary.avgSelTipLen = ncl::calc_mean_dbl(selectedTipLength);
	pdSummary.medSelTipLen = ncl::sort_to_get_median(selectedTipLength);
	pdSummary.avgUnselTipLen = ncl::calc_mean_dbl(unselectedTipLength);
	pdSummary.medUnselTipLen = ncl::sort_to_get_median(unselectedTipLength);
	pdSummary.calculateStatisticsFromCumulatives();
	return pdSummary;
	}


