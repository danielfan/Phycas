#include "phycas/force_include.h"
#include <iomanip>
#include <cmath>

#include <boost/function.hpp>

#include "phycas/trees/split_manager.hpp"
#include "phycas/trees/tree_id.hpp"
#include "ncl/output/nxs_table.hpp"
#include "phycas/trees/tree_node.hpp"
#include "phycas/trees/tree.hpp"
#include "ncl/taxa/nxs_taxa_manager.hpp"
#include "phycas/trees/tree_inl.hpp"

using std::vector;
using std::ostream;
using std::ios;
using std::string;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::sqrt;
	using std::log;
#endif

// Split IDs can get very long, and are generally unintelligible to users, 
// so unless you are debugging, keep SHOW_SPLITID macro undefined
//
#undef SHOW_SPLITID

void SplitManager::BuildMajorityRuleTreeID(
  TreeID& id,		/* this id will be flushed and then filled again based on splits currently stored */
  double cutoff) const	/* only splits with this probability or higher will be included (default: 0.5) */
	{
	assert(cutoff >= 0.5);

	id.Clear();

	SplitMap::const_iterator i = smap.begin();
	for( ; i != smap.end(); ++i ) 
		{
		if ((!i->second.IsTrivial()) && i->second.GetPosteriorProbability(total) > cutoff)
				id.AddSplit(i->first);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function is called by Tree::PostorderTraverse for every node visited. SplitManager calls the PostorderTraverse 
|	member function in its RecordAllSplits function, so this function's job is to record the split defined by the node 
|	`nd'. Saves current split relative frequency in nd's plotData.pct variable, but only if plotData has been
|	previously allocated for nd.
*/
bool SplitManager::RecordSplit(
  TreeNode *nd)	/**< the node for which the split should be recorded */
	{
	if (nd->IsRoot())
		{
		return true;
		}

	bool is_trivial = false;
	if (nd->IsShootTip()) 
		{
		is_trivial = true;

		// Splits for tips just have one bit set, namely the one corresponding to their node number
		//
		nd->split.Clear();
		nd->split.SetBit(nd->nodeNum);
		}

	// If this tip is the leftmost tip of its parent, then clear the split of the parent and copy this tip's 
	// split to the parent's split. Subsequent children of this tip's parent will add their splits but will
	// not clear the parent's split.
	//
	if (nd->IsFirstPostorderChild())
		nd->par->split.Clear();
	nd->par->split.CombineWith(nd->split);

	// Add this split
	//
	Split x = nd->split;
	double splitPr = AddSplit(x, nd->GetEdgeLen(), is_trivial, false);

	// Record the current proportion of this split out of all splits in the split manager. This is useful for 
	// display purposes.
	//
	if (nd->IsPlotData())
		nd->SetPct(100.0 * splitPr);

	// Return true to indicate that traversal should continue
	//
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls PostorderTraverse member function of `tree' to walk the nodes. A functor supplied to PostorderTraverse calls 
|	the RecordSplit member function (of SplitManager) for each node visited.
*/
void SplitManager::RecordAllSplits(
  Tree &tree)	/**< Tree to traverse */
	{
	IncrementTotal();

	//POL 29 Dec 2003
	// PostorderTraverse used to take a TreeNodeListener* argument, then call the TreeNodeListener object's
	// NodeWork virtual function. This restricted each TreeNodeListener to just one NodeWork function, and the name
	// of the function gave no clue about the job it performed. Thus, the TreeNodeListener system was replaced 
	// by a functor system. Here, a functor is constructed and passed to the PostorderTraverse function. The functor
	// object is a wrapper around a member function of this SplitManager object, namely the SplitManager::RecordSplit 
	// member function, which takes a TreeNode* and returns a bool. Searching for PostorderTraverse will uncover 
	// the other instances where TreeNodeListener was formerly used. 
	//
	// The substitution of functors for TreeNodeListener seems to slow things down a tiny bit (run times are 
	// about a third of a percent longer). In one test, the new way took 659.348 sec whereas the old way required
	// 657.165 sec.
	//
	tree.PostorderTraverse(boost::function1<bool, TreeNode*>(std::bind1st(std::mem_fun(&SplitManager::RecordSplit), this))); 
	}

bool SplitManager::ComputeEntropy(double &obsEntropy)
	{
	unsigned nTaxa = Split::splitNTax;
	assert(nTaxa > 3);

	obsEntropy = 0.0;

	unsigned numSplits = GetNumSplits();
	if (numSplits == 0)
		return false;

	// Leave out the trivial splits
	//
	numSplits -= nTaxa;

	assert(total > 0L);
	double dtotal = (double)total;
	long numTrivial = 0L;
	double sum = 0.0;
	SplitMap::iterator p;
	for (p = smap.begin(); p != smap.end(); ++p)
		{
		if ((*p).second.IsTrivial())
			{
			++numTrivial;
			continue;
			}
		const unsigned frq = p->second.GetFrequency();
		assert(frq > 0);
		const double dfrq = (double) frq;
		assert(dfrq > 0);
		sum += (dfrq * log(dfrq));

		if (total - frq > 0L)
			sum += ((dtotal - dfrq)*log(dtotal - dfrq));
		}

	assert(numTrivial == (long) nTaxa);

	obsEntropy  = (double)numSplits * log(dtotal);
	obsEntropy -= (sum / dtotal);

	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Fills supplied vector with numbers of splits of each complexity, starting with the number of trivial splits
|	complexity 1) and ending with the number of splits having complexity Split::ntax/2. Returns length of the vector v.
|	This is not a cheap function to call, and is primarily intended for debugging purposes.
*/
int SplitManager::EnumSplitsByComplexity(vector<int> &v)
	{
	int vlen = (int)(Split::splitNTax / 2);

	v.erase(v.begin(), v.end());

	int k;
	for (k = 0; k < vlen; ++k)
		v.push_back(0);

	SplitMap::const_iterator p;
	for (p = smap.begin(); p != smap.end(); ++p)
		{
		const Split &x = (*p).first;
		unsigned c = x.CalcComplexity();
		assert(c > 0);
		v[c-1]++;
		}

	return vlen;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Visits nodes in post-order sequence, recomputing each node's split and the tree's TreeID. For each node visited,
|	sets branch length to the mean length of that node's split stored in the SplitManager. Also sets the split
|	proportion, which is stored in the `pct' variable in the node's `plotData' structure. The RefreshPlotData function
|	is called for the node to create the `plotData' structure if it does not exist.
*/
void SplitManager::SetBrlensFromSplits(Tree& tree) const
	{
	assert(Split::splitNTax == tree.GetNLeaves());
	tree.id.Clear();
	Length len;

	for (TreeNode* nd = tree.GetFirstPostorder(); nd != NULL; nd = nd->GetNextPostorder())
		{
		if (nd->IsShootTip()) 
			{
			nd->split.Clear();
			nd->split.SetBit(nd->nodeNum);
			Split x = nd->split;
			len.f = SplitMeanLength(x);
			nd->SetEdgeLen(len);
			nd->RefreshPlotData();
			nd->SetPct(100.0 * SplitProportion(x));

			if (nd->IsFirstPostorderChild())
				nd->par->split.Clear();
			nd->par->split.CombineWith(nd->split);
			}
		else if (!nd->IsRoot()) 
			{
			if (nd->IsFirstPostorderChild())
				nd->par->split.Clear();
			nd->par->split.CombineWith(nd->split);

			// nd->split has already been constructed by its child
			// nodes; all that is left is to add it to this Phylogeny's
			// id and to sm
			//
			Split x = nd->split;
			len.f = SplitMeanLength(x);
			tree.id.AddSplit(x);
			nd->SetEdgeLen(len);
			nd->RefreshPlotData();
			nd->SetPct(100.0 * SplitProportion(x));
			}
		}
	tree.hasEdgelens = true;
	}

unsigned SplitManager::SplitFrequency(const Split& s) const
	{
	const SplitMap::const_iterator p = smap.find( s );
	return ( p == smap.end() ? 0 : p->second.GetFrequency());
	}

double SplitManager::SplitProportion(const Split& s) const
	{
	SplitMap::const_iterator p = smap.find( s );
	return ( p == smap.end() ? 0.0 : p->second.GetPosteriorProbability(total));
	}

double SplitManager::SplitMeanLength(const Split& s) const
	{
	SplitMap::const_iterator p = smap.find( s );
	return (p == smap.end() ? 0.0 : p->second.GetMeanEdgeLength());
	}


void SplitManager::SaveSplitsAsNexusFile(
  ostream &splitsf) const
	{
	//@POL need to check for file existence!
	//
	//ofstream splitsf(fn.c_str());
	splitsf.setf(ios::fixed, ios::floatfield);
	splitsf.setf(ios::showpoint);
	splitsf << "#NEXUS\n\n";

	//@POL at the moment this (crudely) just spits out all the splits
	// that have a posterior probability greater than the specified cutoff,
	// but needs to be more sophisticated (e.g. calculating a weakly
	// compatible split system incorporating the highest weighted splits)
	//
	splitsf << "begin st_splits;\n";
	splitsf << "  dimensions ntax=" << Split::splitNTax << " nsplits=" << (unsigned)smap.size() << ";\n";
	splitsf << "  format labels weights;\n";
	splitsf << "  [properties\n";
	splitsf << "    fit=97.7665\n";
	splitsf << "    weakly compatible\n";
	splitsf << "    cyclic\n";
	splitsf << "  ;]\n";
	splitsf << "  [cycle 1 2 4 7 8 6 5 3;]\n";
	splitsf << "  matrix\n";

	unsigned k = 0;
	for (SplitMap::const_iterator i = smap.begin(); i != smap.end(); ++i) 
		{
		const double enumSplitsCutoff = 0.9; //@POL replace this enumSplitsCutoff with user setting
		const double prob = i->second.GetPosteriorProbability(total);

		if (prob > enumSplitsCutoff)
			{
			splitsf << (++k) << '\t';
			splitsf << prob << '\t';
			unsigned last = Split::splitNTax;
			for (unsigned j = 0; j < last; ++j)
				{
				if (i->first.IsBitSet(j))
					splitsf << (j+2) << ' ';
				}
			splitsf << ",\n";
			}
		}
	splitsf << ";\n";
	splitsf << "end;\n";
	}



double SplitManager::GetStdDev() const
	{
	double sd = 0.0;
	const long numSplits = (long)smap.size();
	if( numSplits > 1 ) 
		{
		double ss = 0.0;
		double s = 0.0;
		
		for(SplitMap::const_iterator i = smap.begin(); i != smap.end(); ++i ) 
			{
			const unsigned k = i->second.GetFrequency();
			s += (double) k;
			ss += (double) k * (double) k;
			}
		double n = (double)numSplits;
		double d = ss - ( s * s / n);
		double v = d / (n - 1.0);
		sd = sqrt(v);
		}
	return sd;
	}

bool SplitManager::CheckFull()
	{
	if (GetNumSplits() >= split_limit)
		{
		if (autoSplitInc)
			{
			while (GetNumSplits() >= split_limit)
				split_limit += incrSplitLimitBy;
			}
		else
			split_limit_exceeded = true;
		}
	else
		split_limit_exceeded = false;
	return split_limit_exceeded;
	}


void SplitManager::TaxaChanged(BaseTaxaManager *t, TaxaChangeType chgtype)
	{
	// NxsTaxaListener::ManagerIsDying calls this function with t == NULL
	// and chgtype == NxsTaxaListener::kTaxaMgrDestroyed. In this case, 
	// we need not do anything.
	//
	assert(t != NULL || chgtype == NxsTaxaListener::kTaxaMgrDestroyed);

	if (chgtype != NxsTaxaListener::kTaxaMgrDestroyed)
		{
		Split::CalcStatics(t->GetNumActive());
		FlushSplitsMap(); //@temp should compress
		}
	}

void SplitManager::FlushSplitsMap()
	{
	split_limit_exceeded = false;
	total = 0L;
	smap.clear();
	}

SplitManager::SplitManager(NxsTaxaManager *taxM) 
  :NxsTaxaListener(taxM),
  split_limit(500), 
  incrSplitLimitBy(100), 
  split_limit_exceeded(false), 
  autoSplitInc(true),
  total(0),
  taxaMgr(taxM)
  	{
	if (taxM != NULL)
		Split::CalcStatics(taxM->GetNumActive());
	}

SplitManager::~SplitManager()
	{
	FlushSplitsMap();
	}
	
// Returns current proportion of this split out of total splits
// a negative return value means proportion not able to be computed
// perhaps because split limit has been exceeded.
// This was done so that split proportions could be polled with
// very little overhead as a tree;s splits are being recorded.
double SplitManager::AddSplit(const Split& s, Length edgeLength, bool tip_node, bool increment_total /* = true */ )
	{
	SplitMap::iterator p = smap.find( s );
	if(p == smap.end())
		{
		//	This split is not in the map
		//
		if (split_limit_exceeded)// return -1 to indicate invalid split proportion
			return -1.0; 
		
		// Haven't yet reached split_limit, but split is not already
		// in the map, so add a new element to map
		//
		if( increment_total )
			++total;
		smap[s] = SplitInfo(tip_node, 1, edgeLength.f);
		CheckFull();
		return 1.0 / (double)total;
		}
	if( increment_total )
		++total;

	// Update existing map element
	//
	SplitInfo& si = (*p).second;
	si.frequency++;
	si.sum_lengths += edgeLength.f;
	
	return si.GetPosteriorProbability(total);
	}


#include "phycas/trees/tree_manip.hpp"
void SplitManager::BuildMajorityRuleTree(Tree &tree, double cutoffProportion, unsigned ntax) const
	{
	TreeID tid;
	BuildMajorityRuleTreeID(tid, cutoffProportion);
	TreeManip tmanip(&tree);
	tmanip.SimpleBuildTreeFromID(ntax, tid);

	//POL 22jun2005
	// Might be safer to use BuildTreeFromID rather than SimpleBuildTreeFromID
	// Need to think through the consequences of excluding taxa
	//tmanip.BuildTreeFromID(taxaMgr, tid);
	}

