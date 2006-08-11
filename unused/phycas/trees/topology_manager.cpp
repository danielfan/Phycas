#include "phycas/force_include.h"
#include <cmath>
#include <iomanip>
#include "phycas/trees/topology_manager.hpp"
#include "phycas/trees/tree.hpp"
#include "ncl/taxa/nxs_taxa_manager.hpp"
#include "ncl/output/nxs_table.hpp"
#include "phycas/misc/msg_exception.hpp"
using std::make_pair;
using std::ostream;
using std::setw;
using std::pair;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::log;
#endif

void TopologyManager::AddTopology(const Tree& tree)
	{
	AddTopology(tree.GetID());
	}

void TopologyManager::TaxaChanged(
  BaseTaxaManager *taxaMgr, 
  TaxaChangeType)
	{
	nTaxa = taxaMgr->GetNumActive();
	FlushTopoMap(); //@should only clear after a few taxa change types
	}

TopologyManager::TopologyManager(
  NxsTaxaManager *tm) 
  :NxsTaxaListener(tm),
  topo_limit(1000U), 
  incrTopoLimitBy(100U), 
  topo_limit_exceeded(false), 
  autoTopoInc(true),
  totalAddAttempts(0U),
  sortVectDirty(true),
  taxaMgr(tm) //POL 22jun2005
 	{
	nTaxa = 0;
	if (tm != NULL)
		nTaxa = tm->GetNumActive();
	}

TopologyManager::~TopologyManager()
	{
	FlushTopoMap();
	}


void TopologyManager::CreateSortVect() const
	{
	// Create sortmultimap
	//
	unsigned topNum = 0;
	TopoSortMultimap sortmultimap;
	TopoMap::const_iterator j = tmap.begin();
	for (; j != tmap.end(); ++j, ++topNum)
		sortmultimap.insert(make_pair(j->second.frequency, topNum));
	
	// Create sortvect;
	//
	sortvect.clear();
	sortvect.reserve(sortmultimap.size());
	TopoSortMultimap::reverse_iterator jj = sortmultimap.rbegin();
	for (; jj != sortmultimap.rend(); ++jj)
		sortvect.push_back(jj->second);
	sortVectDirty = false;
	}
	
// Returns index of element with ith highest frequency. First call slowest because vector has to be sorted; 
// after that, requires only random access of vector element.
// 
unsigned TopologyManager::GetSorted(unsigned i) const
	{
	assert( i < GetNumTopologiesStored());
	if (sortVectDirty)
		CreateSortVect();
	return sortvect[i];
	}

void TopologyManager::AddTopology(const TreeID& t, double wt /* = -1.0*/)
	{
	sortVectDirty = true;
	TopoMap::iterator p = tmap.find( t );
	++totalAddAttempts;
	if( p == tmap.end() && topo_limit_exceeded )
		{
		// Too many topologies already stored, from now on 
		// we increment totalAddAttempts but do not add any more
		// elements to map
		//
		return;
		}
	if( p == tmap.end() ) 
		{
		// Haven't yet reached topo_limit, so add a 
		// new element to map
		//
		tmap[t] = TopoInfo(wt, 1);
		CheckFull();
		}
	else 
		p->second.frequency++;
	}

bool TopologyManager::CheckFull()
	{
	if (GetNumTopologiesStored() >= topo_limit)
		{
		if (autoTopoInc)
			{
			while (GetNumTopologiesStored() >= topo_limit)
				topo_limit += incrTopoLimitBy;
			}
		else
			topo_limit_exceeded = true;
		}
	else
		topo_limit_exceeded = false;
	return topo_limit_exceeded;
	}

void TopologyManager::FlushTopoMap()
	{
	sortVectDirty = true;
	topo_limit_exceeded = false;
	totalAddAttempts = 0;
	tmap.erase( tmap.begin(), tmap.end() );
	}

void TopologyManager::SetTopologyLimit( unsigned limit )
	{
	topo_limit = limit;
	CheckFull();
	}

void TopologyManager::SetTopoAutoInc(bool autoInc, unsigned incr)
	{
	autoTopoInc = autoInc;
	incrTopoLimitBy = incr;
	}

void TopologyManager::SaveTopoMap( ostream& out, int width /* = 12 */ ) const
	{
	out << setw(width);
	out << "\nTopology frequencies:\n";
	out << setw(width) << "Freq." << "  Phylogeny ID\n";
	
	unsigned numTopologies = (unsigned)tmap.size();
	if( numTopologies > 0U ) 
		{
		for(TopoMap::const_iterator i = tmap.begin(); i != tmap.end(); ++i ) 
			{
			out << setw(width) << i->second.frequency << "  ";
			out << i->first.CreateIdRepresentation() << '\n';
			}
		}
	else
		out << "  No toplogies have been stored\n";
		
	out << "Total topologies = " << GetNumTopologiesStored() << '\n';
	out << "Total distinct   = " << numTopologies << '\n';
	}

// TreeIDs can get very long, and are generally unintelligible to users, 
// so unless you are debugging, keep SHOW_TREEID macro undefined
//
#undef SHOW_TREEID


unsigned TopologyManager::ShowTopoMapAsTable( NxsTable& table, bool show_wts /* = false */, unsigned width /* = 12 */, unsigned precision /* = 6 */ ) const
	{
	unsigned numTopologies = (unsigned)tmap.size();
	if( numTopologies == 0 )
		return 0;
		
	assert( totalAddAttempts > 0U );
	
	NxsTableCell::def_width     = width;
	NxsTableCell::def_precision = precision;

	table.Reset();

	table.AddString("Topology");
	table.SetLeftMargin();
	table.AddString("Freq.");
	if (show_wts)
		table.AddString("Weight");
	table.AddString("Prob.");
#	if defined(SHOW_TREEID)
		table.AddString("Id");
#	endif
	table.SetRightMargin();
	table.HyphenLine();
	table.SetTopMargin();

	unsigned k = 0;
	
	for(TopoMap::const_iterator i = tmap.begin(); i != tmap.end(); ++i ) 
		{
		table.AddUInt(++k);
		const double frq = i->second.frequency;
		const double prob = frq / (double)totalAddAttempts;
		table.AddDouble( frq );
		if (show_wts)
			table.AddDouble(i->second.weight);
		table.AddDouble( prob );
#		if defined(SHOW_TREEID)
			table.AddString( i->first.CreateIdRepresentation() );
#		endif
		}

	table.SetBottomMargin();
		
	return numTopologies;
	}

unsigned TopologyManager::GetFrequency(unsigned which) const
	{
	assert( which < GetNumTopologiesStored() );
	unsigned k = 0;
	for(TopoMap::const_iterator i = tmap.begin(); i != tmap.end(); ++i, ++k)
		{
		if (k == which)
			return i->second.frequency;
		}
	assert(0);
	THROW_MSG_EXCEPTION("Requesting Frequency from the ToplogyManager with an illegal index");
	}

double TopologyManager::TopologyWeight(const TreeID& t) const
	{
	TopoMap::const_iterator p = tmap.find(t);
	return (p != tmap.end() ? p->second.weight : 0.0);
	}

unsigned TopologyManager::TopologyFrequency(const TreeID& t ) const
	{
	TopoMap::const_iterator p = tmap.find( t );
	return (p != tmap.end() ? p->second.frequency : 0);
	}

pair<TreeID, TopoInfo> TopologyManager::GetTopology(unsigned which) const
	{
	unsigned k = 0;
	assert(which < GetNumTopologiesStored());
	TopoMap::const_iterator i;
	for(i = tmap.begin(); i != tmap.end(); ++i, ++k)
		{
		if (k == which)
			return (*i);
		}
	assert(0);
	THROW_MSG_EXCEPTION("Requesting iterator from the ToplogyManager with an illegal index");
	}

TreeID TopologyManager::GetID( unsigned which ) const
	{
	assert( which < GetNumTopologiesStored() );
	unsigned k = 0;
	for(TopoMap::const_iterator i = tmap.begin(); i != tmap.end(); ++i, ++k)
		{
		if (k == which)
			return i->first;
		}
	assert(0);
	THROW_MSG_EXCEPTION("Requesting TreeID from the ToplogyManager with an illegal index");
	}

// Note: involves copying TreeID object twice (once into which, then again
// upon function return), but this seems preferable to doing a second loop
// Assumes TopologyManager has at least one topology stored.
//
TreeID TopologyManager::GetIDofBest() const
	{
	assert(GetNumTopologiesStored() > 0);

	TopoMap::const_iterator i = tmap.begin();
	TopoMap::const_iterator best = i++;
	for(; i != tmap.end(); ++i ) 
		{
		if( best->second.frequency < i->second.frequency ) 
			best = i;
		}

	return best->first;
	}

bool TopologyManager::ComputeEntropy(double &obsEntropy, double &maxEntropy, double &KLDistance) const 
	{
	obsEntropy = 0.0;
	maxEntropy = 0.0;
	KLDistance = 0.0;

	unsigned numTopos = GetNumTopologiesStored();
	if (numTopos == 0L)
		return false;

	double sumFrqLnFrq = 0.0;
	double N = 0.0;
	for (unsigned i = 0; i < numTopos; ++i)
		{
		unsigned k = GetSorted(i);
		unsigned frq = GetFrequency(k);
		assert(frq > 0L);
		N += (double)frq;

		double dfrq = (double)frq;
		sumFrqLnFrq += (dfrq * log(dfrq));
		}

	// Maximum entropy is average entropy one would calculate if every sampled tree was distinct
	// Note: this is not the absolute maximum entropy possible - for that N would equal the number
	// of possible unrooted trees for the number of taxa under consideration. Here, N is equal to
	// the number of sampled trees.
	//
	// maxEntropy = Sum_{i=1}^N p_i \ln(1/p_i)
	//            = Sum_{i=1}^N (1/N) \ln(N)
	//            = \ln(N)
	//
	maxEntropy = log(N);

	// Observed entropy uses distribution now stored
	//
	// obsEntropy = Sum_{i=1}^N p_i \ln(1/p_i)
	//            = Sum_{i=1}^N (n_i/N) \ln(N/n_i)
	//            = (1/N) Sum_{i=1}^N n_i [\ln(N) - \ln(n_i)]
	//            = \ln(N) (1/N) Sum_{i=1}^N n_i - (1/N) Sum_{i=1}^N n_i \ln(n_i)
	//            = \ln(N) - (1/N) Sum_{i=1}^N n_i \ln(n_i)
	//            = maxEntropy - (sumFrqLnFrq / N)
	// 
	obsEntropy = maxEntropy - (sumFrqLnFrq / N);

	// Kullback-Liebler distance measures how far obsEntropy is from maxEntropy
	// In the notation of Burnham and Anderson, we are considering the observed
	// distribution to be the "true" model and the distribution under maximum entropy
	// as the "approximating" model. Thus, the distance we are computing is the 
	// KL distance from the expected distribution under maximum entropy to the observed
	// distribution.
	//
	// KLDistance = Sum_{i=1}^N p_i \ln(\frac{p_i}{1/N})
	//            = Sum_{i=1}^N (n_i/N) [\ln(n_i/N) - \ln(1/N)]
	//            = [ (1/N) Sum_{i=1}^N n_i \ln(n_i/N) ] - [ (1/N) Sum_{i=1}^N n_i \ln(1/N) ]
	//            = [ (1/N) Sum_{i=1}^N n_i (\ln(n_i) - \ln(N) ] - [ (1/N) Sum_{i=1}^N n_i (0 - \ln(N)) ]
	//            = (1/N) Sum_{i=1}^N n_i \ln(n_i) - (1/N) Sum_{i=1}^N n_i \ln(N) + (1/N) Sum_{i=1}^N n_i \ln(N))
	//            = (sumFrqLnFrq / N) - \ln(N) + \ln(N))
	//            = sumFrqLnFrq / N
	//            = maxEntropy - obsEntropy
	//
	KLDistance = maxEntropy - obsEntropy;

	return true;
	}
