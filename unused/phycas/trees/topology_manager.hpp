// topomgr.h  implements class used to manage map of tree topologies
// Copyright (C) 2000 Paul O. Lewis.
//

#ifndef __TOPOMGR_H
#define __TOPOMGR_H

#include "ncl/taxa/nxs_taxa_listener.hpp"
#include <boost/shared_ptr.hpp>
class Tree;
class TreeID;
class NxsTaxaManager;
class NxsTable;
class TopoInfo
	{
	friend class TopologyManager;

		double		weight;
		unsigned	frequency;
			
	public:
		TopoInfo() : weight(1.0), frequency(0){}
		TopoInfo(double w, unsigned f) : weight(w), frequency(f){}
	};

typedef std::map< TreeID, TopoInfo > TopoMap;
typedef std::multimap< unsigned, unsigned > TopoSortMultimap;
typedef std::vector<unsigned> VecTopoSort;

class PhoFloor;

class SplitManager;
typedef boost::shared_ptr<SplitManager>	SplitManagerShPtr;

class TopologyManager: public NxsTaxaListener
	{
		unsigned					topo_limit;
		unsigned					incrTopoLimitBy;
		bool 						topo_limit_exceeded;
		bool 						autoTopoInc;
		unsigned					nTaxa;				/* current number of active taxa */
		unsigned 					totalAddAttempts;	/* #times a topology was added as workspace to determine order sorted by frequency */
		TopoMap 					tmap;
		SplitManagerShPtr			splitMgr;			/* split manager to use when assigning edge lengths to tree topologies */
		mutable bool				sortVectDirty;	/* used as workspace to determine order sorted by frequency */
		mutable VecTopoSort			sortvect;	/* used as workspace to determine order sorted by frequency */
		NxsTaxaManager				*taxaMgr; //POL 22jun2005
		
	public:
									TopologyManager(NxsTaxaManager *);
									~TopologyManager();
		
		unsigned 					GetFrequency(unsigned which) const ;
		unsigned 					GetNumTopologiesStored() const;
		TreeID 						GetID( unsigned which ) const;
		TreeID 						GetIDofBest() const;
		std::pair<TreeID, TopoInfo>		GetTopology(unsigned i) const;

		unsigned 					GetSorted(unsigned i) const;
		bool 						IsFull() const ;
		double	 					TopologyWeight(const TreeID& t) const;
		unsigned 					TopologyFrequency(const TreeID& t) const;
		
		void						AddTopology(const TreeID& t, double wt = -1.0);
		void 						AddTopology(const Tree & tree);
		bool 						CheckFull();
		void						SetSplitMgr(SplitManagerShPtr sm);
		void 						SetTopoAutoInc(bool autoInc, unsigned incr);
		void 						SetTopologyLimit( unsigned limit );
		void 						FlushTopoMap();
		
		bool 						ComputeEntropy(double &obsEntropy, double &maxEntropy, double &KLDistance) const;
		
		void 						SaveTopoMap(std::ostream& out, int width = 12) const;
		template<class T> unsigned	ShowTopoMap(T &outS, bool show_wts = false, unsigned width = 12, unsigned precision = 6) const;
		unsigned					ShowTopoMapAsTable(NxsTable & t, bool show_wts = false, unsigned width = 12, unsigned precision = 6) const;
			
		void						TaxaChanged(BaseTaxaManager *taxaMgr, TaxaChangeType);

	private:
		void						CreateSortVect() const;
	};

//typedef boost::shared_ptr<TopologyManager> TopologyManagerShPtr;


inline void TopologyManager::SetSplitMgr(SplitManagerShPtr sm)
	{
	splitMgr = sm;
	}

inline unsigned TopologyManager::GetNumTopologiesStored() const
	{
	return (unsigned)tmap.size();
	}

inline bool TopologyManager::IsFull() const
	{
	return topo_limit_exceeded;
	}

#define SAVE_TOPOLGY_TREEFILE
#include "ncl/output/nxs_table_cell.hpp"
#include "phycas/trees/tree_manip.hpp"
#include "phycas/trees/tree_id.hpp"
#include "phycas/trees/split_manager.hpp"
#include "phycas/trees/tree.hpp"
#include "ncl/misc/string_extensions.hpp"

template<class T> unsigned TopologyManager::ShowTopoMap(T &outStream, bool show_wts /* = false */, unsigned width /* = 12 */, unsigned precision /* = 6 */ ) const
	{
	unsigned numTopologies = GetNumTopologiesStored();
	if( numTopologies < 1L )
		return 0L;
		
	NxsTableCell::def_width     = width;
	NxsTableCell::def_precision = precision;

	std::string message;
	RightJustifyString(message, "Topology", 12);
	RightJustifyString(message, "Freq.", 12);
	if (show_wts)
		RightJustifyString(message, "Weight", 12);
	RightJustifyString(message, "Prob.", 12);
#	if defined(ADHOC_FOR_ENTROPY_PAPER)
		message << "  Newick";
#	endif
	outStream << message << ncl::endl;
	
#if defined(SAVE_TOPOLGY_TREEFILE)
	std::ofstream treef("topologies.tre");
	treef << "#nexus\n\n";
	treef << "begin trees;\n";
#endif

	long k = 0L;
	//TopoMap::const_iterator i = tmap.begin();
	//for(; i != tmap.end(); ++i)
	for (unsigned j = 0; j < numTopologies; ++j)
		{
		unsigned topologyIndex = GetSorted(j);
		std::pair<TreeID, TopoInfo> i = GetTopology(topologyIndex);
		long frq = (long)(i.second.frequency);
		double prob = (double)frq / (double)totalAddAttempts;
		double enumTopolCutoff = 0.0; //@POL replace this enumTopolCutoff with user setting
		if (prob < enumTopolCutoff)
			continue;
		message.clear();
		RightJustifyLong(message, ++k, 12);
		RightJustifyLong(message, frq, 12);
		if (show_wts)
			RightJustifyDbl(message, i.second.weight, 12, 5); //@ casting double to long
		RightJustifyDbl(message, prob, 12, 3);
#		if defined(SAVE_TOPOLGY_TREEFILE)
			Tree tree;
			TreeID treeId = i.first;
			TreeManip tmanip(&tree);
			tmanip.SimpleBuildTreeFromID(nTaxa, treeId);
			if (splitMgr)
				splitMgr->SetBrlensFromSplits(tree);
			std::string newick;
			tree.AppendNewickRepresentation(newick, true, true);
			treef << "  utree " << frq << " = " << newick << ';' << std::endl;
#		endif
#		if defined(ADHOC_FOR_ENTROPY_PAPER)
			std::string newick;
			Tree tree;
			TreeID treeId = i.first;
			tree.BuildTreeFromID(treeId);
			floor->splitMgr.SetBrlensFromSplits(tree);
			tree.MakeNewick(newick, false);
			message << "  " << newick;
#		endif
	outStream << message << ncl::endl;

		//@POL TODO: should inform user if no topologies were frequent
		// enough to be displayed
		//
		//@POL TODO: should allow user to cancel this if 
		// the list of trees turns out to be very long
		}
		
#	if defined(SAVE_TOPOLGY_TREEFILE)
		treef << "end;" << std::endl;
		treef.close();
#	endif

	return numTopologies;
	}

#endif
