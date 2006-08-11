// splitmgr.h  implements class used to manage split objects
// Copyright (C) 2000 Paul O. Lewis.
//

#ifndef PHO_SPLITMGR_H
#define PHO_SPLITMGR_H

#include <boost/shared_ptr.hpp>
#include "phycas/trees/split.hpp"
#include "ncl/taxa/nxs_taxa_listener.hpp"
#include "ncl/misc/string_extensions.hpp"
#include "ncl/output/nxs_output_stream.hpp"

class NxsTable;
class SplitManager;
class TreeNode;
class SplitInfo
	{
	public:
		SplitInfo();
		SplitInfo(bool triv, unsigned nTimesSeen, double sumEdges);
		
		void Clear();
		const bool IsTrivial() const 
			{
			return isTrivial;
			}
		const double GetLength() const 
			{
			return sum_lengths;
			}
		const UInt GetFrequency() const 
			{
			return frequency;
			}
		const double GetPosteriorProbability(const UInt totalNSteps) const
			{
			return (double)GetFrequency()/(double) totalNSteps;
			}
		const double GetMeanEdgeLength() const
			{
			return (GetFrequency() > 0 ? (GetLength()/(double)GetFrequency()) : 0.0);
			}
	private:
		double 		sum_lengths;
		unsigned 	frequency;
		bool 		isTrivial;
		friend class SplitManager;
	};
	
inline void SplitInfo::Clear()
	{
	sum_lengths = 0.0;
	frequency = 0;
	isTrivial = false;
	}
	
inline SplitInfo::SplitInfo()
	{
	Clear();
	}

inline SplitInfo::SplitInfo(
  bool triv,
  unsigned nTimesSeen, 
  double sumEdges)
  	:sum_lengths(sumEdges), 
  	frequency(nTimesSeen),
  	isTrivial(triv)
	{
	}
	

typedef std::map< Split, SplitInfo > SplitMap;

class SplitManager : public NxsTaxaListener
	{
	public:
		SplitManager(NxsTaxaManager *);
		~SplitManager();
		
		double 		AddSplit(const Split& s, Length edgeLen, bool tip_node, bool increment_total = true );
		bool 		CheckFull();
		bool 		IsFull() const;
		unsigned 	GetNumSplits() const;
		double 		GetStdDev() const;
		UInt		GetNumTreesRecorded() const
			{
			return total;
			}
		unsigned 	SplitFrequency(const Split& s) const;
		double		SplitProportion(const Split& s) const;
		double 		SplitMeanLength(const Split& s) const;
		void 		FlushSplitsMap();
		void 		IncrementTotal();
		void 		SetSplitAutoInc(bool autoInc, unsigned incr);
		void 		SetSplitLimit( unsigned limit );
		
		template<typename OutStream> void ShowSplitsMap( OutStream& out ) const; //@should be template or typedef
		template<typename OutStream> unsigned	OutputSplitsMap(OutStream &, double lower_bound, bool show_trivial) const; //@should be template or typedef
		void 		SaveSplitsAsNexusFile(std::ostream & splitsf) const;
		//unsigned	AddSplitsMapToTable(NxsTable* t, bool showtrivials, double cutoff, unsigned width = 12, unsigned precision = 6) const;
		
		void 		RecordAllSplits(Tree &tree);
		void 		SetBrlensFromSplits(Tree& tree) const;
		void 		BuildMajorityRuleTreeID(TreeID& id, double cutoff = 0.5) const;
		bool 		ComputeEntropy(double &obsEntropy);
		int	 		EnumSplitsByComplexity(std::vector<int> &v);
		
		const SplitMap &GetSplitsMap() const
			{
			return smap;
			}
		bool RecordSplit(TreeNode *nd);
		
		void TaxaChanged(BaseTaxaManager *, TaxaChangeType);
		
		void BuildMajorityRuleTree(Tree & tree, double cutoffProportion, unsigned ntax) const;
		SplitMap::const_iterator begin() const
			{
			return smap.begin();
			}
		SplitMap::const_iterator end() const
			{
			return smap.end();
			}
		
	private:
		unsigned 		split_limit;
		unsigned 		incrSplitLimitBy;
		bool 			split_limit_exceeded;
		bool 			autoSplitInc;
		unsigned 		total;
		SplitMap 		smap;	
		NxsTaxaManager	*taxaMgr; //POL 22jun2005
	};

//typedef boost::shared_ptr<SplitManager> SplitManagerShPtr;

inline void SplitManager::SetSplitAutoInc(bool autoInc, unsigned incr)
	{
	autoSplitInc = autoInc;
	incrSplitLimitBy = incr;
	}

inline void SplitManager::SetSplitLimit( unsigned limit )
	{
	split_limit = limit;
	CheckFull();
	}

inline unsigned SplitManager::GetNumSplits() const
	{
	return (unsigned) smap.size();
	}

inline void SplitManager::IncrementTotal()
	{
	++total;
	CheckFull();
	}

inline bool SplitManager::IsFull() const
	{
	return split_limit_exceeded;
	}


template<typename OutStream> void SplitManager::ShowSplitsMap( OutStream& out ) const
	{
	out << "\nSplit frequencies:\n";

	long numSplits = (long)smap.size();
	if( numSplits > 0 ) 
		{
		std::string s;
		for(SplitMap::const_iterator i = smap.begin(); i != smap.end(); ++i ) 
			{
			out << i->second.frequency << "  ";
			s.clear();
			i->first.CreateAndAppendPatternRepresentation( &s );
			out << s;
			s.clear();
			s << i->first;
			out << s << '\n';
			}
		}
	else
		out << "  No splits have been stored\n";
	out << "Total splits   = " << total << '\n';
	out << "Total distinct = " << numSplits << '\n';
	out << "Std. dev. = " << GetStdDev() << '\n';
	}

template<typename OutStream> unsigned SplitManager::OutputSplitsMap(OutStream &out, double lower_bound, bool show_trivial) const
	{
	//@@@ need to do this with a table
	if( smap.size() < 1 )
		return 0;
		
	assert( total > 0L );
	unsigned nt = Split::splitNTax;
	unsigned repColWid = nt;
	if (repColWid <= 6)
		repColWid = 7;
	bool repTooWide = false;
	char repshow[21];
	if (repColWid >= 20)
		{
		repColWid = 20;
		repTooWide = true;
		}
	
	assert(nt < 1024);
	std::string message;
	RightJustifyString(message, "Split", 12);	//@POL Mark, how do I access the column width set by the user?
	RightJustifyString(message, "Pattern", repColWid + 2);
	RightJustifyString(message, "Length", 12);
	RightJustifyString(message, "Prior", 12);
	RightJustifyString(message, "Freq.", 12);
	RightJustifyString(message, "Prob.", 12);
	out << message << ncl::endl;

	long k = 1L;
	for (SplitMap::const_iterator i = smap.begin(); i != smap.end(); ++i, ++k) 
		{
		const SplitInfo &si = (*i).second;
		unsigned frq = si.frequency;

		//@POL needs work
		//  It is possible for a non-trivial split to have posterior probability 1.0
		//	Also, whether or not to display trivial splits should be up to the user
		//
		if (!show_trivial && si.isTrivial)
			continue;
		double enumSplitsCutoff = lower_bound; //@POL replace this enumSplitsCutoff with user setting
		double prob = (double)frq / (double)total;
		if (prob <= enumSplitsCutoff)
			continue;

		// Split number
		//
		message.clear();
		RightJustifyLong(message, k, 12);	//@POL replace 12 with column width

		// Split representation
		//
		std::string s;
		i->first.CreateAndAppendPatternRepresentation(&s);
		if (repTooWide)
			{
			unsigned m = 0;
			const char *const begStr = s.c_str();
			const char *p = begStr;
			while (*p != '\0')
				{
				repshow[m++] = *p++;
				if (k == repColWid || *p == '\0')
					{
					if (*p == '\0')
						{
						while (m < repColWid)
							repshow[m++] = ' ';	
						}
					repshow[m] = '\0';
					k = 0;
					if ((unsigned)(p - begStr) > repColWid)
						message << "\n      ";
					RightJustifyString(message, repshow, repColWid + 2);
					}
				}
			}
		else
			RightJustifyString(message, s, repColWid + 2);

		// Mean length of split
		//
		double mean_length = (frq > 0L ? (si.sum_lengths / (double)frq ) : 0.0);
		RightJustifyDbl(message, mean_length, 12, 5);	//@POL replace 12 with column width and 5 with precision

		// Prior probability of split
		//
		double priorprob = (*i).first.CalcPriorProb();
		RightJustifyDbl(message, priorprob, 12, 5);

		// Frequency of split
		//
		RightJustifyDbl(message, frq, 12, 5);

		// Posterior probability of split
		//
		RightJustifyDbl(message, prob, 12, 5);

		out << message << ncl::endl;

		//@POL TODO: should allow user to cancel this if 
		// the list of splits turns out to be very long
		}
		
	return (unsigned)smap.size();
	}

#endif
