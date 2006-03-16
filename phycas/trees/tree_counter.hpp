#ifndef PHO_TREECOUNTER_H
#define PHO_TREECOUNTER_H

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
class TreeCounter
	{
	public:
							TreeCounter();
							~TreeCounter();

		double				GetRootedCount(unsigned n, unsigned m);
		double				GetUnrootedCount(unsigned n, unsigned m);

		double				GetSaturatedRootedCount(unsigned n);
		double				GetSaturatedUnrootedCount(unsigned n);

		double				GetTotalRooted(unsigned n);
		double				GetTotalUnrooted(unsigned n);

	private:

		void				RecalculateCounts(unsigned n);
		void				RecalcNext(unsigned n);
		void				RecalcPrev(unsigned n);

		unsigned			ntax;		/* TreeCounter is currently holding counts for this number of taxa */
		std::vector<double> counts;		/* vector of length n-1 holding counts for trees having all possible numbers of internal nodes */
	};

inline TreeCounter::TreeCounter() 
	{
	ntax = 0;
	}

inline TreeCounter::~TreeCounter() 
	{
	counts.clear();
	}

#endif
