#ifndef PHYCAS_PHYLODIVERSITY_SUMMARY_HPP
#define PHYCAS_PHYLODIVERSITY_SUMMARY_HPP

#define NEW_PD_WAY


class PhylodiversitySummary
	{
	public:
		bool	isValid; // false if first taxon is in the selected set
		unsigned treeNumber;
		double cumTotal;
		double cumI;
#		if defined(NEW_PD_WAY)
			double cumE;
#		else
			double cumEmin;
			double cumEmax;
#		endif
		double avgSelTipLen;
		double avgUnselTipLen;
		double pEI, pIT, pET, medSelTipLen, medUnselTipLen;
		PhylodiversitySummary()
			:isValid(true),
			cumTotal(0.0),
			cumI(0.0),
#			if defined(NEW_PD_WAY)
				cumE(0.0)
#			else
				cumEmin(0.0),
				cumEmax(0.0)
#			endif
			{}
		void calculateStatisticsFromCumulatives();
	};

#endif
