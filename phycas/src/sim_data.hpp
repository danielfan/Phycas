/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#if ! defined(SIM_DATA_HPP)
#define SIM_DATA_HPP

#include "phycas/src/states_patterns.hpp"

typedef	PatternMapType	SimPatternMapType;

namespace phycas
{

class Tree;
typedef std::vector<std::string>				StringVect;

/*----------------------------------------------------------------------------------------------------------------------
|	Serves as a container for data simulated by the TreeLikelihood::simulate function. Has member functions for saving 
|	the stored data to a nexus file.
*/
class SimData
	{
	friend class TreeLikelihood;

	public:

									SimData();


		void						clear();
		void						zeroCounts();
		void						appendCountsToFile(std::string fn, bool binary);
		void						debugAppendCountsToFile(std::string row_name, std::string fn);
		std::string					createMapleTuples(unsigned row, unsigned cutoff);
		std::vector<std::string>	getPatterns(std::vector<std::string> symbols);
		void						resetPatternLength(unsigned ntaxa);
		void						wipePattern();
		void						setState(unsigned pos, int8_t state);
		void						insertPattern(PatternCountType count);

        void                        buildBinVector(unsigned nstates);
		std::vector<double>		    getBinnedCounts();

        double						calct(unsigned nstates);
		double						calctBinned(unsigned nstates);
		void						includeDataFrom(SimData &);
		unsigned					getPatternLength();
		VecStateList &				getCurrPattern();
		PatternCountType			getTotalCount();
		unsigned					getNUniquePatterns();
		void						addDataTo(SimData & other, PatternCountType mult);
		void						setNumAdditions(unsigned n);
		unsigned					getNumAdditions();
		void						setTotalCount(PatternCountType total);
		void						multBy(PatternCountType factor);
		void						divideBy(PatternCountType factor);

		std::string					patternTable(const StringVect & state_symbols);
		void						saveToNexusFile(const std::string filename, const StringVect & taxon_names, const std::string datatype, const StringVect & state_symbols);

		const static int8_t			missing_state;			/**< The value that represents missing data */

		const SimPatternMapType &	getSimPatternMap() const;

	private:

		void						insertPatternToRunningAverage(PatternCountType count, PatternCountType p);

	private:

		PatternCountType			total_count;			/**< The number of patterns inserted since sim_pattern_map was last cleared (note that this is not sim_pattern_map.size(), but instead equals the sum of all the counts) */
		unsigned					pattern_length;			/**< Number of taxa, used to reserve elements in new pattern vectors */
		SimPatternMapType			sim_pattern_map;		/**< Keys are patterns, values are pattern counts */
		VecStateList				tmp_pattern;			/**< Workspace for building up a pattern */
		std::string					outstr;					/**< Workspace for building up a tabular representation of `sim_pattern_map' (used by showPatternMap function) */
        std::vector<double>         binv;                   /**< Stores binned counts if calctBinned is called; otherwise, will be an empty vector. */
	};

typedef boost::shared_ptr<SimData>	SimDataShPtr;

} // namespace phycas

#include "sim_data.inl"

#endif

