#if ! defined(SIM_DATA_HPP)
#define SIM_DATA_HPP

#include "pyphy/common/states_patterns.hpp"

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
		void						resetPatternLength(unsigned ntaxa);
		void						wipePattern();
		void						setState(unsigned pos, int8_t state);
		void						insertPattern(PatternCountType count);

		double						calct(unsigned nstates);
		void						includeDataFrom(SimData &);
		unsigned					getPatternLength();
		VecStateList &				getCurrPattern();
		PatternCountType			getTotalCount();
		void						addDataTo(SimData & other, PatternCountType mult);
		void						divideBy(PatternCountType total);

		std::string					patternTable(const StringVect & state_symbols);
		void						saveToNexusFile(const std::string filename, const StringVect & taxon_names, const std::string datatype, const StringVect & state_symbols);

		const static int8_t			missing_state;			/**< The value that represents missing data */

		const PatternMapType &		getSimPatternMap() const;

	private:

		PatternCountType			total_count;			/**< The number of patterns inserted since sim_pattern_map was last cleared (note that this is not sim_pattern_map.size(), but instead equals the sum of all the counts) */
		unsigned					pattern_length;			/**< Number of taxa, used to reserve elements in new pattern vectors */
		PatternMapType				sim_pattern_map;		/**< Keys are patterns, values are pattern counts */
		VecStateList				tmp_pattern;			/**< Workspace for building up a pattern */
		std::string					outstr;					/**< Workspace for building up a tabular representation of `sim_pattern_map' (used by showPatternMap function) */
	};

typedef boost::shared_ptr<SimData>	SimDataShPtr;

} // namespace phycas

#include "sim_data.inl"

#endif

