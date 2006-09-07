//#include "phycas/force_include.h"
#include <cmath>
#include <algorithm>
#include <functional>
#include <boost/format.hpp>
#include "phypy/src/cipres/ConfigDependentHeaders.h"	// int8_t typedef
#include "phypy/src/basic_tree.hpp"
#include "phypy/src/tip_data.hpp"
#include "phypy/src/sim_data.hpp"
#include "phypy/src/pyphy_string.hpp"

const int8_t phycas::SimData::missing_state = -1;	// GCC 3.2.3 (Red Hat Linux 3.2.3-20) requires initializing here rather than in header
using std::exp;

class LookupStateSymbol : public std::unary_function<int8_t, char>
	{
	public:
		LookupStateSymbol(const std::vector<std::string> & state_symbols) : symbols(state_symbols) {}

		char operator()(int8_t i)
			{
			PHYCAS_ASSERT((unsigned)i < (unsigned)symbols.size());
			return symbols[i][0];
			}

	private:
		const std::vector<std::string> & symbols;
	};

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates the t function used in the Gelfand-Ghosh measure for the patterns currently stored in this object. The
|	t function looks like this for multinomial data and 4 taxa:
|>
|	sum_{i=1}^{256} (x_i + eps)/(n + 1) log[(x_i + eps)/(n + 1)]
|>
|	where 256 is the total number of possible data patterns (=nstates^ntaxa), n is the number of characters, and eps 
|	equals 1/256 (the inverse of the total number of possible patterns).
*/
double SimData::calct(unsigned nstates)
	{
	// m is the number of distinct patterns
	double m					= (double)sim_pattern_map.size();
	//unused:  double log_m				= std::log(m);

	// s is the number of character states
	double s					= (double)nstates;
	double log_s				= std::log(s);

	// n is the number of characters (i.e. sum of counts of all distinct patterns)
	double n					= (double)total_count;
	double n_plus_one			= n + 1.0;
	double log_n_plus_one		= std::log(n_plus_one);

	// ntaxa is the number of taxa in the tree
	double ntaxa				= (double)pattern_length;
	double ntaxa_times_log_s	= ntaxa*log_s;

	// epsilon is the inverse of the total number of possible patterns
	double epsilon				= std::exp(-ntaxa_times_log_s);

	// added_on is the sum of all terms in t in which the pattern count equals zero
	double log_term				= -ntaxa_times_log_s - log_n_plus_one;
	double added_on				= (1.0 - m*epsilon)*log_term/n_plus_one;

	// now compute the sum of all terms in t in which the pattern count is greater than zero
	double sum					= 0.0;
	for (PatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
		{
		double count				= (double)it->second;
		double count_plus_epsilon	= count + epsilon;
		double first				= count_plus_epsilon/n_plus_one;
		double second				= std::log(count_plus_epsilon) - log_n_plus_one;
		double this_term			= first*second;
		sum += this_term;
		}

	double t = sum + added_on;

	return t;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns simulated data stored in `sim_pattern_map' as a string in the form of a two-column table. The first column
|	is labeled "Count" and the second is labeled "Pattern". The Count column shows the number of times its associated 
|	pattern was inserted using the insertPattern function. The Pattern column shows a representation of the pattern 
|	itself, using symbols for states provided in the `state_symbols' argument. The `state_symbols' argument should be 
|	a vector of single-character strings supplying a symbol to represent each state that might show up in any pattern.
|	Assumes that no state in any pattern stored in `sim_pattern_map' is greater than or equal to the length of the 
|	`state_symbols' vector (because states are used as indices into `state_symbols').
*/
std::string SimData::patternTable(
  const StringVect & state_symbols)		/**< is a vector of strings representing states (e.g. {"A", "C", "G", "T"}). Note that each state symbol should be a string of length 1 (i.e. a single character) */	
	{
	PHYCAS_ASSERT(state_symbols.size() > 0);

	outstr.clear();
	if (sim_pattern_map.empty())
		{
		outstr = "Sorry, no patterns are stored";
		}
	else
		{
		outstr = "     Count  Pattern";

		for (PatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
			{
			// Output the count first
			outstr << str(boost::format("\n%10.1f") % it->second) << "  ";

			// Now output the pattern
			std::transform(it->first.begin(), it->first.end(), std::back_inserter(outstr), LookupStateSymbol(state_symbols));
			}
		}

	return outstr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Saves simulated data stored in `sim_pattern_map' to a file named `filename'. If the file already exists, it will be
|	overwritten without warning. If the file does not yet exist, it will be created. The file written will be a valid
|	NEXUS data file suitable for executing in phylogenetic analysis software that reads the NEXUS file format. The
|	supplied `taxon_names' will be used in the matrix command of the NEXUS file to specify the names of the taxa. 
|	Assumes that the number of elements in `taxon_names' equals `pattern_length'. The `datatype' argument should be the
|	correct NEXUS datatype (e.g. "dna", "standard") for the data simulated. The symbols used for the states are 
|	supplied in the `state_symbols' vector. Each element of this vector should be a single-character string. Assumes
|	that no state in any pattern stored in `sim_pattern_map' is greater than or equal to the length of the 
|	`state_symbols' vector (because states are used as indices into `state_symbols').
*/
void SimData::saveToNexusFile(
  const std::string filename,			/**< is the name of the file to create containing the simulated data stored in this object */
  const StringVect & taxon_names,		/**< is a vector containing the taxon names to use in the saved file */
  const std::string datatype,			/**< is a string to be used as the NEXUS datatype (e.g. "dna" or "standard") */	
  const StringVect & state_symbols)		/**< is a vector of strings representing states (e.g. {"A", "C", "G", "T"}). Note that each state symbol should be a string of length 1 (i.e. a single character) */	
	{
	PHYCAS_ASSERT(state_symbols.size() > 0);
	PHYCAS_ASSERT(taxon_names.size() == pattern_length);

	// Find length of longest string in taxon_names vector; this is used later for formatting purposes
	// The 2 is included in case apostrophes are needed when taxon names are output in the matrix command
	unsigned length_of_longest_name = 2 + (unsigned)std::max_element(taxon_names.begin(), taxon_names.end(), StringLengthLess())->length();
	
	std::ofstream outf(filename.c_str());

	outf << "#nexus" << "\n\n";
	outf << "begin data;" << "\n";
	outf << str(boost::format("  dimensions ntax=%d nchar=%d;") % pattern_length % (unsigned)total_count) << "\n";
	outf << "  format datatype=" << datatype << ";" << "\n";
	outf << "  matrix" << "\n";

	// Create a format string to use with boost::format that left-justifies (the "-" flag) the 
	// taxon names in a field of width length_of_longest_name 
	std::string fmtstr = str(boost::format("    %%-%ds") % length_of_longest_name);

	for (unsigned i = 0; i < pattern_length; ++i)
		{
		if (taxon_names[i].find(' ') != std::string::npos)
			{
			std::string s = "'";
			s += taxon_names[i];
			s += "'";
			outf << str(boost::format(fmtstr) % s) << "  ";
			}
		else
			{
			outf << str(boost::format(fmtstr) % taxon_names[i]) << "  ";
			}

		for (PatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
			{
			// The first member of it is the state, which must be converted from its coded form
			// to the standard symbol for this data type (e.g. A, C, G or T for DNA characters)
			int8_t j = (*it).first[i];

			PHYCAS_ASSERT(j < (int8_t)state_symbols.size());
			char s = state_symbols[j][0];	// use the first (and hopefully only) character in the string at position j

			// The second member of it is the pattern count
			unsigned n = (unsigned)it->second;	//@POL assuming counts not fractional
			for (unsigned k = 0; k < n; ++k)
				{
				outf << s;
				}
			}
		outf << "\n";
		}

	outf << "  ;" << "\n";
	outf << "end;" << "\n";

	outf.close();
	}

}	// namespace phycas
