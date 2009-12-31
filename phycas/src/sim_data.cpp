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

//#include "phycas/force_include.h"
#include <cmath>
#include <algorithm>
#include <functional>
#include <boost/format.hpp>
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/cond_likelihood.hpp"
#include "phycas/src/cond_likelihood_storage.hpp"
#include "phycas/src/tip_data.hpp"
#include "phycas/src/sim_data.hpp"
#include "phycas/src/phycas_string.hpp"

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
|	Returns the vector of counts (`binv') that is created in the calctBinned function. The table below shows the 
|   identity of each bin for the case of 4 states. The length of `binv' is 2k - 1, where k is the number of states.
|   The number of states can thus be obtained as (sz + 1)/2, where sz is the length of `binv'.
|>
|	Count   Bin description                 
|   ----------------------------------------
|	n_0     all patterns containing only A  
|	n_1     all patterns containing only C  
|	n_2     all patterns containing only G  
|	n_3     all patterns containing only T  
|	n_4     patterns containing any 2 states
|	n_5     patterns containing any 3 states
|	n_6     patterns containing any 4 states
|   ----------------------------------------
|>
*/
std::vector<double> SimData::getBinnedCounts()
	{
    if (binv.empty())
        {
        buildBinVector(4);  //POL: not good to assume 4 states here - SimData object should know number of states
        }
    return binv;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Builds up the `binv' vector by classifying stored site patterns into the following categories (bins):
|>
|	Count   Bin description                 
|   ----------------------------------------
|	n_0     all patterns containing only A  
|	n_1     all patterns containing only C  
|	n_2     all patterns containing only G  
|	n_3     all patterns containing only T  
|	n_4     patterns containing any 2 states
|	n_5     patterns containing any 3 states
|	n_6     patterns containing any 4 states
|   ----------------------------------------
|>
|	Warning: this function currently assumes no missing data! A missing state is treated as if it were one of the
|	other states in the pattern. For example, the count for the pattern 'AAACA?GAA' would be stuffed into the 3-state
|	bin, with the implicit assumption being that the ? equals either A, C or G.
*/
void SimData::buildBinVector(
  unsigned nstates)   /**< is the number of states */
	{
#if POLPY_OLDWAY	// simulations not yet working with partitioning
	unsigned nbins = 2*nstates - 1;
    binv.clear();
	binv.resize(nbins, 0.0);
	std::set<int8_t> state_set;
	for (SimPatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
		{
		//POL for speed, should classify patterns as they are added to SimData in the first place!
		const VecStateList & p = it->first;
		//std::cerr << "\npattern is";
		state_set.clear();
		unsigned last_state = UINT_MAX;
		for (VecStateList::const_iterator pit = p.begin(); pit != p.end(); ++pit)
			{
			int8_t curr_state = *pit;
			int cs = (int)curr_state;
			if (cs >= 0 && cs < (int)nstates)
				{
				// state is not a gap (-1), missing (nstates), or an ambiguity code (> nstates), so add to set
				//std::cerr << " " << (int)curr_state;
				state_set.insert(curr_state);
				last_state = cs;
				}
			}
		unsigned sz = (unsigned)state_set.size();
		//std::cerr << "\n  num. states in pattern = " << sz;
		PHYCAS_ASSERT(sz > 0);
		PHYCAS_ASSERT(sz <= nstates);

		double this_count = (double)(it->second);
		if (sz == 1)
			{
			// pattern had only one state, so add pattern count to appropriate constant site bin
			binv[last_state] += this_count;
			}
		else
			{
			// pattern had sz states, so add pattern count to appropriate variable site bin
			binv[nstates + sz - 2] += this_count;
			}
		}
#endif
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates the binned t function used in the Gelfand-Ghosh measure for the patterns currently stored in this 
|	object. The binned t function looks like this for multinomial data, 4 DNA states and 4 taxa:
|>
|	Count   Bin description                     t calculation
|   ----------------------------------------------------------------------------------------------
|	n_0     all patterns containing only A      t_0 = (n_0 + eps)/(n + 1) log[(n_0 + eps)/(n + 1)]
|	n_1     all patterns containing only C      t_1 = (n_1 + eps)/(n + 1) log[(n_1 + eps)/(n + 1)]
|	n_2     all patterns containing only G      t_2 = (n_2 + eps)/(n + 1) log[(n_2 + eps)/(n + 1)]
|	n_3     all patterns containing only T      t_3 = (n_3 + eps)/(n + 1) log[(n_3 + eps)/(n + 1)]
|	n_4     patterns containing any 2 states    t_4 = (n_4 + eps)/(n + 1) log[(n_4 + eps)/(n + 1)]
|	n_5     patterns containing any 3 states    t_5 = (n_5 + eps)/(n + 1) log[(n_5 + eps)/(n + 1)]
|	n_6     patterns containing any 4 states    t_6 = (n_6 + eps)/(n + 1) log[(n_6 + eps)/(n + 1)]
|   ----------------------------------------------------------------------------------------------
|   n = n_0 + n_1 + ... + n_6                   t = t_0 + t_1 + ... + t_6
|>
|	where n_i is the count of the number of characters placed into bin i, n is the total number of characters, and 
|	eps equals 1/7 (the inverse of the total number of bins). The total number of bins is (2*nstates - 1), regardless
|	of the number of taxa, which allows this approach to Gelfand-Ghosh to scale nicely to large problems. The drawback,
|	of course, is that information about frequencies of individual rare patterns is not used.
|
|	Warning: this function currently assumes no missing data! A missing state is treated as if it were one of the
|	other states in the pattern. For example, the count for the pattern 'AAACA?GAA' would be stuffed into the 3-state
|	bin, with the implicit assumption being that the ? equals either A, C or G.
*/
double SimData::calctBinned(unsigned nstates)
	{
	// m is the number of distinct patterns
	//double m					= (double)sim_pattern_map.size();

	// s is the number of character states
	//double s					= (double)nstates;
	//double log_s				= std::log(s);

	// n is the number of characters (i.e. sum of counts of all distinct patterns)
	double n					= (double)total_count;
	double n_plus_one			= n + 1.0;
	double log_n_plus_one		= std::log(n_plus_one);

	// epsilon is the number of bins (nstates constant bins plus nstates - 1 additional bins for patterns with 2, 3, ..., nstates states)
	unsigned nbins				= 2*nstates - 1;
	double epsilon				= 1.0/(double)nbins;

	// classify patterns and build up bin, the vector of bin counts
    buildBinVector(nstates);

    //std::ofstream f("bins.txt", std::ios::out | std::ios::app);
	//f.setf(std::ios::floatfield, std::ios::fixed);
	//f.setf(std::ios::showpoint);
    //f << "\nStarting calculation of t:\n";
    //f << "\n--------------------------\n" << std::endl;
    //f << "n                  = " << n << "\n";
    //f << "n_plus_one         = " << n_plus_one << "\n";
    //f << "log_n_plus_one     = " << log_n_plus_one << "\n";
    //f << "nbins              = " << nbins << "\n";
    //f << "epsilon            = " << epsilon << "\n";
    //f << std::endl;

	// now accumulate t by summing over bins
	double t = 0.0;
	for (std::vector<double>::iterator vit = binv.begin(); vit != binv.end(); ++vit)
		{
		double count				= (*vit);
		double count_plus_epsilon	= count + epsilon;
		double first				= count_plus_epsilon/n_plus_one;
		double second				= std::log(count_plus_epsilon) - log_n_plus_one;
		double this_term			= first*second;
		t += this_term;

        //f << "count              = " << count << "\n";
        //f << "count_plus_epsilon = " << count_plus_epsilon << "\n";
        //f << "first              = " << first << "\n";
        //f << "second             = " << second << "\n";
        //f << "this_term          = " << this_term << "\n";
        //f << "t                  = " << t << "\n";
        //f << std::endl;
		}
	//f.close();

	return t;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates the t function used in the Gelfand-Ghosh measure for the patterns currently stored in this object. The
|	t function looks like this for multinomial data and 4 taxa:
|>
|	sum_{i=1}^{256} (x_i + eps)/(n + 1) log[(x_i + eps)/(n + 1)]
|>
|	where 256 is the total number of possible data patterns (=nstates^ntaxa), x_i is the count of the number of 
|	characters having pattern i, n is the total number of characters, and eps equals 1/256 (the inverse of the total 
|	number of possible patterns).
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
	for (SimPatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
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

		for (SimPatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
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
	std::cerr << "taxon_names size = " << taxon_names.size() << std::endl;
	std::cerr << "pattern_length   = " << pattern_length << std::endl;
	std::cerr << "taxon_names: ";
	std::copy(taxon_names.begin(), taxon_names.end(), std::ostream_iterator<std::string>(std::cerr, "|"));
	
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

#if 1
        // Spit out characters in the order in which they were simulated. While this is a nice feature, 
        // it currently requires storing the data twice (patternVect and sim_pattern_map)
        unsigned nchar = (unsigned)patternVect.size();
        for (unsigned k = 0; k < nchar; ++k)
            {
            unsigned j = (unsigned)patternVect[k][i];
            char s = state_symbols[j][0];
            outf << s;
            }
#else
		for (SimPatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
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
#endif
		outf << "\n";
		}

	outf << "  ;" << "\n";
	outf << "end;" << "\n";

	outf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds `tmp_pattern' to `sim_pattern_map' then passes `missing_state' to the wipePattern() function to fill 
|	`tmp_pattern' with invalid values. 
|	
|	Important! Unlike insertPattern, this function does not update total_count to reflect the new sum of pattern counts
|	over all patterns. The calling routine is responsible for updating total_count.
*/
void SimData::insertPatternToRunningAverage(
#if POLPY_NEWWAY
  pattern_count_t count,	/**< is the number of times this pattern was seen */
  pattern_count_t p)		/**< is the inverse of the number of datasets added */
#else // old way
  PatternCountType count,	/**< is the number of times this pattern was seen */
  PatternCountType p)		/**< is the inverse of the number of datasets added */
#endif
	{
	// In debug build, check to make sure there are no elements in tmp_pattern that still equal
	// missing_state (if assert trips, it means that not all elements of tmp_pattern have been 
	// replaced with actual states)
	PHYCAS_ASSERT(std::find(tmp_pattern.begin(), tmp_pattern.end(), missing_state) == tmp_pattern.end());

	// insert the pattern stored in tmp_pattern into sim_pattern_map
	// Add tmp_pattern to sim_pattern_map if it has not yet been seen, otherwise increment the count 
	// for this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
	SimPatternMapType::iterator lowb = sim_pattern_map.lower_bound(tmp_pattern);
	if (lowb != sim_pattern_map.end() && !(sim_pattern_map.key_comp()(tmp_pattern, lowb->first)))
		{
		// Pattern is already in sim_pattern_map, so just modify its count
#if POLPY_NEWWAY
		pattern_count_t curr = lowb->second;
		pattern_count_t one_minus_p = (pattern_count_t)(1.0 - p);
		pattern_count_t new_value = curr*one_minus_p + count*p;
#else // old way
		PatternCountType curr = lowb->second;
		PatternCountType one_minus_p = (PatternCountType)(1.0 - p);
		PatternCountType new_value = curr*one_minus_p + count*p;
#endif
		lowb->second = new_value;

		std::ofstream f("smuteye2.txt", std::ios::out | std::ios::app);
		f << curr << '\t' << one_minus_p << '\t' << count << '\t' << p << '\t' << new_value << '\t' << "found" << std::endl;
		f.close();
		}
	else
		{
		// tmp_pattern has not yet been stored in sim_pattern_map
		sim_pattern_map.insert(lowb, SimPatternMapType::value_type(tmp_pattern, count*p));

		std::ofstream f("smuteye2.txt", std::ios::out | std::ios::app);
		f << 0.0 << '\t' << (1-p) << '\t' << count << '\t' << p << '\t' << (count*p) << '\t' << "new" << std::endl;
		f.close();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds `tmp_pattern' to `sim_pattern_map' then passes `missing_state' to the wipePattern() function to fill 
|	`tmp_pattern' with invalid values. Finally, increments total_count to reflect the new total number of counts 
|	over all patterns.
*/
#if POLPY_NEWWAY
void SimData::insertPattern(
  pattern_count_t count)	/**< is the count to be associated with the pattern now stored in `tmp_pattern' */
#else // old way
void SimData::insertPattern(PatternCountType count)
#endif
	{
	// In debug build, check to make sure there are no elements in tmp_pattern that still equal
	// missing_state (if assert trips, it means that not all elements of tmp_pattern have been 
	// replaced with actual states)
	PHYCAS_ASSERT(std::find(tmp_pattern.begin(), tmp_pattern.end(), missing_state) == tmp_pattern.end());

	// insert the pattern stored in tmp_pattern into sim_pattern_map
	// Add tmp_pattern to sim_pattern_map if it has not yet been seen, otherwise increment the count 
	// for this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
	SimPatternMapType::iterator lowb = sim_pattern_map.lower_bound(tmp_pattern);
	if (lowb != sim_pattern_map.end() && !(sim_pattern_map.key_comp()(tmp_pattern, lowb->first)))
		{
		// Pattern is already in sim_pattern_map, so just modify its count
		lowb->second += count;
		}
	else
		{
		// tmp_pattern has not yet been stored in sim_pattern_map
		sim_pattern_map.insert(lowb, SimPatternMapType::value_type(tmp_pattern, count));
		}

    patternVect.push_back(tmp_pattern);
	total_count += count;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds data currently stored in `sim_pattern_map' to the patterns already in `other'. The value of `mult' is used to
|	modify the counts before they are added to `other'; that is, the count of each pattern added to `other' is the 
|	original count multiplied by `mult'. Normally, `mult' would be specified to be 1.0, but in some cases it is
|	necessary to build up SimData objects that represent averages of other SimData objects, and it is these situations
|	where `mult' is handy. Assumes that `mult' is positive and non-zero. Assumes that `pattern_length' for this SimData
|	object is identical to the `pattern_length' of `other'.
*/
#if POLPY_NEWWAY
void SimData::addDataTo(
  SimData & other, 			/**< is the SimData object that will receive the data currently contained in this SimData object */
  pattern_count_t mult)		/**< is the factor multiplied by each pattern's count before pattern is stored in `other' */
#else // old way
void SimData::addDataTo(SimData & other, PatternCountType mult)
#endif
	{
#if POLPY_OLDWAY	// simulations not yet working with partitioning
	PHYCAS_ASSERT(mult > 0.0);

	// If this object has no patterns, return immediately
	if (total_count == 0.0)
		return;

	// If other is empty, then it most likely needs to be initialized
	if (other.getTotalCount() == 0)
		{
		other.resetPatternLength(pattern_length);
		}
	PHYCAS_ASSERT(pattern_length == other.getPatternLength());

	for (SimPatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
		{
		PatternCountType count = it->second;
		
		// GetCurrPattern returns a workspace for building up a pattern to be added to other
		VecStateList & other_pattern = other.getCurrPattern();

		// Copy pattern represented by *it to other's workspace
		std::copy(it->first.begin(), it->first.end(), other_pattern.begin());

		// Add the pattern in other's workspace to other's pattern map
		PatternCountType mult_count = mult*count;
		other.insertPattern(mult_count);
		}
#endif
	}

#if 0
/*----------------------------------------------------------------------------------------------------------------------
|	Adds data currently stored in `sim_pattern_map' to the patterns already in `other'. The value of `mult' is used to
|	modify the counts before they are added to `other'; that is, the count of each pattern added to `other' is the 
|	original count multiplied by `mult'. Normally, `mult' would be specified to be 1.0, but in some cases it is
|	necessary to build up SimData objects that represent averages of other SimData objects, and it is these situations
|	where `mult' is handy. Assumes that `mult' is positive and non-zero. Assumes that `pattern_length' for this SimData
|	object is identical to the `pattern_length' of `other'. This function maintains a running average. It depends on 
|	'other' keeping track of how many times it has received new data (the 'num_additions' data member is used for this).
|	
|	Here is an example of how the running average works. Suppose there are only two possible patterns, and the following 
|	pairs of pattern counts are generated by performing three posterior predictive simulations:
|>
|	[x1, y1] are the counts for the two patterns in posterior predictive simulated dataset number 1
|	[x2, y2] are the counts for the two patterns in posterior predictive simulated dataset number 2
|	[x3, y3] are the counts for the two patterns in posterior predictive simulated dataset number 3
|>
|	In each case, the sum of x and y is n, so total_count is 3n after all three have been added. The strategy before 
|	would be to add the x values and y values, then divide by 3 at the end, so the average dataset would be:
|>
|	[(x1+x2+x3)/3,  (y1+y2+y3)/3]
|>
|	Waiting to do the division until the end can lead to overflow (see the first entry in the BUGS file for details), 
|	however, so instead this function maintains a running average as follows:
|>
|	After adding pair 1:  p = 1/1, 1-p = 0
|	                      x = x*(1-p) + x1*p = x1
|	                      y = y*(1-p) + y1*p = y1
|	                      total_count = total_count*(1-p) + (x1 + y1)*p = x1 + y1
|	After adding pair 2:  p = 1/2, 1-p = 1/2
|	                      x = x*(1-p) + x2*p = (x1 + x2)/2
|	                      y = y*(1-p) + y2*p = (y1 + y2)/2 
|	                      total_count = total_count*(1-p) + (x2 + y2)*p = (x1 + x2 + y1 + y2)/2
|	After adding pair 3:  p = 1/3, 1-p = 2/3
|	                      x = x*(1-p) + x3*p = (x1 + x2 + x3)/3
|	                      y = y*(1-p) + y3*p = (y1 + y2 + y3)/3 
|	                      total_count = total_count*(1-p) + (x3 + y3)*p = (x1 + x2 + x3 + y1 + y2 + y3)/3
|>
*/
void SimData::addToRunningAverage(
  SimData & other, 
  PatternCountType mult)
	{
	PHYCAS_ASSERT(mult > 0.0);
	PHYCAS_ASSERT(total_count > 0.0);

	// If other is empty, then it most likely needs to be initialized
	if (other.getTotalCount() == 0)
		{
		other.resetPatternLength(pattern_length);
		}
	PHYCAS_ASSERT(pattern_length == other.getPatternLength());

	// calculate p, the weight to be used for this addition
	unsigned nadd = other.getNumAdditions();
	PatternCountType numer = (PatternCountType)1.0;
	PatternCountType denom = (PatternCountType)(1 + nadd);
	PatternCountType p = numer/denom;

	PatternCountType sum = 0.0;
	for (SimPatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
		{
		PatternCountType count = it->second;
		
		// GetCurrPattern returns a workspace for building up a pattern to be added to other
		VecStateList & other_pattern = other.getCurrPattern();

		// Copy pattern represented by *it to other's workspace
		std::copy(it->first.begin(), it->first.end(), other_pattern.begin());

		// Add the pattern now in other's tmp_pattern workspace to other's pattern map
		PatternCountType mult_count = mult*count;
		sum += mult_count;
		other.insertPatternToRunningAverage(mult_count, p);
		}

	// Update total_count
	PatternCountType curr_total = other.getTotalCount();
	PatternCountType one_minus_p = (PatternCountType)(1.0 - p);
	PatternCountType new_total = curr_total*one_minus_p + sum*p;
	other.setTotalCount(new_total);
	other.setNumAdditions(nadd + 1);

	if (nadd > 1)
		{
		std::exit(0);
		}
	}
#endif // #if 0

}	// namespace phycas
