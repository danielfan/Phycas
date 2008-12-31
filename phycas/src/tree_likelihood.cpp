/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
|																			  |
|  This program is free software; you can redistribute it and/or modify		  |
|  it under the terms of the GNU General Public License as published by		  |
|  the Free Software Foundation; either version 2 of the License, or		  |
|  (at your option) any later version.										  |
|																			  |
|  This program is distributed in the hope that it will be useful,			  |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of			  |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the			  |
|  GNU General Public License for more details.								  |
|																			  |
|  You should have received a copy of the GNU General Public License along	  |
|  with this program; if not, write to the Free Software Foundation, Inc.,	  |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.				  |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//#include "phycas/force_include.h"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/tip_data.hpp"
#include "phycas/src/internal_data.hpp"
#include "phycas/src/sim_data.hpp"
#include "phycas/src/phycas_string.hpp"
#include "phycas/src/basic_lot.hpp"
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <numeric>
#include "phycas/src/edge_iterators.hpp"
#include "phycas/src/univents.hpp"

// formerly in tree_likelihood.inl
#include "phycas/src/edge_endpoints.hpp"

static int8_t codon_state_codes[] =
	{
	0,	// 0 AAA
	1,	// 1 AAC
	2,	// 2 AAG
	3,	// 3 AAT
	4,	// 4 ACA
	5,	// 5 ACC
	6,	// 6 ACG
	7,	// 7 ACT
	8,	// 8 AGA
	9,	// 9 AGC
	10, // 10 AGG
	11, // 11 AGT
	12, // 12 ATA
	13, // 13 ATC
	14, // 14 ATG
	15, // 15 ATT
	16, // 16 CAA
	17, // 17 CAC
	18, // 18 CAG
	19, // 19 CAT
	20, // 20 CCA
	21, // 21 CCC
	22, // 22 CCG
	23, // 23 CCT
	24, // 24 CGA
	25, // 25 CGC
	26, // 26 CGG
	27, // 27 CGT
	28, // 28 CTA
	29, // 29 CTC
	30, // 30 CTG
	31, // 31 CTT
	32, // 32 GAA
	33, // 33 GAC
	34, // 34 GAG
	35, // 35 GAT
	36, // 36 GCA
	37, // 37 GCC
	38, // 38 GCG
	39, // 39 GCT
	40, // 40 GGA
	41, // 41 GGC
	42, // 42 GGG
	43, // 43 GGT
	44, // 44 GTA
	45, // 45 GTC
	46, // 46 GTG
	47, // 47 GTT
	61, // 48 TAA stop
	48, // 49 TAC
	61, // 50 TAG stop
	49, // 51 TAT
	50, // 52 TCA
	51, // 53 TCC
	52, // 54 TCG
	53, // 55 TCT
	61, // 56 TGA stop
	54, // 57 TGC
	55, // 58 TGG
	56, // 59 TGT
	57, // 60 TTA
	58, // 61 TTC
	59, // 62 TTG
	60	// 63 TTT
	};

namespace phycas
{

Univents & getUniventsRef(TreeNode &nd)
	{
	return (nd.IsTip() ? nd.GetTipData()->getUniventsRef() : nd.GetInternalData()->getUniventsRef());
	}

const Univents & getUniventsConstRef(const TreeNode &nd) 
	{
	return (nd.IsTip() ? nd.GetTipData()->getUniventsConstRef() : nd.GetInternalData()->getUniventsConstRef());
	}

// **************************************************************************************
// ***** Former TreeLikelihood inlines (begin) ******************************************
// **************************************************************************************

/*----------------------------------------------------------------------------------------------------------------------
|	TreeLikelihood constructor.
*/
TreeLikelihood::TreeLikelihood(
  ModelShPtr mod)		/**< is the substitution model */
  :
  likelihood_root(0),
  store_site_likes(false),
  no_data(false),
  nTaxa(0),
  num_patterns(0),
  num_states(mod->getNStates()),
  num_rates(mod->getNRatesTotal()),
  model(mod), 
  rate_means(mod->getNRatesTotal(), 1.0), //POL_BOOKMARK rate_means vector init
  rate_probs(mod->getNRatesTotal(), 1.0), 
  nevals(0),
  debugging_now(false)
  ,using_unimap(false),
  univentProbMgr(mod),
  sMat(NULL),
  sMatValid(false)
	{
#if POLPY_NEWWAY	//CLAShPtr
	cla_pool = CondLikelihoodStorageShPtr(new CondLikelihoodStorage());
#endif
	mod->recalcRatesAndProbs(rate_means, rate_probs);
	underflow_policy.setTriggerSensitivity(50);
	underflow_policy.setCorrectToValue(10000.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	TreeLikelihood destructor.
*/
TreeLikelihood::~TreeLikelihood()
	{
	std::cerr << "\n>>>>> TreeLikelihood dying..." << std::endl;
	if (sMat != NULL)
		DeleteTwoDArray<unsigned>(sMat);
	} 

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the vector of site likelihoods (data member `site_likelihood') computed in 
|   TreeLikelihood::harvestLnLFromValidEdge.
*/
const std::vector<double> & TreeLikelihood::getSiteLikelihoods() const
    {
    return site_likelihood;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the vector of pattern counts (data member `pattern_counts').
*/
const std::vector<double> & TreeLikelihood::getPatternCounts() const
    {
    return pattern_counts;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the data member `charIndexToPatternIndex'.
*/
const std::vector<unsigned> & TreeLikelihood::getCharIndexToPatternIndex() const
    {
    return charIndexToPatternIndex;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `store_site_likes'.
*/
bool TreeLikelihood::storingSiteLikelihoods() const
    {
    return store_site_likes;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the current value of the data member `store_site_likes' to true if `yes' is true, and false if `yes' is false.
*/
void TreeLikelihood::storeSiteLikelihoods(
  bool yes)    /**< is the value to which `store_site_likes' should be set */
    {
    store_site_likes = yes;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Swaps InternalData data member `state_time' and edge lengths for the supplied nodes `nd1' and `nd2'. Assumes 
|	`nd1' and `nd2' are both internal nodes.
*/
void TreeLikelihood::swapInternalDataAndEdgeLen(
  TreeNode * nd1,	/**< is the first of two nodes whose InternalData structures and edge lengths are to be swapped */
  TreeNode * nd2)	/**< is the second of two nodes whose InternalData structures and edge lengths are to be swapped */
	{
	PHYCAS_ASSERT(nd1->IsInternal());
	PHYCAS_ASSERT(nd1->GetInternalData() != NULL);
	PHYCAS_ASSERT(nd2->IsInternal());
	PHYCAS_ASSERT(nd2->GetInternalData() != NULL);

	// Swap edge lengths
	double edgelen = nd1->GetEdgeLen();
	nd1->SetEdgeLen(nd2->GetEdgeLen());
	nd2->SetEdgeLen(edgelen);

	if (using_unimap)
		nd1->GetInternalData()->swapUnivents(nd2->GetInternalData());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `using_unimap'.
*/
bool TreeLikelihood::isUsingUnimap()
	{
	return using_unimap;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Specifies whether the TreeLikelihood object will use uniformized mapping likelihoods or the standard Felsenstein-
|	style integrated likelihood.
*/
void TreeLikelihood::useUnimap(
  bool yes_or_no)	/**< is either true (to assume uniformized mapping) or false (to use the integrated likelihood) */
	{
	using_unimap = yes_or_no;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string representation of the current value of `sMat' for debugging purposes.
*/
std::string TreeLikelihood::debugShowSMatrix()
	{
	unsigned i, j, v, total, trace;
	total = 0;
	trace = 0;
	std::vector<unsigned> rowsum(num_states, 0);
	std::vector<unsigned> colsum(num_states, 0);
	std::string s = str(boost::format(" %8s") % " ");
	for (j = 0; j < num_states; ++j)
		{
		s += boost::str(boost::format(" %8d") % j);
		}
	s += boost::str(boost::format(" %8s\n") % std::string("sum"));
	for (i = 0; i < num_states; ++i)
		{
		s += boost::str(boost::format(" %8d") % i);
		for (j = 0; j < num_states; ++j)
			{
			v = sMat[i][j];
			s += boost::str(boost::format(" %8d") % v);
			total += v;
			rowsum[i] += v;
			colsum[j] += v;
			if (i == j)
				trace += v;
			}
		s += boost::str(boost::format(" %8d\n") % rowsum[i]);
		}

	// print out row of column sums, then total at the end
	s += boost::str(boost::format(" %8s") % std::string("sum"));
	for (j = 0; j < num_states; ++j)
		{
		s += boost::str(boost::format(" %8d") % colsum[j]);
		}
	s += boost::str(boost::format(" %8d\n") % total);
	s += boost::str(boost::format("trace	 = %d\n") % trace);
	s += boost::str(boost::format("nunivents = %d\n") % nunivents);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recalculates the two-dimensional matrix `sMat' of numbers of all 16 possible univent transitions over all sites.
|	A no-op if using_unimap is false.
*/
void TreeLikelihood::recalcSMatrix(
  TreeShPtr t)	/**< is the tree to use */
	{
	if (isUsingUnimap())
		{
		// Make sure sMat exists
		if (sMat == NULL)
			sMat = NewTwoDArray<unsigned>(num_states, num_states); 

		// Zero every element of sMat
		for (unsigned i = 0; i < num_states; ++i)
			{
			for (unsigned j = 0; j < num_states; ++j)
				sMat[i][j] = 0;
			}

		// Loop over nodes in the tree
		preorder_iterator nd = t->begin();	// skip the tip root node
		for (++nd; nd != t->end(); ++nd)
			{
			const Univents & u = getUniventsConstRef(*nd);
			PHYCAS_ASSERT(u.isValid());
			const std::vector<StateMapping> & v = u.getVecEventsVecConstRef();
			const std::vector<int8_t> & states_vec =  u.getEndStatesVecConstRef();
			std::vector<int8_t>::const_iterator statesIt = states_vec.begin();

			for (std::vector<StateMapping>::const_iterator sit = v.begin(); sit != v.end(); ++sit, ++statesIt)
				{
				PHYCAS_ASSERT(statesIt != states_vec.end());
				const StateMapping & stlist = (*sit);
				if (!stlist.empty())
					{
					int8_t prev_state = *statesIt;
					const StateMapping::const_iterator endIt = stlist.end();
					for (StateMapping::const_iterator it = stlist.begin(); it != endIt; ++it)
						{
						const int8_t new_state = *it;
						sMat[prev_state][new_state] += 1;
						prev_state = new_state;
						}
					}
				}
			}
		} // if using unimap
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The node `slider' is moved along the combined edges of `slider' and `other'. If `fraction' is positive, `slider' 
|	gains edge length at the expense of `other'. If `fraction' is negative, `other' gains edge length at the expense of
|	`slider'. This function not only changes edge lengths, but also takes care of transferring any univents that lie 
|	along the affected path. If `slider' is a sibling of other, then `slider' is somewhat of a misnomer since it does 
|	not actually slide; instead, the parent of both `slider' and `other' does the sliding.
*/
void TreeLikelihood::slideNode(
  double fraction, 
  TreeNode * slider, 
  TreeNode * other)
	{
	bool other_is_child = (other->GetParent() == slider);
	PHYCAS_ASSERT(other_is_child || slider->GetParent() == other->GetParent());
	if (other_is_child)
		{
		if (fraction < 0.0)
			{
			// Fraction is negative, which means that other gains edge length at its base at the 
			// expense of slider, which loses edge length at its end
			//						   |
			//						   |
			//	-----------------------*--------------------------------*
			//				 <---------^slider							^other
			//	<----xnew---><------------------ynew-------------------->
			//	<-----------x----------><--------------y---------------->
			//				 <--delta-->
			//
			//PHYCAS_ASSERT(which_case == 2 || which_case == 5 || which_case == 11);
			double x	  = slider->GetEdgeLen();
			double y	  = other->GetEdgeLen();
			double delta  = -x*fraction;
			double xnew	  = x - delta;
			double ynew	  = y + delta;
			slider->SetEdgeLen(xnew);
			other->SetEdgeLen(ynew);
			if (isUsingUnimap())
				{
				PHYCAS_ASSERT(false);
#				if 0 // this code has not been revisited since the June, 2008 Durham meeting
				double cutoff = 1.0 - fraction;
				StateTimeListVect & slider_vect = (slider->IsInternal() ? slider->GetInternalData()->state_time : slider->GetTipData()->state_time);
				StateTimeListVect & other_vect	= (other->IsInternal()	? other->GetInternalData()->state_time : other->GetTipData()->state_time);
				PHYCAS_ASSERT(num_patterns == slider_vect.size());
				PHYCAS_ASSERT(num_patterns == other_vect.size());

				// Iterate over sites
				StateTimeListVect::iterator other_site_it  = other_vect.begin();
				StateTimeListVect::iterator slider_site_it = slider_vect.begin();
				for (; slider_site_it != slider_vect.end(); ++slider_site_it, ++other_site_it)
					{
					// Grab state_time vectors for current site from both slider and other
					StateTimeList & other_stlist  = (*other_site_it);
					StateTimeList & slider_stlist = (*slider_site_it);

					// Iterate over slider's state_time vector starting from the end, copying
					// elements to the front of other as needed
					unsigned nmoved = 0;
					StateTimeList::reverse_iterator rit, slider_first, slider_last;
					slider_first = slider_stlist.rbegin() + 1;
					slider_last	 = slider_first;
					for (rit = slider_first; rit != slider_stlist.rend() - 1; ++rit)
						{
						StateTimePair & p = (*rit);
						float time = p.second;
						if (time > cutoff)
							{
							other_stlist.insert(other_stlist.begin() + 1, p);
							++nmoved;
							++slider_last;
							}
						else
							break;
						}

					// If any elements were copied to other, erase them now from slider
					if (nmoved > 0)
						{
						// Gymnastics needed because erase does not work with reverse_iterators
						StateTimeList::iterator first = (slider_last+1).base();
						StateTimeList::iterator last  = (slider_first+1).base();
						slider_stlist.erase(++first, ++last);
						}

					// Recalculate the times in slider
					StateTimeList::iterator it;
					for (it = slider_stlist.begin() + 1; it != slider_stlist.end() - 1; ++it)
						{
						StateTimePair & p = (*it);
						double time = p.second;
						double new_time = time*x/xnew;
						p.second = (float)new_time;
						}

					// Recalculate the newly-added times at beginning of other
					StateTimeList::iterator one_beyond_newly_added = other_stlist.begin() + nmoved + 1;
					for (it = other_stlist.begin() + 1; it != one_beyond_newly_added; ++it)
						{
						StateTimePair & p = (*it);
						double time = p.second;
						double new_time = (time*x - xnew)/ynew;
						p.second = (float)new_time;
						}

					// Recalculate the original times in other
					for (it = one_beyond_newly_added; it != other_stlist.end() - 1; ++it)
						{
						StateTimePair & p = (*it);
						double time = p.second;
						double new_time = (delta + time*y)/ynew;
						p.second = (float)new_time;
						}
					}	// loop over sites
#				endif //0 June 2008, Durham
				} // if unimap
			}	// fraction negative
		else
			{
			// Fraction is positive, so slider gains edge length at its end at the expense of other, 
			// which loses edge length at its base
			//						   |
			//						   |
			//	-----------------------*--------------------------------*
			//					 slider^-------------->					^other
			//	<------------------xnew---------------><------ynew------>
			//	<----------x----------><--------------y----------------->
			//						   <-----delta---->
			//
			//PHYCAS_ASSERT(which_case == 1 || which_case == 4 || which_case == 10);
			double x	  = slider->GetEdgeLen();
			double y	  = other->GetEdgeLen();
			double delta  = y*fraction;
			double xnew	  = x + delta;
			double ynew	  = y - delta;
			slider->SetEdgeLen(xnew);
			other->SetEdgeLen(ynew);
			if (isUsingUnimap())
				{
				PHYCAS_ASSERT(false);
#				if 0 // this code has not been revisited since the June, 2008 Durham meeting
				double cutoff = fraction;
				StateTimeListVect & slider_vect = (slider->IsInternal() ? slider->GetInternalData()->state_time : slider->GetTipData()->state_time);
				StateTimeListVect & other_vect	= (other->IsInternal() ? other->GetInternalData()->state_time : other->GetTipData()->state_time);
				PHYCAS_ASSERT(num_patterns == slider_vect.size());
				PHYCAS_ASSERT(num_patterns == other_vect.size());

				// Iterate over sites
				StateTimeListVect::iterator other_site_it  = other_vect.begin();
				StateTimeListVect::iterator slider_site_it = slider_vect.begin();
				for (; slider_site_it != slider_vect.end(); ++slider_site_it, ++other_site_it)
					{
					// Grab state_time vectors for current site from both slider and other
					StateTimeList & other_stlist  = (*other_site_it);
					StateTimeList & slider_stlist = (*slider_site_it);

					// Iterate over other's state_time vector starting from the beginning, copying
					// elements to the end of slider as needed
					unsigned nmoved = 0;
					StateTimeList::iterator it, other_first, other_last;
					other_first = other_stlist.begin() + 1;
					other_last = other_first;
					for (it = other_first; it != other_stlist.end() - 1; ++it)
						{
						StateTimePair & p = (*it);
						float time = p.second;
						if (time < cutoff)
							{
							slider_stlist.insert(slider_stlist.end() - 1, p);
							++nmoved;
							++other_last;
							}
						else
							break;
						}

					// If any elements were copied to slider, erase them now from other
					if (nmoved > 0)
						{
						other_stlist.erase(other_first, other_last);
						}

					// Recalculate the original times in slider
					unsigned offset = (unsigned)slider_stlist.size() - nmoved - 1;
					StateTimeList::iterator first_newly_added = slider_stlist.begin() + offset;
					for (it = slider_stlist.begin() + 1; it != first_newly_added; ++it)
						{
						StateTimePair & p = (*it);
						double time = p.second;
						double new_time = time*x/xnew;
						p.second = (float)new_time;
						}

					// Recalculate the new times added to the end of slider
					for (it = first_newly_added; it != slider_stlist.end() - 1; ++it)
						{
						StateTimePair & p = (*it);
						double time = p.second;
						double new_time = (x + time*y)/xnew;
						p.second = (float)new_time;
						}

					// Recalculate the times in other
					for (it = other_stlist.begin() + 1; it != other_stlist.end() - 1; ++it)
						{
						StateTimePair & p = (*it);
						double time = p.second;
						double new_time = (time*y - delta)/ynew;
						p.second = (float)new_time;
						}
					}	// loop over sites
#				endif // 0 // this code has not been revisited since the June, 2008 Durham meeting
				} // if unimap
			}	// fraction positive
		}	// other is child
	else
		{
		// other is the sibling, not the child, of slider
		if (fraction < 0.0)
			{
			// Fraction is negative, which means that other gains edge length at its base at the 
			// expense of slider, which loses edge length at its base
			//
			//	   +-----------x-----------* slider
			//	---+
			//	   +----------------y-------------* other
			//
			//				<-----delta---->
			//
			//	   +--xnew--* slider
			//	---+
			//	   +----------------------ynew----------------------* other
			//				 
			//PHYCAS_ASSERT(which_case == 8);
			double x	  = slider->GetEdgeLen();
			double y	  = other->GetEdgeLen();
			double delta  = -x*fraction;
			double xnew	  = x - delta;
			double ynew	  = y + delta;
			slider->SetEdgeLen(xnew);
			other->SetEdgeLen(ynew);
			if (isUsingUnimap())
				{
				PHYCAS_ASSERT(false);
#				if 0 // this code has not been revisited since the June, 2008 Durham meeting
				double cutoff = fraction;
				StateTimeListVect & slider_vect = (slider->IsInternal() ? slider->GetInternalData()->state_time : slider->GetTipData()->state_time);
				StateTimeListVect & other_vect	= (other->IsInternal()	? other->GetInternalData()->state_time : other->GetTipData()->state_time);
				PHYCAS_ASSERT(num_patterns == slider_vect.size());
				PHYCAS_ASSERT(num_patterns == other_vect.size());

				// Iterate over sites
				StateTimeListVect::iterator other_site_it  = other_vect.begin();
				StateTimeListVect::iterator slider_site_it = slider_vect.begin();
				for (; slider_site_it != slider_vect.end(); ++slider_site_it, ++other_site_it)
					{
					// Grab state_time vectors for current site from both slider and other
					StateTimeList & other_stlist  = (*other_site_it);
					StateTimeList & slider_stlist = (*slider_site_it);

					// Iterate over slider's state_time vector starting from the beginning, copying
					// elements to the front of other as long as the time is less than cutoff
					unsigned nmoved = 0;
					StateTimeList::iterator it, slider_first, slider_last;
					slider_first = slider_stlist.begin() + 1;
					slider_last	 = slider_first;
					for (it = slider_first; it != slider_stlist.end() - 1; ++it)
						{
						StateTimePair & p = (*it);
						float time = p.second;
						if (time < cutoff)
							{
							other_stlist.insert(other_stlist.begin() + 1, p);
							++nmoved;
							++slider_last;
							}
						else
							break;
						}

					// If any elements were copied to other, erase them now from slider
					if (nmoved > 0)
						{
						slider_stlist.erase(slider_first, slider_last);
						}

					// Recalculate the times in slider
					for (it = slider_stlist.begin() + 1; it != slider_stlist.end() - 1; ++it)
						{
						StateTimePair & p = (*it);
						double time = p.second;
						double new_time = (time*x - delta)/xnew;
						p.second = (float)new_time;
						}

					// Recalculate the newly-added times at base of other
					StateTimeList::iterator one_beyond_newly_added = other_stlist.begin() + nmoved + 1;
					for (it = other_stlist.begin() + 1; it != one_beyond_newly_added; ++it)
						{
						StateTimePair & p = (*it);
						double time = p.second;
						double new_time = (delta - time*x)/ynew;
						p.second = (float)new_time;
						}

					// Recalculate the original times in other
					for (it = one_beyond_newly_added; it != other_stlist.end() - 1; ++it)
						{
						StateTimePair & p = (*it);
						double time = p.second;
						double new_time = (delta + time*y)/ynew;
						p.second = (float)new_time;
						}
					}	// loop over sites
#				endif //  0 // this code has not been revisited since the June, 2008 Durham meeting
				} // if unimap
			}	// fraction negative
		else
			{
			// Fraction is positive, so slider gains edge length at its base at the expense of other, 
			// which loses edge length at its base
			//
			//	   +-----------x-----------* slider
			//	---+
			//	   +----------------y-------------* other
			//
			//							   <-----delta----->
			//
			//	   +------------------xnew------------------* slider
			//	---+
			//	   +---------ynew------* other
			//				 
			//PHYCAS_ASSERT(which_case == 7);
			double x	  = slider->GetEdgeLen();
			double y	  = other->GetEdgeLen();
			double delta  = y*fraction;
			double xnew	  = x + delta;
			double ynew	  = y - delta;
			slider->SetEdgeLen(xnew);
			other->SetEdgeLen(ynew);
			if (isUsingUnimap())
				{
				PHYCAS_ASSERT(false);
#				if 0 // this code has not been revisited since the June, 2008 Durham meeting
				double cutoff = fraction;
				StateTimeListVect & slider_vect = (slider->IsInternal() ? slider->GetInternalData()->state_time : slider->GetTipData()->state_time);
				StateTimeListVect & other_vect	= (other->IsInternal() ? other->GetInternalData()->state_time : other->GetTipData()->state_time);
				PHYCAS_ASSERT(num_patterns == slider_vect.size());
				PHYCAS_ASSERT(num_patterns == other_vect.size());

				// Iterate over sites
				StateTimeListVect::iterator other_site_it  = other_vect.begin();
				StateTimeListVect::iterator slider_site_it = slider_vect.begin();
				for (; slider_site_it != slider_vect.end(); ++slider_site_it, ++other_site_it)
					{
					// Grab state_time vectors for current site from both slider and other
					StateTimeList & other_stlist  = (*other_site_it);
					StateTimeList & slider_stlist = (*slider_site_it);

					// Iterate over other's state_time vector starting from the beginning, 
					// copying elements to the beginning of slider as needed
					unsigned nmoved = 0;
					StateTimeList::iterator it, other_first, other_last;
					other_first = other_stlist.begin() + 1;
					other_last	= other_first;
					for (it = other_first; it != other_stlist.end(); ++it)
						{
						StateTimePair & p = (*it);
						float time = p.second;
						if (time <= cutoff)
							{
							slider_stlist.insert(slider_stlist.begin() + 1, p);
							++nmoved;
							++other_last;
							}
						else
							break;
						}

					// If any elements were copied to slider, erase them now from other
					if (nmoved > 0)
						{
						other_stlist.erase(other_first, other_last);
						}

					// Recalculate the newly-added times at the base of slider's edge
					StateTimeList::iterator one_beyond_newly_added = slider_stlist.begin() + nmoved + 1;
					for (it = slider_stlist.begin() + 1; it != one_beyond_newly_added; ++it)
						{
						StateTimePair & p = (*it);
						double time = p.second;
						double new_time = (delta - time*y)/xnew;
						p.second = (float)new_time;
						}

					// Recalculate the original times at the end of slider's edge
					for (it = one_beyond_newly_added; it != slider_stlist.end() - 1; ++it)
						{
						StateTimePair & p = (*it);
						double time = p.second;
						double new_time = (delta + time*x)/xnew;
						p.second = (float)new_time;
						}

					// Recalculate the times in other
					for (it = other_stlist.begin() + 1; it != other_stlist.end() - 1; ++it)
						{
						StateTimePair & p = (*it);
						double time = p.second;
						double new_time = (time*y - delta)/ynew;
						p.second = (float)new_time;
						}
					}	// loop over sites
#				endif //  0 // this code has not been revisited since the June, 2008 Durham meeting
				} // if unimap
			}	// fraction positive
		}	// other is sibling
	}	// TreeLikelihood::slideNode

/*----------------------------------------------------------------------------------------------------------------------
|	Refreshes the univent mapping for all sites over the entire tree. This function will wipe out all stored states 
|   and times on the edges of the tree and create a fresh set compatible with the tip states.
*/
void TreeLikelihood::fullRemapping(
  TreeShPtr t,		        /**< is the tree to use for the mapping */
  LotShPtr rng,             /**< is the random number generator to use for the mapping */
  bool doSampleUnivents)    /**< is True if ... */
	{
	t->renumberInternalNodes(t->GetNTips());
	univentProbMgr.recalcUMat();

	// Reset the matrix of uniformized transition counts
	if (sMat == NULL)
		sMat = NewTwoDArray<unsigned>(num_states, num_states); 
	for (unsigned i = 0; i < num_states; ++i)
		{
		for (unsigned j = 0; j < num_states; ++j)
			sMat[i][j] = 0;
		}
	sMatValid = false;
	nunivents = 0;

	TreeNode * root_tip = t->GetFirstPreorder();
	PHYCAS_ASSERT(root_tip->IsTipRoot());
	TreeNode * subroot = root_tip->GetLeftChild();

	//useAsLikelihoodRoot(subroot);
	//invalidateAwayFromNode(*subroot);
	//invalidateBothEnds(subroot);
	useAsLikelihoodRoot(NULL);

	// Work down the tree in postorder fashion updating conditional likelihoods
	// (This code stolen from TreeLikelihood::calcLnLFromNode.)

	// The functor below will return true if the conditional likelihood arrays pointing away from the
	// focal_node are up-to-date, false if they need to be recomputed
	NodeValidityChecker valid_functor = boost::bind(&TreeLikelihood::isValid, this, _1, _2);

	// The iterator below will visit nodes that need their CLAs updated centripetally 
	// (like a postorder traversal but also coming from below the focal node). Each node 
	// visited is guaranteed by valid_functor to need its CLA updated.
	effective_postorder_edge_iterator iter(subroot, valid_functor);
	effective_postorder_edge_iterator iter_end;
	for (; iter != iter_end; ++iter)
		refreshCLA(*iter->first, iter->second);

	// Must now refresh the CLA of the focal node (subroot) because this one was not
	// recalculated by the effective_postorder_edge_iterator
	refreshCLA(*subroot, root_tip);

	// Work up the tree in preorder fashion choosing states for internal nodes along the way
	TreeNode * nd = subroot;
	SquareMatrix p_mat_trans_scratch(num_states, 0.0);
	double ** p_mat_trans_scratch_ptr = p_mat_trans_scratch.GetMatrix();
	for (; nd != NULL; nd = nd->GetNextPreorder())
		{
		Univents & nd_univents = getUniventsRef(*nd);
		const double nd_edge_len = nd->GetEdgeLen();
		if (nd->IsInternal())
			{
			// Choose states for all sites at this internal node
			InternalData * nd_data =  nd->GetInternalData();
			const LikeFltType * cla = nd_data->getChildCondLikePtr()->getCLA();
			ConstPMatrices pmatrices = nd_data->getConstPMatrices();
			double const * const * pmatrix = pmatrices[0]; // index is 0 because assuming only one rate 
			if (nd == subroot)
				{
				//
				// Choose states and mappings for the subroot node
				//
				// The subroot is a special case for several reasons:
				// 1. its conditional likelihood arrays do not take into account its parent, 
				//	  which is the tip serving as the root of the tree. 
				// 2. this is the first node to get assigned a state, hence we must use the 
				//	  equilibrium freqencies as the prior instead of the transition probability
				//	  from the state below
				// 3. because the subroot serves as the likelihood root, we must keep track
				//	  of the frequency of each state assigned to this node in order to compute
				//	  the likelihood (these frequencies are stored in the vector
				//	  TreeLikelihood::obs_state_counts)
				std::vector<LikeFltType> prob(num_states*num_patterns);
				const std::vector<double> & freqs = model->getStateFreqs();

				// Gather arrays needed from the root_tip
				double						   root_tip_edge_len  = subroot->GetEdgeLen();
				const TipData &				   root_tip_data	  = *(root_tip->GetTipData());
				double * * *				   root_tip_p		  = root_tip_data.getMutableTransposedPMatrices();
				calcPMatTranspose(root_tip_p, root_tip_data.getConstStateListPos(),	 root_tip_edge_len);
				const double * const * const * root_tip_tmatrix	  = root_tip_data.getConstTransposedPMatrices();
				const int8_t *				   root_tip_codes	  = root_tip_data.getConstStateCodes();
				unsigned offset = 0;
				for (unsigned j = 0; j < num_patterns; ++j)
					{
					double total = 0.0;
					for (unsigned k = 0; k < num_states; ++k)
						{
						const double conditional_likelihood = *cla++;
						const unsigned root_tip_state = (unsigned)root_tip_codes[j];
						const double transition_prob = root_tip_tmatrix[0][k][root_tip_state];	  // note: first index is 0 because assuming no rate heterogeneity
						const double unnorm_prob = freqs[k]*conditional_likelihood*transition_prob;
						total += unnorm_prob;
						prob[offset+k] = unnorm_prob;
						}
					for (unsigned k = 0; k < num_states; ++k)
						prob[offset+k] /= total; 
					offset += num_states;
					}
				obs_state_counts.resize(num_states);
				univentProbMgr.sampleRootStates(nd_univents, &prob[0], *rng.get(), true, &obs_state_counts[0]);
				if (doSampleUnivents)
					{
					fillTranspose(p_mat_trans_scratch_ptr, root_tip_tmatrix[0], num_states);
					univentProbMgr.sampleUnivents(nd_univents, root_tip_edge_len, root_tip_codes, const_cast<const double * const*>(p_mat_trans_scratch_ptr), *rng.get(), sMat);
					}
				}
			else
				{
				//
				// Choose states and mappings for non-subroot internal node
				//
				TreeNode * par = nd->GetParent();
				PHYCAS_ASSERT(par != NULL);
				PHYCAS_ASSERT(par->IsInternal());
				Univents & ndP_univents = getUniventsRef(*par);
				const std::vector<int8_t> & par_states_vec = ndP_univents.getEndStatesVecRef();
				const int8_t * par_states_ptr = &par_states_vec[0];
				univentProbMgr.sampleDescendantStates(nd_univents, pmatrix, cla, par_states_ptr, *rng.get());
				if (doSampleUnivents)
					univentProbMgr.sampleUnivents(nd_univents, nd_edge_len,  par_states_ptr, pmatrix, *rng.get(), sMat);
				}
			}
		else if (doSampleUnivents)
			{
			// Choose mappings for tip node
			TreeNode * par = nd->GetParent();
			PHYCAS_ASSERT(par != NULL);
			PHYCAS_ASSERT(par->IsInternal());
			Univents & ndP_univents = getUniventsRef(*par);
			const std::vector<int8_t> & par_states_vec = ndP_univents.getEndStatesVecRef();
			const int8_t * par_states_ptr = &par_states_vec[0];
			TipData * nd_data =	 nd->GetTipData();
			ConstPMatrices tmatrices = nd_data->getTransposedPMatrices();
			fillTranspose(p_mat_trans_scratch_ptr, tmatrices[0], num_states);
			univentProbMgr.sampleUnivents(nd_univents, nd_edge_len, par_states_ptr, const_cast<const double **>(p_mat_trans_scratch_ptr), *rng.get(), sMat);
			}
		}
	if (doSampleUnivents)
		sMatValid = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the setTriggerSensitivity function of the data member `underflow_policy' to set the number of edges that must
|	be traversed before taking action to prevent underflow.
*/
void TreeLikelihood::setUFNumEdges(
  unsigned nedges)	/**< is the number of edges to traverse before taking action to prevent underflow */
	{
	underflow_policy.setTriggerSensitivity(nedges);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of bytes allocated for each CLA. This equals sizeof(LikeFltType) times the product of the number of
|	patterns, number of rates and number of states. Calls corresponding function of data member `cla_pool' to get the
|	value returned.
*/
unsigned TreeLikelihood::bytesPerCLA() const
	{
#if POLPY_NEWWAY	//CLAShPtr
	return cla_pool->bytesPerCLA();
#else
	return cla_pool.bytesPerCLA();
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of CondLikelihood objects created since the `cla_pool' data member was constructed, or since the
|	last call to the function clearStack of `cla_pool', which resets the value to zero. Calls corresponding function of
|	data member `cla_pool' to get the value returned.
*/
unsigned TreeLikelihood::numCLAsCreated() const
	{
#if POLPY_NEWWAY	//CLAShPtr
	return cla_pool->numCLAsCreated();
#else
	return cla_pool.numCLAsCreated();
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current number of CondLikelihood objects stored in `cla_pool'. The total number of CLAs currently
|	checked out to the tree can be obtained as TreeLikelihood::numCLAsCreated() minus TreeLikelihood::numCLAsStored().
|	Calls corresponding function of data member `cla_pool' to get the value returned.
*/
unsigned TreeLikelihood::numCLAsStored() const
	{
#if POLPY_NEWWAY	//CLAShPtr
	return cla_pool->numCLAsStored();
#else
	return cla_pool.numCLAsStored();
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of `likelihood_root'. See TreeLikelihood::useAsLikelihoodRoot for more information about
|	the meaning of the likelihood root.
*/
TreeNode * TreeLikelihood::getLikelihoodRoot()
	{
	return likelihood_root;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of the node currently serving as the likelihood root. If `likelihood_root' is NULL, returns -1
|	instead to indicate that no node is currently designated as the likelihood root. This function was written primarily
|	for use by the TreeViewer.py application for debugging purposes.
*/
int TreeLikelihood::getLikelihoodRootNodeNum() const
	{
	return (likelihood_root ? (int)likelihood_root->GetNodeNumber() : -1);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a shared pointer to the CondLikelihood object that would be used to compute the likelihood conditional on a 
|	particular state being assigned to `focal_nd' when `avoid' is the likelihood root. This function obtains the correct
|	shared pointer, but if that pointer does not point to an object, it goes to `cla_pool' to get another one.
*/
CondLikelihoodShPtr getCondLikePtr(
  TreeNode * focal_nd,	/**< is the focal node */
  TreeNode * avoid)		/**< is the focal node neighbor (node closer to likelihood root) */
	{
	EdgeEndpoints e(focal_nd, avoid);
	return getCondLikePtr(e);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a shared pointer to the CondLikelihood object that would be used to compute the likelihood conditional on a 
|	particular state being assigned to `focal_nd' when `avoid' is the likelihood root. 
*/
ConstCondLikelihoodShPtr getValidCondLikePtr(
  const TreeNode * focal_nd,	/**< is the focal node */
  const TreeNode * avoid)		/**< is the focal node neighbor (node closer to likelihood root) */
	{
	ConstEdgeEndpoints e(focal_nd, avoid);
	return getValidCondLikePtr(e);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a shared pointer to the CondLikelihood object that would be used to compute the likelihood conditional on a
|	particular state being assigned to the focal node of `edge' when `edge' focal neighbor is the likelihood root. This
|	function obtain the correct shared pointer, but if that pointer does not point to a CondLikelihood object, it goes 
|	to `cla_pool' to get one.
*/
CondLikelihoodShPtr getCondLikePtr(
  EdgeEndpoints edge) /**< is the edge specifying the focal node and the focal node neighbor */
	{
	TreeNode * actual_child = edge.getActualChild();
	if (actual_child == edge.getFocalNode())
		{
		// focal node F is a child of N, the focal neighbor (likelihood root is somewhere below N)
		//
		//		\	/
		//		 \ /
		//		  F
		//	 \	 / <-- filial CLA of focal node is returned
		//	  \ /
		//	   N
		//	  / 
		//
		PHYCAS_ASSERT(actual_child->IsInternal());
		InternalData * child_internal_data = actual_child->GetInternalData();
		PHYCAS_ASSERT(child_internal_data != NULL);
		return child_internal_data->getChildCondLikePtr();
		}

	// focal neighbor N is a child of the focal node F
	PHYCAS_ASSERT(actual_child == edge.getFocalNeighbor());
	if (actual_child->IsInternal())
		{
		// focal neighbor N is an internal node (likelihood root is somewhere above N)
		//
		//		\	/
		//		 \ /
		//		  N
		//	 \	 /
		//	  \ / <-- parental CLA of focal neighbor is returned
		//	   F
		//	  / 
		//
		InternalData * child_internal_data = actual_child->GetInternalData();
		PHYCAS_ASSERT(child_internal_data != NULL);
		return child_internal_data->getParentalCondLikePtr();		
		}

	// focal neighbor N is a tip node (N equals the likelihood root in this case)
	//
	//		  N
	//	 \	 /
	//	  \ / <-- parental CLA of focal neighbor is returned
	//	   F
	//	  / 
	//
	TipData * child_tip_data = actual_child->GetTipData();
	PHYCAS_ASSERT(child_tip_data != NULL);
	return child_tip_data->getParentalCondLikePtr();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a pointer to the CondLikelihood object that would be used to compute the likelihood conditional on a
|	particular state being assigned to the focal node of `edge' when `edge' focal neighbor is the likelihood root.
|	This is a const version of the corresponding getCondLikePtr function, to be used when it can be assumed that the
|	conditional likelihood arrays are up-to-date.
*/
ConstCondLikelihoodShPtr getValidCondLikePtr(
  ConstEdgeEndpoints edge) /**< is the edge comprising the focal node and focal node neighbor */
	{
#if 1
	const TreeNode * actual_child = edge.getActualChild();
	if (actual_child == edge.getFocalNode())
		{
		// focal node F is a child of N, the focal neighbor (likelihood root is somewhere below N)
		//
		//		\	/
		//		 \ /
		//		  F
		//	 \	 / <-- filial CLA of focal node is returned
		//	  \ /
		//	   N
		//	  / 
		//
		PHYCAS_ASSERT(actual_child->IsInternal());
		const InternalData * child_internal_data = actual_child->GetInternalData();
		PHYCAS_ASSERT(child_internal_data != NULL);
		return child_internal_data->getValidChildCondLikePtr();
		}

	// focal neighbor N is a child of the focal node F
	PHYCAS_ASSERT(actual_child == edge.getFocalNeighbor());
	if (actual_child->IsInternal())
		{
		// focal neighbor N is an internal node (likelihood root is somewhere above N)
		//
		//		\	/
		//		 \ /
		//		  N
		//	 \	 /
		//	  \ / <-- parental CLA of focal neighbor is returned
		//	   F
		//	  / 
		//
		const InternalData * child_internal_data = actual_child->GetInternalData();
		PHYCAS_ASSERT(child_internal_data != NULL);
		return child_internal_data->getValidParentalCondLikePtr();		
		}

	// focal neighbor N is a tip node (N equals the likelihood root in this case)
	//
	//		  N
	//	 \	 /
	//	  \ / <-- parental CLA of focal neighbor is returned
	//	   F
	//	  / 
	//
	const TipData * child_tip_data = actual_child->GetTipData();
	PHYCAS_ASSERT(child_tip_data != NULL);
	return child_tip_data->getValidParentalCondLikePtr();
#else
	const TreeNode * c = edge.getActualChild();
	if (edge.getFocalNode() == c)
		{
		PHYCAS_ASSERT(c->IsInternal());
		const InternalData * childInternalData = c->GetInternalData();
		PHYCAS_ASSERT(childInternalData != NULL);
		return childInternalData->getValidChildCondLike();
		}
	// moving up the tree in calculations (root to leaves).
	PHYCAS_ASSERT(c == edge.getFocalNeighbor());
	if (c->IsInternal())
		{
		const InternalData * childInternalData = c->GetInternalData();
		PHYCAS_ASSERT(childInternalData != NULL);
		return childInternalData->getValidParentalCondLike();		
		}
	const TipData * childTipData = c->GetTipData();
	PHYCAS_ASSERT(childTipData != NULL);
	return childTipData->getValidParentalCondLike();
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the value of the `no_data' data member to true.
*/
void TreeLikelihood::setNoData()
	{
	no_data = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the value of the `no_data' data member to false.
*/
void TreeLikelihood::setHaveData()
	{
	no_data = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the `num_patterns' data member.
*/
void TreeLikelihood::setNPatterns(
  unsigned nPatterns)	/**< is the number of patterns */
	{
	num_patterns = nPatterns;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the `model' data member, replacing the model defined in the constructor.
*/
void TreeLikelihood::replaceModel(
  ModelShPtr m)
	{
	model = m;
	recalcRelativeRates();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the values of the `num_states' and `num_rates' data members according to the model, then calls the 
|	recalcRatesAndProbs function of the model to force recalculation of the `rate_means' and `rate_probs' vectors. 
|	Should be called after changing the number of rate categories, the gamma shape parameter, or the pinvar parameter of
|	the model. Note that trees on which likelihoods need to be calculated also need to be re-equipped by calling 
|	prepareForLikelihood if the number of rate categories changes (not done by this function). 
*/
void TreeLikelihood::recalcRelativeRates()
	{
	num_states = model->getNStates();
	num_rates = model->getNRatesTotal();
	model->recalcRatesAndProbs(rate_means, rate_probs); //POL_BOOKMARK recalcRatesAndProbs call
	likelihood_rate_site.resize(num_rates*num_patterns, 0.0);
#if POLPY_NEWWAY	//CLAShPtr
	if (!no_data)
		cla_pool->setCondLikeDimensions(num_patterns, num_rates, num_states);
#else
	if (!no_data)
		cla_pool.setCondLikeDimensions(num_patterns, num_rates, num_states);
#endif
	underflow_policy.setDimensions(num_patterns, num_rates, num_states);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `cla_pool' data member.
*/
#if POLPY_NEWWAY	//CLAShPtr
const CondLikelihoodStorageShPtr TreeLikelihood::getCLAStorage() const
#else
const CondLikelihoodStorage & TreeLikelihood::getCLAStorage() const
#endif
	{
	return cla_pool;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `nTaxa' data member.
*/
unsigned TreeLikelihood::getNTaxa() const
	{
	return nTaxa;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `num_patterns' data member.
*/
unsigned TreeLikelihood::getNPatterns() const
	{
	return num_patterns;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `num_rates' data member. Assumes that the length of the `rate_means' vector 
|	equals the value of the `num_rates' data member.
*/
unsigned TreeLikelihood::getNRatesTotal() const
	{
	PHYCAS_ASSERT(rate_means.size() == num_rates);
	return num_rates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `num_states' data member.
*/
unsigned TreeLikelihood::getNStates() const
	{
	return num_states;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns a copy of the (shared_ptr) data member `model'.
*/
ModelShPtr TreeLikelihood::getModel() const
	{
	return model;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `state_list'.
*/
const VecStateList & TreeLikelihood::getStateList() const
	{
	return state_list;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `state_list_pos'.
*/
const VecStateListPos & TreeLikelihood::getStateListPos() const
	{
	return state_list_pos;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `rate_means'.
*/
const std::vector<double> & TreeLikelihood::getRateMeans() const
	{
	return rate_means; //POL_BOOKMARK TreeLikelihood::getRateMeans
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `rate_probs'.
*/
const std::vector<double> & TreeLikelihood::getRateProbs() const
	{
	return rate_probs;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `category_boundaries'.
*/
std::vector<double> TreeLikelihood::getCategoryLowerBoundaries() const
	{
	std::vector<double> tmp_means;
	std::vector<double> returned_boundaries;
	model->recalcGammaRatesAndBoundaries(tmp_means, returned_boundaries);
	return returned_boundaries;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls TreeLikelihood::simulateImpl specifying that the transition probabilities should be recalculated (or 
|	calculated the first time) before beginning simulations. After TreeLikelihood::simulateFirst is called once, 
|	TreeLikelihood::simulate can be called many more times to generate more data sets using the same transition 
|	probabilities.
*/
void TreeLikelihood::simulateFirst(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar)
	{
	simulateImpl(sim_data, t, rng, nchar, true);	// true means recalculate transition probabilities
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls TreeLikelihood::simulateImpl specifying that the transition probabilities should NOT be recalculated. Only
|	call this function after calling TreeLikelihood::simulateFirst at least once, otherwise the transition probabilities
|	will contain garbage.
*/
void TreeLikelihood::simulate(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar)
	{
	simulateImpl(sim_data, t, rng, nchar, false);	// false means do not recalculate transition probabilities
	}

// **************************************************************************************
// ***** Former TreeLikelihood inlines (end) ********************************************
// **************************************************************************************

/*----------------------------------------------------------------------------------------------------------------------
|	Specifies the (internal) node to use as the likelihood root (it will be stored in the `likelihood_root' data member). 
|	The likelihood root is separate from the actual root of the tree, and specifies the node used when harvesting the 
|	log-likelihood from surrounding conditional likelihood arrays (CLAs). It behooves one to set the likelihood root to 
|	that node requiring the fewest CLA recalculations. Specifying NULL for the likelihood root will result in the 
|	unconditional recalculation of all CLAs in the entire tree, and subsequently the likelihood root will be set to the
|	subroot node of the tree (the only child of the root).
*/
void TreeLikelihood::useAsLikelihoodRoot(
  TreeNode * nd)
	{
	likelihood_root = nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Updates the conditional likelihood array for `nd' in a direction away from `avoid'. All adjacent nodes (other than 
|	`avoid') are assumed to be valid.
*/
void TreeLikelihood::refreshCLA(TreeNode & nd, const TreeNode * avoid)
	{
	if (nd.IsTip())
		return;

	PHYCAS_ASSERT(avoid != NULL);
	//@MTH-NESCENT this const_cast should be safe because we aren't relying on const-ness of CLA's but it needs to be revisited
	//@POL-NESCENT Mark, if we are going to immediately cast away the const on avoid, why do we need to make it const in the first place?
	CondLikelihoodShPtr ndCondLike = getCondLikePtr(&nd, const_cast<TreeNode *>(avoid)); 
	TreeNode * parent = nd.GetParent();
	TreeNode * lChild = nd.GetLeftChild();
	PHYCAS_ASSERT(parent != NULL);
	PHYCAS_ASSERT(lChild != NULL);

	// The first neighbor can either be the parent or the leftmost child. The second neighbor must be a child, but
	// which child depends on the first neighbor. Third and subsequent neighbors are always next sibs.
	TreeNode * firstNeighbor = NULL;
	TreeNode * secondNeighbor = NULL;
	double firstEdgeLen = 0.0;
	const bool movingTowardLeaves = !(parent == avoid);
	if (movingTowardLeaves)
		{
		// A child of nd is closer to the likelihood root than nd
		firstNeighbor = parent;
		firstEdgeLen = nd.GetEdgeLen();
		secondNeighbor = (lChild == avoid ? lChild->GetRightSib() : lChild);
		}
	else
		{
		// The parent of nd is closer to the likelihood root than nd
		firstNeighbor = lChild;
		firstEdgeLen = lChild->GetEdgeLen();
		secondNeighbor = lChild->GetRightSib();
		}
	
	if (firstNeighbor->IsTip())
		{
		TipData & firstTD = *(firstNeighbor->GetTipData());
		calcPMatTranspose(firstTD.getTransposedPMatrices(), firstTD.getConstStateListPos(), firstEdgeLen);
		if (secondNeighbor->IsTip())
			{
			// 1. both neighbors are tips
			TipData & secondTD = *(secondNeighbor->GetTipData());
			calcPMatTranspose(secondTD.getTransposedPMatrices(), secondTD.getConstStateListPos(), secondNeighbor->GetEdgeLen());
			calcCLATwoTips(*ndCondLike, firstTD, secondTD);
			}
		else
			{
			// 2. first neighbor is a tip, but second is an internal node
			InternalData & secondID = *(secondNeighbor->GetInternalData());
			calcPMat(secondID.getPMatrices(), secondNeighbor->GetEdgeLen());
			CondLikelihoodShPtr secCL = getCondLikePtr(secondNeighbor, &nd);
			ConstPMatrices secPMat = secondID.getConstPMatrices();
			calcCLAOneTip(*ndCondLike, firstTD, secPMat, *secCL);
			}
		}
	else
		{
		InternalData & firstID = *(firstNeighbor->GetInternalData());
		calcPMat(firstID.getPMatrices(), firstEdgeLen);
		const CondLikelihood & firCL = *getCondLikePtr(firstNeighbor, &nd);
		ConstPMatrices firPMat = firstID.getConstPMatrices();	
		if (secondNeighbor->IsTip())
			{
			// 3. first neighbor internal node, but second is a tip
			TipData & secondTD = *(secondNeighbor->GetTipData());
			calcPMatTranspose(secondTD.getTransposedPMatrices(), secondTD.getConstStateListPos(), secondNeighbor->GetEdgeLen());
			calcCLAOneTip(*ndCondLike, secondTD, firPMat, firCL);
			}
		else
			{
			// 4. both neighbors are internal nodes
			InternalData & secondID = *(secondNeighbor->GetInternalData());
			calcPMat(secondID.getPMatrices(), secondNeighbor->GetEdgeLen());
			const CondLikelihood & secCL = *getCondLikePtr(secondNeighbor, &nd);
			ConstPMatrices secPMat = secondID.getConstPMatrices();
			calcCLANoTips(*ndCondLike, firPMat, firCL, secPMat, secCL);
			}
		}

	// Deal with possible polytomy in which secondNeighbor has siblings
	for (TreeNode * currNd = secondNeighbor->GetRightSib(); currNd != NULL; currNd = currNd->GetRightSib())
		{
		if (currNd != avoid)
			{
			if (currNd->IsTip())
				{
				TipData & currTD = *(currNd->GetTipData());
				calcPMatTranspose(currTD.getTransposedPMatrices(), currTD.getConstStateListPos(), currNd->GetEdgeLen());
				conditionOnAdditionalTip(*ndCondLike, currTD);
				}
			else
				{
				InternalData & currID = *(currNd->GetInternalData());
				calcPMat(currID.getPMatrices(), currNd->GetEdgeLen());
				const CondLikelihood & currCL = *getCondLikePtr(currNd, &nd);
				ConstPMatrices currPMat = currID.getConstPMatrices();
				conditionOnAdditionalInternal(*ndCondLike, currPMat, currCL);
				}
			}
		}

#if 0
	// Turn this section on for debugging purposes only! Walks through conditional likelihood arrays
	// for this node looking for any conditional likelihood that is negative. Negative cond. likes
	// have shown up in the past (18 Oct 2007) for the GTR model, which results in a NaN (represented
	// as -1.#IND in the VC compiler when the log of the site likelihood is taken.
	LikeFltType * cla = ndCondLike->getCLA();
	for (unsigned i = 0; i < num_rates; ++i)
		{
		for (unsigned j = 0; j < num_patterns; ++j)
			{
			for (unsigned k = 0; k < num_states; ++k)
				{
				if (*cla++ < 0.0)
					{
					std::cerr << "*** Doh! cla negative for rate = " << i << ", pattern = " << j << ", state = " << k << " ***" << std::endl;
					std::exit(0);
					}
				}
			}
		}

#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the conditional likelihood of refNd is up to date for calculations centered at some effective root 
|	node (neighborCloserToEffectiveRoot will be a node adjacent to refNd, but closer than refNd to the the effective 
|	root). Called in a context in which neighborCloserToEffectiveRoot is requesting that all of its neighbors update 
|	their likelihood temporaries.
*/
bool TreeLikelihood::isValid(const TreeNode * refNd, const TreeNode * neighborCloserToEffectiveRoot)
	{
	PHYCAS_ASSERT(refNd != NULL);
	PHYCAS_ASSERT(neighborCloserToEffectiveRoot != NULL);
	PHYCAS_ASSERT(refNd != neighborCloserToEffectiveRoot);
	if (refNd->IsTip())
		{
		// tip nodes always return true because they must be the child of neighborCloserToEffectiveRoot
		// and they do not hold filial CLAs
		return true;
		}
	else
		{
		// refNd is internal
		if (refNd->GetParentConst() == neighborCloserToEffectiveRoot)
			{
			// refNd is the child of neighborCloserToEffectiveRoot
			// If refNd has a filial CLA, then return true because that means refNd is valid
			const InternalData * id = refNd->GetInternalData();
			return (id->childWorkingCLA);			
			}
		else
			{
			// neighborCloserToEffectiveRoot is the child of refNd
			// If neighborCloserToEffectiveRoot has a parental CLA, then return true because 
			// that means refNd is valid
			const InternalData * id = neighborCloserToEffectiveRoot->GetInternalData();
			return (id->parWorkingCLA);			
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Unconditionally invalidates parental (and, if ref_nd is internal, filial) conditional likelihood arrays and also
|	removes cached conditional likelihood arrays if present. This function always returns false so it can be used in 
|	conjunction with effective_postorder_edge_iterator to invalidate every CLA in the entire tree and ensure that there
|	are also no cached CLAs as well.
*/
bool TreeLikelihood::invalidateBothEndsDiscardCache(TreeNode * ref_nd, TreeNode * /* unused */)
	{
	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();
		if (td != NULL)
			{
			// Invalidate the parental CLAs if they exist
			if (td->parWorkingCLA)
				{
#if POLPY_NEWWAY	//CLAShPtr
				cla_pool->putCondLikelihood(td->parWorkingCLA);
#else
				cla_pool.putCondLikelihood(td->parWorkingCLA);
#endif
				td->parWorkingCLA.reset();
				}
			// Remove cached parental CLAs if they exist
			if (td->parCachedCLA)
				{
#if POLPY_NEWWAY	//CLAShPtr
				cla_pool->putCondLikelihood(td->parCachedCLA);
#else
				cla_pool.putCondLikelihood(td->parCachedCLA);
#endif
				td->parCachedCLA.reset();
				}
			}
		}
	else
		{
		// Internal nodes have both parental and filial CLAs
		InternalData * id = ref_nd->GetInternalData();
		if (id != NULL)
			{
			// Invalidate the parental CLAs if they exist
			if (id->parWorkingCLA)
				{
#if POLPY_NEWWAY	//CLAShPtr
				cla_pool->putCondLikelihood(id->parWorkingCLA);
#else
				cla_pool.putCondLikelihood(id->parWorkingCLA);
#endif
				id->parWorkingCLA.reset();
				}
			// Remove cached parental CLAs if they exist
			if (id->parCachedCLA)
				{
#if POLPY_NEWWAY	//CLAShPtr
				cla_pool->putCondLikelihood(id->parCachedCLA);
#else
				cla_pool.putCondLikelihood(id->parCachedCLA);
#endif
				id->parCachedCLA.reset();
				}

			// Invalidate the filial CLAs if they exist
			if (id->childWorkingCLA)
				{
#if POLPY_NEWWAY	//CLAShPtr
				cla_pool->putCondLikelihood(id->childWorkingCLA);
#else
				cla_pool.putCondLikelihood(id->childWorkingCLA);
#endif
				id->childWorkingCLA.reset();
				}
			// Remove cached filial CLAs if they exist
			if (id->childCachedCLA)
				{
#if POLPY_NEWWAY	//CLAShPtr
				cla_pool->putCondLikelihood(id->childCachedCLA);
#else
				cla_pool.putCondLikelihood(id->childCachedCLA);
#endif
				id->childCachedCLA.reset();
				}
			}
		}

	return false;
	}

#if POLPY_NEWWAY	//CLAShPtr
CondLikelihoodStorageShPtr TreeLikelihood::getCondLikelihoodStorage()
#else
CondLikelihoodStorage & TreeLikelihood::getCondLikelihoodStorage()
#endif
	{
	return cla_pool;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Unconditionally invalidate parental (and, if ref_nd is internal, filial) conditional likelihood arrays. Always 
|	returns false so it can be used in conjunction with effective_postorder_edge_iterator to invalidate every CLA in
|	the entire tree.
*/
bool TreeLikelihood::invalidateBothEnds(TreeNode * ref_nd, TreeNode * /* unused */)
	{
	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();

		// Invalidate the parental CLAs if they exist
		if (td->parWorkingCLA)
			{
#if POLPY_NEWWAY	//CLAShPtr
			if (td->parCachedCLA)
				cla_pool->putCondLikelihood(td->parCachedCLA);
#else
			if (td->parCachedCLA)
				cla_pool.putCondLikelihood(td->parCachedCLA);
#endif
			td->parCachedCLA = td->parWorkingCLA;
			td->parWorkingCLA.reset();
			}
		}
	else
		{
		// Internal nodes have both parental and filial CLAs
		InternalData * id = ref_nd->GetInternalData();

		// Invalidate the parental CLAs if they exist
		if (id->parWorkingCLA)
			{
#if POLPY_NEWWAY	//CLAShPtr
			if (id->parCachedCLA)
				cla_pool->putCondLikelihood(id->parCachedCLA);
#else
			if (id->parCachedCLA)
				cla_pool.putCondLikelihood(id->parCachedCLA);
#endif
			id->parCachedCLA = id->parWorkingCLA;
			id->parWorkingCLA.reset();
			}

		// Invalidate the filial CLAs if they exist
		if (id->childWorkingCLA)
			{
#if POLPY_NEWWAY	//CLAShPtr
			if (id->childCachedCLA)
				cla_pool->putCondLikelihood(id->childCachedCLA);
#else
			if (id->childCachedCLA)
				cla_pool.putCondLikelihood(id->childCachedCLA);
#endif
			id->childCachedCLA = id->childWorkingCLA;
			id->childWorkingCLA.reset();
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function used in conjunction with effective_postorder_edge_iterator to invalidate the appropriate conditional 
|	likelihood arrays (CLAs) of all nodes starting with a focal node. For nodes that are descendants of the focal node
|	(descendant means that a node can be found using only left child and right sib pointers starting from the focal
|	node), it is the parental CLAs that are invalidated. For a node that is an ancestor of the focal node, it is the 
|	filial CLAs that are invalidated. For all other nodes (e.g. on independent lineages derived from an ancestral node), 
|	it is the parental CLAs that are invalidated. This pattern of invalidation ensures that the likelihood will be 
|	correctly computed using any node in the tree as the likelihood root. For any likelihood root n, it now appears as 
|	though an edge length just on the other side of the focal node from n has been changed, necessitating the 
|	recalculation of all CLAs starting from that point back toward the likelihood root. The return value indicates
|	whether or not the iterator should continue into the subtree defined by `refNd' (away from 
|	`neighborCloserToEffectiveRoot'). This function always returns false because the goal is to invalidate every node
|	that needs to be invalidated.
*/
bool TreeLikelihood::invalidateNode(TreeNode * ref_nd, TreeNode * neighbor_closer_to_likelihood_root)
	{
	if (ref_nd->IsTip() && !ref_nd->IsTipRoot())
		{
		TipData * td = ref_nd->GetTipData();
		if (!td->parWorkingCLA)
			return false;
#if POLPY_NEWWAY	//CLAShPtr
		if (td->parCachedCLA)
			cla_pool->putCondLikelihood(td->parCachedCLA);
#else
		if (td->parCachedCLA)
			cla_pool.putCondLikelihood(td->parCachedCLA);
#endif
		td->parCachedCLA = td->parWorkingCLA;
		td->parWorkingCLA.reset();
		}
	else
		{
		if (ref_nd->GetParent() == neighbor_closer_to_likelihood_root)
			{
			// ref_nd is the actual child
			InternalData * id = ref_nd->GetInternalData();
			if (!id->parWorkingCLA)
				return false;
#if POLPY_NEWWAY	//CLAShPtr
			if (id->parCachedCLA)
				cla_pool->putCondLikelihood(id->parCachedCLA);
#else
			if (id->parCachedCLA)
				cla_pool.putCondLikelihood(id->parCachedCLA);
#endif
			id->parCachedCLA = id->parWorkingCLA;
			id->parWorkingCLA.reset();
			}
		else
			{
			// neighbor_closer_to_likelihood_root is the actual child
			if (neighbor_closer_to_likelihood_root->IsTip())
				return false;
			// Either ref_nd is not a tip, or it is the tip serving as the root node
			InternalData * id = neighbor_closer_to_likelihood_root->GetInternalData();
			if (!id->childWorkingCLA)
				return false;
#if POLPY_NEWWAY	//CLAShPtr
			if (id->childCachedCLA)
				cla_pool->putCondLikelihood(id->childCachedCLA);
#else
			if (id->childCachedCLA)
				cla_pool.putCondLikelihood(id->childCachedCLA);
#endif
			id->childCachedCLA = id->childWorkingCLA;
			id->childWorkingCLA.reset();
			}
		}
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Invalidates all conditional likelihood arrays (CLAs) pointing toward the supplied `focal_node'. That is, the
|	parental CLA of a child of `focal_node' would be invalidated, as would be the filial CLA for the parent of the
|	`focal_node'. The invalidation propagates throughout the tree from `focal_node' so that regardless of the node
|	serving as the likelihood root, the correct CLAs will be recalculated the next time the log-likelihood is 
|	recomputed. If the edge of `focal_node' has been changed, the function invalidateBothEnds(`focal_node') should also
|	be called to invalidate the one remaining CLA not invalidated by this function.
*/
void TreeLikelihood::invalidateAwayFromNode(
  TreeNode & focal_node)		/**< invalidation of conditional likelihood arrays will proceed outward from this node */
	{
	NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::invalidateNode, this, _1, _2);
	effective_postorder_edge_iterator(&focal_node, validFunctor); // need only construct unnamed iterator object
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Unconditionally discards cached parental (and, if ref_nd is internal, filial) conditional likelihood arrays if
|	present. This function always returns false so it can be used in conjunction with effective_postorder_edge_iterator
|	to discard every cached CLA in the entire tree (something that needs to be done when any move is accepted).
*/
bool TreeLikelihood::discardCacheBothEnds(TreeNode * ref_nd, TreeNode * /* unused */)
	{
	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();

		// Remove cached parental CLAs if they exist
		if (td->parCachedCLA)
			{
#if POLPY_NEWWAY	//CLAShPtr
			cla_pool->putCondLikelihood(td->parCachedCLA);
#else
			cla_pool.putCondLikelihood(td->parCachedCLA);
#endif
			td->parCachedCLA.reset();
			}
		}
	else
		{
		// Internal nodes have both parental and filial CLAs
		InternalData * id = ref_nd->GetInternalData();

		// Remove cached parental CLAs if they exist
		if (id->parCachedCLA)
			{
#if POLPY_NEWWAY	//CLAShPtr
			cla_pool->putCondLikelihood(id->parCachedCLA);
#else
			cla_pool.putCondLikelihood(id->parCachedCLA);
#endif
			id->parCachedCLA.reset();
			}

		// Remove cached filial CLAs if they exist
		if (id->childCachedCLA)
			{
#if POLPY_NEWWAY	//CLAShPtr
			cla_pool->putCondLikelihood(id->childCachedCLA);
#else
			cla_pool.putCondLikelihood(id->childCachedCLA);
#endif
			id->childCachedCLA.reset();
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores parental conditional likelihood arrays from cache and stores the existing conditional likelihood arrays 
|	for later use. This function always returns false so it can be used in conjunction with 
|	effective_postorder_edge_iterator to restore every parental CLA from cache in the entire tree. If `ref_nd' 
|	has a working parental CLA but no corresponding cached parental CLA, then the working parental CLA is discarded. 
|	This is done because this working parental CLA must have been calculated while the tree was in a proposed state, 
|	and thus would be invalid if the tree were reverted.
*/
bool TreeLikelihood::restoreFromCacheParentalOnly(TreeNode * ref_nd, TreeNode * /* unused */)
	{
	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();

		// Store working CLA in any case
#if POLPY_NEWWAY	//CLAShPtr
		if (td->parWorkingCLA)
			cla_pool->putCondLikelihood(td->parWorkingCLA);
#else
		if (td->parWorkingCLA)
			cla_pool.putCondLikelihood(td->parWorkingCLA);
#endif
		td->parWorkingCLA.reset();

		// Move cached to working if cached exists
		if (td->parCachedCLA)
			{
			td->parWorkingCLA = td->parCachedCLA;
			td->parCachedCLA.reset();
			}
		}
	else
		{
		InternalData * id = ref_nd->GetInternalData();

		// Store working CLA in any case
#if POLPY_NEWWAY	//CLAShPtr
		if (id->parWorkingCLA)
			cla_pool->putCondLikelihood(id->parWorkingCLA);
#else
		if (id->parWorkingCLA)
			cla_pool.putCondLikelihood(id->parWorkingCLA);
#endif
		id->parWorkingCLA.reset();

		// Move cached to working if cached exists
		if (id->parCachedCLA)
			{
			id->parWorkingCLA = id->parCachedCLA;
			id->parCachedCLA.reset();
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores parental (and, if ref_nd is internal, filial) conditional likelihood arrays from cache and stores the
|	existing conditional likelihood arrays for later use. This function always returns false so it can be used in 
|	conjunction with effective_postorder_edge_iterator to restore every CLA from cache in the entire tree. If `ref_nd' 
|	has a working CLA but no corresponding cached CLA, then the working CLA is discarded. This is done because this 
|	working CLA must have been calculated while the tree was in a proposed state, and thus would be invalid if the tree
|	were reverted.
*/
bool TreeLikelihood::restoreFromCacheBothEnds(TreeNode * ref_nd, TreeNode * /* unused */)
	{
	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();

		// Store working CLA in any case
#if POLPY_NEWWAY	//CLAShPtr
		if (td->parWorkingCLA)
			cla_pool->putCondLikelihood(td->parWorkingCLA);
#else
		if (td->parWorkingCLA)
			cla_pool.putCondLikelihood(td->parWorkingCLA);
#endif
		td->parWorkingCLA.reset();

		// Move cached to working if cached exists
		if (td->parCachedCLA)
			{
			td->parWorkingCLA = td->parCachedCLA;
			td->parCachedCLA.reset();
			}
		}
	else
		{
		// Internal nodes have both parental and filial CLAs
		InternalData * id = ref_nd->GetInternalData();

		// Restore parental CLA from cache

		// Store working CLA in any case
#if POLPY_NEWWAY	//CLAShPtr
		if (id->parWorkingCLA)
			cla_pool->putCondLikelihood(id->parWorkingCLA);
#else
		if (id->parWorkingCLA)
			cla_pool.putCondLikelihood(id->parWorkingCLA);
#endif
		id->parWorkingCLA.reset();

		// Move cached to working if cached exists
		if (id->parCachedCLA)
			{
			id->parWorkingCLA = id->parCachedCLA;
			id->parCachedCLA.reset();
			}

		// Restore filial CLA from cache

		// Store working CLA in any case
#if POLPY_NEWWAY	//CLAShPtr
		if (id->childWorkingCLA)
			cla_pool->putCondLikelihood(id->childWorkingCLA);
#else
		if (id->childWorkingCLA)
			cla_pool.putCondLikelihood(id->childWorkingCLA);
#endif
		id->childWorkingCLA.reset();

		// Move cached to working if cached exists
		if (id->childCachedCLA)
			{
			id->childWorkingCLA = id->childCachedCLA;
			id->childCachedCLA.reset();
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function used in conjunction with effective_postorder_edge_iterator to restore the conditional likelihood arrays 
|	(CLAs) of all nodes starting with the supplied focal node `ref_nd'. For nodes that are descendants of the focal node
|	(descendant means that a node can be found using only left child and right sib pointers starting from the focal
|	node), it is the parental CLAs that are restored. For a node that is an ancestor of the focal node, it is the 
|	filial CLAs that are restored. For all other nodes (e.g. on independent lineages derived from an ancestral node), 
|	it is the parental CLAs that are restored. This function always returns false because the goal is to restore every 
|	node that needs to be invalidated. If a node is has a working CLA but no corresponding cached CLA, then the working
|	CLA is discarded. This is done because this working CLA must have been calculated while the tree was in a proposed
|	state, and thus would be invalid if the tree were reverted.
*/
bool TreeLikelihood::restoreFromCacheNode(TreeNode * ref_nd, TreeNode * neighbor_closer_to_likelihood_root)
	{
	if (ref_nd->IsTip() && !ref_nd->IsTipRoot())
		{
		TipData * td = ref_nd->GetTipData();

		// Store the working CLA in any case
#if POLPY_NEWWAY	//CLAShPtr
		if (td->parWorkingCLA)
			cla_pool->putCondLikelihood(td->parWorkingCLA);
#else
		if (td->parWorkingCLA)
			cla_pool.putCondLikelihood(td->parWorkingCLA);
#endif
		td->parWorkingCLA.reset();

		// Move cached to working if there is a cached CLA
		if (td->parCachedCLA)
			{
			td->parWorkingCLA = td->parCachedCLA;
			td->parCachedCLA.reset();
			}
		}
	else
		{
		if (ref_nd->GetParent() == neighbor_closer_to_likelihood_root)
			{
			InternalData * id = ref_nd->GetInternalData();

			// Store the working CLA in any case
#if POLPY_NEWWAY	//CLAShPtr
			if (id->parWorkingCLA)
				cla_pool->putCondLikelihood(id->parWorkingCLA);
#else
			if (id->parWorkingCLA)
				cla_pool.putCondLikelihood(id->parWorkingCLA);
#endif
			id->parWorkingCLA.reset();

			// Move cached to working if there is a cached CLA
			if (id->parCachedCLA)
				{
				id->parWorkingCLA = id->parCachedCLA;
				id->parCachedCLA.reset();
				}
			}
		else
			{
			// neighbor_closer_to_likelihood_root is the actual child
			if (neighbor_closer_to_likelihood_root->IsTip())
				return false;

			// Either ref_nd is not a tip, or it is the tip serving as the root node
			InternalData * id = neighbor_closer_to_likelihood_root->GetInternalData();

			// Store the working CLA in any case
#if POLPY_NEWWAY	//CLAShPtr
			if (id->childWorkingCLA)
				cla_pool->putCondLikelihood(id->childWorkingCLA);
#else
			if (id->childWorkingCLA)
				cla_pool.putCondLikelihood(id->childWorkingCLA);
#endif
			id->childWorkingCLA.reset();

			// Move cached to working if there is a cached CLA
			if (id->childCachedCLA)
				{
				id->childWorkingCLA = id->childCachedCLA;
				id->childCachedCLA.reset();
				}
			}
		}
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores all conditional likelihood arrays (CLAs) pointing toward the supplied `focal_node' from cache. That is, the
|	parental CLA of a child of `focal_node' would be restored, as would be the filial CLA for the parent of the 
|	`focal_node'. It may be necessary to call the function invalidateBothEnds(`focal_node') also to restore the one 
|	remaining CLA not restored by this function.
*/
void TreeLikelihood::restoreFromCacheAwayFromNode(
  TreeNode & focal_node)		/**< restoration of cached conditional likelihood arrays will proceed outward from this node */
	{
	NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::restoreFromCacheNode, this, _1, _2);
	effective_postorder_edge_iterator(&focal_node, validFunctor); // need only construct unnamed iterator object
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Discards all cached conditional likelihood arrays (CLAs) pointing toward the supplied `focal_node'. That is, the
|	cached parental CLA of a child of `focal_node' would be discarded, as would be the cached filial CLA for the parent 
|	of the `focal_node'. None of the working CLAs are affected. It may be necessary to call the function 
|	discardCacheBothEnds(`focal_node') also to discard the cache from the two remaining CLAs not affected by this 
|	function.
*/
void TreeLikelihood::discardCacheAwayFromNode(
  TreeNode & focal_node)		/**< restoration of cached conditional likelihood arrays will proceed outward from this node */
	{
	NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::discardCacheBothEnds, this, _1, _2);
	effective_postorder_edge_iterator(&focal_node, validFunctor); // need only construct unnamed iterator object
	}

#if 0
class CalcTransitionMatrixForOneNode : public std::unary_function<TreeNode &, void>
	{
	private:
		TreeLikelihood & treelike;

	public:
		CalcTransitionMatrixForOneNode(TreeLikelihood & tl) : treelike(tl) {}
		void operator()(TreeNode & nd)
			{
#			error do not use unless root node special case is taken into account
			if (!nd.IsTipRoot())
				{
				double edge_len = nd.GetEdgeLen();
				if (nd.IsTip())
					{
					TipData & ndTD = *(nd.GetTipData());
					treelike.calcTMatForSim(ndTD, edge_len);
					}
				else
					{
					InternalData & ndID = *(nd.GetInternalData());
					treelike.calcPMat(ndID, edge_len);
					}
				}
			}
	};
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Simulates data for `nchar' characters using the current model and edge lengths. Transition matrices are recomputed
|	if `refresh_probs' is true; if `refresh_probs' is false, assumes transition probability matrices are up-to-date. 
|	Only the `num_states' primary states are generated, and the state created for each node is stored initially in the 
|	`state' data member of the TipData or InternalData structure associated with the tip or internal node, respectively.
|	This function expects these data structures to be in place (the TreeLikelihood::prepareForSimulation and 
|	TreeLikelihood::prepareForLikelihood functions both perform this task). As states are generated for tip nodes, they
|	are copied into a temporary pattern inside `sim_data', and this pattern is then inserted into a pattern map inside 
|	`sim_data' when it is completed. Assumes that `nchar' is greater than zero and that the shared pointers `sim_data' 
|	and `rng' actually point to real objects.
|	
|	The following makes use of T matrices, which are transposed, augmented transition matrices. These are transposed 
|	because the "from" states form the columns rather than the rows. They are augmented because, ordinarily, there are
|	additional rows corresponding to ambiguities observed in some tip nodes. With simulated data, however, there are 
|	never any ambiguities, so in this case T matrices are nothing more than transposed transition matrices. 
|	Important: note that T matrices are only used in the TipData structures; InternalData structures store normal 
|	untransposed transition probability matrices in which the rows form the "from" states. 
|	
|	Here is an example of a T matrix created using the JC model and an edge length equal to 0.1:
|>	
|				|--------------- from state -------------|
|					0		   1		  2			 3
|					A		   C		  G			 T
|	t  0   A	 0.90638	0.03121	   0.03121	  0.03121
|	o  1   C	 0.03121	0.90638	   0.03121	  0.03121
|	   2   G	 0.03121	0.03121	   0.90638	  0.03121
|	s  3   T	 0.03121	0.03121	   0.03121	  0.90638
|	t  4   N	 1.00000	1.00000	   1.00000	  1.00000 \
|	a  5 {GT}	 0.06241	0.06241	   0.93759	  0.93759  | These rows not present if prepareForSimulation function
|	t  6 {ACT}	 0.96879	0.96879	   0.09362	  0.96879  | was used to create the TipData structures
|	e  7 {AG}	 0.93757	0.06241	   0.93759	  0.06241 /
|>
|	The `pMatrixTranspose' data member in TipData structures holds the array of T matrices (one T matrix for each 
|	rate category).
*/
void TreeLikelihood::simulateImpl(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar, bool refresh_probs)
	{
	PHYCAS_ASSERT(sim_data);
	PHYCAS_ASSERT(rng);
	PHYCAS_ASSERT(nchar > 0);

	// Recalculate transition probabilities if requested
	if (refresh_probs)
		{
		preorder_iterator nd = t->begin();

		// First preorder node is the root node and represents a special case
		// Its transition matrices must be computed using the "subroot" node's edge length
		// The subroot node's transition matrices need not be calculated 
		TipData & ndTD = *(nd->GetTipData());
		TreeNode * subroot = nd->GetLeftChild();
		PHYCAS_ASSERT(subroot);
		PHYCAS_ASSERT(!subroot->GetRightSib()); //@POL need to create a IsSubroot() member function for TreeNode
		calcTMatForSim(ndTD, subroot->GetEdgeLen());
		++nd;

		// Skip subroot node as its transition matrices are never used and thus do not need to be computed
		++nd;

		// Process the remaining nodes in the tree
		for (; nd != t->end(); ++nd)
			{
			if (nd->IsTip())
				{
				TipData & ndTD = *(nd->GetTipData());
				calcTMatForSim(ndTD, nd->GetEdgeLen());
				}
			else
				{
				InternalData & ndID = *(nd->GetInternalData());
				calcPMat(ndID.getPMatrices(), nd->GetEdgeLen());
				}
			}
		}

	// Create a vector of cumulative state frequencies to use in choosing starting states
	const std::vector<double> & freqs = model->getStateFreqs();
	std::vector<double> cum_freqs(num_states, 0.0);
	std::partial_sum(freqs.begin(), freqs.end(), cum_freqs.begin());

	// Create a vector of cumulative rate probabilities to use in choosing relative rates
	std::vector<double> cum_rate_probs(num_rates, 0.0);
	std::partial_sum(rate_probs.begin(), rate_probs.end(), cum_rate_probs.begin());

	sim_data->resetPatternLength(t->GetNTips());
	sim_data->wipePattern();

	for (unsigned character = 0; character < nchar; ++character)
		{
		// Choose a rate for this character (actually, choose index, the actual rate is rate_means[r])
		unsigned r = 0;
		if (num_rates > 1)
			{
			// warning: removing the if statement will invalidate all examples involving simulated data with rate
			// homogeneity because of the call to rng->Uniform here!
			r = (unsigned)(std::lower_bound(cum_rate_probs.begin(), cum_rate_probs.end(), rng->Uniform(FILE_AND_LINE)) - cum_rate_probs.begin());
			}

		// Generate the starting state
		int8_t j = (unsigned)(std::lower_bound(cum_freqs.begin(), cum_freqs.end(), rng->Uniform(FILE_AND_LINE)) - cum_freqs.begin());

		// Assign starting state to the tip node currently serving as the root of the tree
		preorder_iterator nd = t->begin();
		TipData & rootTD = *(nd->GetTipData());
		rootTD.state = j;

		sim_data->setState(nd->GetNodeNumber(), j);

		// Go ahead and generate the state for the (only) descendant of the root node (the "subroot" node)
		// Note that the root node's T matrix is used for this calculation; the P matrix of the subroot node
		// is never computed
		unsigned parent_state = (unsigned)j;

		// Get the T matrix for the tip node serving as the root
		double * * Tmatrix = rootTD.pMatrixTranspose[r];

		// Choose a uniform random deviate
		double u = rng->Uniform(FILE_AND_LINE);

		// Spin the roulette wheel to choose a state for the subroot node
		double cum = 0.0;
		unsigned i = 0;
		for (; i < num_states; ++i)
			{
			double pr = Tmatrix[i][parent_state];
			//std::cerr << str(boost::format("Tmatrix[%d][%d] = %f") % i % parent_state % pr) << std::endl;
			cum += pr;
			if (u < cum)
				break;
			}

		// Increment iterator so that nd now refers to the subroot (sole descendant of the root)
		++nd;

		// Assign the new state to the subroot node
		InternalData & ndID = *(nd->GetInternalData());
		ndID.state = (int8_t)i;
		//std::cerr << "  Assigning state " << i << " to node " << nd->GetNodeNumber() << std::endl;

		// Walk the remainder of the tree using the preorder sequence, generating data for each node along the way
		for (++nd; nd != t->end(); ++nd)
			{
			// Get state of parent of nd
			TreeNode * parent = nd->GetParent();
			parent_state = UINT_MAX;
			if (parent->IsTip())
				{
				TipData * parentTD = parent->GetTipData();
				parent_state = (unsigned)parentTD->state;
				}
			else
				{
				InternalData * parentID = parent->GetInternalData();
				parent_state = (unsigned)parentID->state;
				}
			PHYCAS_ASSERT(parent_state < num_states);

			if (nd->IsTip())
				{
				// Get the T matrix
				TipData & ndTD = *(nd->GetTipData());
				double * * Tmatrix = ndTD.pMatrixTranspose[r];

				// Choose a uniform random deviate
				double u = rng->Uniform(FILE_AND_LINE);

				// Spin the roulette wheel and assign a state to nd
				double cum = 0.0;
				unsigned i = 0;
				for (; i < num_states; ++i)
					{
					double pr = Tmatrix[i][parent_state];
					//std::cerr << str(boost::format("Tmatrix[%d][%d] = %f") % i % parent_state % pr) << std::endl;
					cum += pr;
					if (u < cum)
						break;
					}
				ndTD.state = (int8_t)i;
				sim_data->setState(nd->GetNodeNumber(), (int8_t)i);
				}
			else
				{
				// Get the T matrix
				InternalData & ndID = *(nd->GetInternalData());
				double * * Pmatrix = ndID.pMatrices[r];

				// Choose a uniform random deviate
				double u = rng->Uniform(FILE_AND_LINE);

				// Spin the roulette wheel and assign a state to nd
				double cum = 0.0;
				unsigned i = 0;
				for (; i < num_states; ++i)
					{
					double pr = Pmatrix[parent_state][i];
					//std::cerr << str(boost::format("Pmatrix[%d][%d] = %f") % parent_state % i % pr) << std::endl;
					cum += pr;
					if (u < cum)
						break;
					}
				ndID.state = (int8_t)i;
				}
			}

		// We are now finished simulating data for one character, so insert the pattern just generated
		// into the pattern map maintained by sim_data; the 1.0 means that the count for this pattern
		// should be incremented by 1
		sim_data->insertPattern(1.0);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Stores all CondLikelihood objects stored in TipData or InternalData structures. This will require all CLAs to be
|	recomputed the next time the likelihood needs to be computed. Returns the subroot node (only child of root node).
*/
TreeNode * TreeLikelihood::storeAllCLAs(
  TreeShPtr t)				/**< is the tree from which all CLAs are to be removed */
	{
	// Start at the root node
	TreeNode * nd = t->GetFirstPreorder();
	PHYCAS_ASSERT(nd);

	// Move to the subroot node
	nd = nd->GetNextPreorder();
	PHYCAS_ASSERT(nd);

	// Invalidate (and do not cache) all CLAs from the tree. This will require all CLAs to be recomputed 
	// when the likelihood is computed using the subroot as the likelihood root. This path should be taken
	// if a parameter is changed that invalidates the entire tree.
	NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::invalidateBothEndsDiscardCache, this, _1, _2);
	effective_postorder_edge_iterator(nd, validFunctor); // constructor does all the work we need
	invalidateBothEndsDiscardCache(nd);

	return nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sanity check to make sure TreeLikelihood::storeAllCLAs really has stored all CLAs. Returns true if any CLAs are 
|	found, false if no CLAs are found (not even in cache positions).
*/
bool TreeLikelihood::debugCheckCLAsRemainInTree(
  TreeShPtr t) const	/**< is the tree to check */
	{
	for (TreeNode * nd = t->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			{
			TipData * td = nd->GetTipData();
			if (td != NULL)
				{
				if (td->parWorkingCLA)
					return true;
				if (td->parCachedCLA)
					return true;
				}
			}
		else
			{
			InternalData * id = nd->GetInternalData();
			if (id != NULL)
				{
				if (id->parWorkingCLA)
					return true;
				if (id->parCachedCLA)
					return true;

				if (id->childWorkingCLA)
					return true;
				if (id->childCachedCLA)
					return true;
				}
			}
		}	// preorder loop over all nodes

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Scans tree for cached CLAs. If any are found, returns true. If no CLAs are currently cached, returns false. False is
|	the expected result, because there should only be cached CLAs in the tree if we are in the middle of performing a 
|	Metropolis-Hastings move. If a move has been proposed, but not yet accepted or rejected, cached CLAs provide a way
|	to return to the pre-move state after rejection without having to recalculate the likelihood.
*/
bool TreeLikelihood::debugCheckForUncachedCLAs(
  TreeShPtr t)	 /**< is the tree to check */
  const
	{
	bool cached_CLAs_found = false;
	for (TreeNode * nd = t->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			{
			TipData * td = nd->GetTipData();
			if (td != NULL)
				{
				if (td->parCachedCLA)
					{
					cached_CLAs_found = true;
					break;
					}
				}
			}
		else
			{
			InternalData * id = nd->GetInternalData();
			if (id != NULL)
				{
				if (id->parCachedCLA)
					{
					cached_CLAs_found = true;
					break;
					}

				if (id->childCachedCLA)
					{
					cached_CLAs_found = true;
					break;
					}
				}
			}
		}	// preorder loop over all nodes
	return cached_CLAs_found;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Scans node for CLAs. If any are found, returns true. If no CLAs are found, returns false.
*/
bool TreeLikelihood::debugCheckCLAsRemainInNode(
  TreeNode * nd)   /**< is the node to check */
  const
	{
	bool CLAs_found = false;
	if (nd->IsTip())
		{
		TipData * td = nd->GetTipData();
		if (td != NULL)
			{
			if (td->parWorkingCLA || td->parCachedCLA)
				CLAs_found = true;
			}
		}
	else
		{
		InternalData * id = nd->GetInternalData();
		if (id != NULL)
			{
			if (id->parWorkingCLA || id->parCachedCLA || id->childWorkingCLA || id->childCachedCLA)
				CLAs_found = true;
			}
		}
	return CLAs_found;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Computes the value of the log-likelihood at substitutional saturation (this occurs when all edges have infinite
|   length). The site log-likelihood in this case is simply the sum of the logs of the equilibrium frequency of each tip
|   state: i.e.,
|>
|   log-like for site j = \sum_{i=1}^{ntax} log(\pi_{s_{ij}})
|>
|   Meaningless if no data has been stored, so assumes that `pattern_map' is not empty.
*/
double TreeLikelihood::calcLogLikeAtSubstitutionSaturation() const
    {
    //@POL needs to take ambiguity into account!
    PHYCAS_ASSERT(!no_data);
    PHYCAS_ASSERT(!pattern_map.empty());

    // Create a vector containing the log of each state frequency
    unsigned nstates = getNStates();
	const double * stateFreq = &model->getStateFreqs()[0]; //PELIGROSO
    std::vector<double> logfreq(nstates, 0.0);
    for (unsigned i = 0; i < nstates; ++i)
        logfreq[i] = log(stateFreq[i]);

    double lnL = 0.0;
    for (PatternMapType::const_iterator pit = pattern_map.begin(); pit != pattern_map.end(); ++pit)
        {
        const VecStateList & pattern = pit->first;
        const PatternCountType & count = pit->second;
        double site_log_like = 0.0;
        for (VecStateList::const_iterator sit = pattern.begin(); sit != pattern.end(); ++sit)
            {
            int8_t s = (*sit);
            PHYCAS_ASSERT(s >= (int8_t)0);
            PHYCAS_ASSERT(s < (int8_t)nstates);
            site_log_like += logfreq[s];
            }
        lnL += count*site_log_like;
        }
    return lnL;
    }

int calcLnLLevel = 0;
/*----------------------------------------------------------------------------------------------------------------------
|	This is the function that needs to be called to recompute the log-likelihood. If `likelihood_root' is not NULL, 
|	the calculation will use that node as the root of the likelihood calculation, and it will be assumed that all
|	conditional likelihood arrays (CLAs) are correctly calculated or have been invalidated if they need to be 
|	recomputed. On the other hand, if `likelihood_root' is NULL, then all CLAs will be invalidated and recomputed
|	(useful if a parameter has been changed that requires recalculation of all CLAs).
*/
double TreeLikelihood::calcLnL(
  TreeShPtr t)
	{
	if (no_data)
		return 0.0;

	// Compute likelihood using likelihood_root if specified
	// Assume that if likelihood_root has been specified, then the necessary 
	// CLA invalidations have already been performed.
	TreeNode * nd = likelihood_root;
	if (nd == NULL)
		{
#if 0
		// If no likelihood_root has been specified, use the subroot node (and
		// invalidate the entire tree to be safe)
		nd = t->GetFirstPreorder();
		PHYCAS_ASSERT(nd);

		// Move to the subroot node
		nd = nd->GetNextPreorder();
		PHYCAS_ASSERT(nd);

		// The subroot node is the new likelihood_root
		likelihood_root = nd;

		// Invalidate (and do not cache) all CLAs from the tree. This will require all CLAs to be recomputed 
		// when the likelihood is computed using the subroot as the likelihood root. This path should be taken
		// if a parameter is changed that invalidates the entire tree.
		NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::invalidateBothEndsDiscardCache, this, _1, _2);
		effective_postorder_edge_iterator(nd, validFunctor); // constructor does all the work we need
		invalidateBothEndsDiscardCache(nd);
#else
		// If no likelihood_root has been specified, invalidate the entire tree to be safe
		nd = storeAllCLAs(t);

		// The subroot node will be the new likelihood_root
		likelihood_root = nd;
#endif
		}

	PHYCAS_ASSERT(nd);
	PHYCAS_ASSERT(nd->IsInternal());

	// Uncomment line below to force recalculation of all CLAs
	//storeAllCLAs(t);

	if (using_unimap)
		{
		PHYCAS_ASSERT(sMatValid);
		return univentProbMgr.calcUnimapLnL(*t, num_patterns, &obs_state_counts[0], sMat);
		}

	// Calculate log-likelihood using nd as the likelihood root
	double lnL = calcLnLFromNode(*nd);
#	if 0 // !defined(NDEBUG)	
		if (calcLnLLevel == 0)
			{
			calcLnLLevel = 1;
			storeAllCLAs(t);
			double lnLRecalc = calcLnL(t);
			PHYCAS_ASSERT(fabs(lnL-lnLRecalc) < 0.000001);
			calcLnLLevel = 0;
			}
#	endif

    //startTreeViewer(t, "lnL = %.5f" % lnL);

	return lnL;
	}

void TreeLikelihood::debugSaveCLAs(TreeShPtr t, std::string fn, bool overwrite)
	{
	std::ofstream tmpf;
	if (overwrite)
		tmpf.open(fn.c_str());
	else
		{
		tmpf.open(fn.c_str(), std::ios::out | std::ios::app);
		}

	for (TreeNode * nd = t->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			{
			tmpf << "\n\nTip node number " << nd->GetNodeNumber();
			if (nd->IsTipRoot())
				tmpf << "\n	 parent = None";
			else
				tmpf << "\n	 parent = " << nd->GetParent()->GetNodeNumber();
			tmpf << "\n	 children = ";
			for (TreeNode * child = nd->GetLeftChild(); child != NULL; child = child->GetRightSib())
				tmpf << child->GetNodeNumber() << " ";

			TipData * td = nd->GetTipData();
			if (td->parentalCLAValid())
				{
				CondLikelihoodShPtr cla = td->getParentalCondLikePtr();
				unsigned sz = cla->getCLASize();
				LikeFltType * arr = cla->getCLA();
				tmpf << "\n	 cla length = " << sz;
				tmpf << "\n	 parental = ";
				for (unsigned i = 0; i < sz; ++i)
					{
					tmpf << str(boost::format("%20.12f ") % arr[i]);
					}
				}
			else
				{
				tmpf << "\n	 parental CLA = None";
				}
			}
		else
			{
			tmpf << "\n\nInternal node number " << nd->GetNodeNumber();

			tmpf << "\n	 parent = " << nd->GetParent()->GetNodeNumber();
			tmpf << "\n	 children = ";
			for (TreeNode * child = nd->GetLeftChild(); child != NULL; child = child->GetRightSib())
				tmpf << child->GetNodeNumber() << " ";

			InternalData * id = nd->GetInternalData();
			if (id->parentalCLAValid())
				{
				CondLikelihoodShPtr cla = id->getParentalCondLikePtr();
				unsigned sz = cla->getCLASize();
				LikeFltType * arr = cla->getCLA();
				tmpf << "\n	 parental CLA length = " << sz;
				tmpf << "\n	 parental CLA = ";
				for (unsigned i = 0; i < sz; ++i)
					{
					tmpf << str(boost::format("%20.12f ") % arr[i]);
					}
				}
			else
				{
				tmpf << "\n	 parental CLA = None";
				}

			if (id->filialCLAValid())
				{
				CondLikelihoodShPtr cla = id->getChildCondLikePtr();
				unsigned sz = cla->getCLASize();
				LikeFltType * arr = cla->getCLA();
				tmpf << "\n	 filial CLA length = " << sz;
				tmpf << "\n	 filial CLA = ";
				for (unsigned i = 0; i < sz; ++i)
					{
					tmpf << str(boost::format("%20.12f ") % arr[i]);
					}
				}
			else
				{
				tmpf << "\n	 filial CLA = None";
				}
			}
		}

	tmpf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the log-likelihood using `focal_node' as the likelihood root (i.e. the root of the likelihood calculation,
|	but not necessarily the root of the tree).
*/
double TreeLikelihood::calcLnLFromNode(
  TreeNode & focal_node)	/**< is the likelihood root (i.e. node around which the likelihood will be computed) */
	{
	double lnL;
	if (no_data)
		lnL =  0.0;
	else
		{
		PHYCAS_ASSERT(!focal_node.IsTip());

		// valid_functor will return true if the conditional likelihood arrays pointing away from the
		// focal_node are up-to-date, false if they need to be recomputed
		NodeValidityChecker valid_functor = boost::bind(&TreeLikelihood::isValid, this, _1, _2);

		// iter will visit nodes that need their CLAs updated centripetally (like a postorder traversal 
		// but also coming from below the focal node). Each node visited is guaranteed by valid_functor 
		// to need its CLA updated.
		effective_postorder_edge_iterator iter(&focal_node, valid_functor);
		effective_postorder_edge_iterator iter_end;
		for (; iter != iter_end; ++iter)
			{
			refreshCLA(*iter->first, iter->second);
			}

		// We have now brought all neighboring CLAs up-to-date, so we can now call harvestLnL to
		// compute the likelihood
		EdgeEndpoints edge(&focal_node, NULL);
		lnL = harvestLnL(edge);
		}
	++nevals;
	return lnL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates the TipData data structure needed to store the data for one tip (the tip corresponding to the supplied
|	`row' index in the data matrix `mat'). Returns a pointer to the newly-created TipData structure. See documentation 
|	for the TipData structure for more explanation.
*/
TipData * TreeLikelihood::allocateTipData(	//POLBM TreeLikelihood::allocateTipData
  unsigned row)		/**< is the row of the data matrix corresponding to the data for this tip node */
	{
	std::map<int8_t, int8_t>					globalToLocal;
	std::vector<unsigned int>					stateListVec;
	std::map<int8_t, int8_t>::const_iterator	foundElement;

	int8_t *									tipSpecificStateCode	= new int8_t[num_patterns];
	//@POL 21-Nov-2005 make tipSpecificStateCode a shared_array or a std:Vector - currently I don't think these are being deleted

	const int8_t								ns						= num_states;
	const int8_t								nsPlusOne				= num_states + 1;
	unsigned									nPartialAmbig			= 0;

	// Loop through all patterns for this row of the matrix. For each global state code encountered,
	// determine which local state code it represents, and build up the tipSpecificStateCode array
	// as we go.
	if (using_unimap)
		{
		unsigned inc_site = 0;
		const unsigned n_char = charIndexToPatternIndex.size();
		for (unsigned site = 0; site < n_char; ++site)
			{
			unsigned pattern_index = charIndexToPatternIndex[site];
			if (pattern_index != UINT_MAX)
				{
				PatternMapType::const_iterator it = pattern_map.begin();
				std::advance(it, pattern_index); //@POL not very efficient for map types, so may be worth using a vector rather than a map for pattern_map when using_unimap is true
				const int8_t globalStateCode = (it->first)[row];

				if (globalStateCode < nsPlusOne)
					{
					// no partial ambiguity, but may be gap state
					tipSpecificStateCode[inc_site] = (globalStateCode < 0 ? ns : globalStateCode);
					}
				else
					{
					// partial ambiguity
					foundElement = globalToLocal.find(globalStateCode);
					if (foundElement == globalToLocal.end())
						{
						// state code needs to be added to map
						globalToLocal[globalStateCode] = nPartialAmbig + nsPlusOne;
						stateListVec.push_back(state_list_pos[globalStateCode]);
						tipSpecificStateCode[inc_site] = nPartialAmbig + nsPlusOne;
						nPartialAmbig++;
						}
					else
						{
						// state code is already in the map
						tipSpecificStateCode[inc_site] = foundElement->second;
						}
					}
				++inc_site;
				}
			//std::cerr << row << ": " << site << " -> " << (int)(tipSpecificStateCode[site]) << std::endl;
			}
		PHYCAS_ASSERT(inc_site == num_patterns);
		}
	else    // not unimap
		{
		unsigned i = 0;
		for (PatternMapType::const_iterator it = pattern_map.begin(); it != pattern_map.end(); ++it, ++i)
			{
			const int8_t globalStateCode = (it->first)[row];

			if (globalStateCode < nsPlusOne)
				{
				// no partial ambiguity, but may be gap state
				tipSpecificStateCode[i] = (globalStateCode < 0 ? ns : globalStateCode);
				}
			else
				{
				// partial ambiguity
				foundElement = globalToLocal.find(globalStateCode);
				if (foundElement == globalToLocal.end())
					{
					// state code needs to be added to map
					globalToLocal[globalStateCode] = nPartialAmbig + nsPlusOne;
					stateListVec.push_back(state_list_pos[globalStateCode]);
					tipSpecificStateCode[i] = nPartialAmbig + nsPlusOne;
					nPartialAmbig++;
					}
				else
					{
					// state code is already in the map
					tipSpecificStateCode[i] = foundElement->second;
					}
				}
			}
		}
	return new TipData( using_unimap,
						num_patterns,
						stateListVec,												// stateListPosVec
						boost::shared_array<const int8_t>(tipSpecificStateCode),	// stateCodesShPtr
						num_rates,													// number of relative rate categories
						num_states,													// number of states in the model
						NULL,
						true,														// managePMatrices
						cla_pool);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates the InternalData data structure needed to store the conditional likelihood arrays and transition matrices
|	needed for likelihood calcuations. Returns a pointer to a newly-constructed InternalData structure. See the 
|	documentation for InternalData for more explanation.
*/
InternalData * TreeLikelihood::allocateInternalData()
	{
	return new InternalData(using_unimap,
							num_patterns,				// number of site patterns
							num_rates,					// number of relative rate categories
							num_states,					// number of model states
							NULL,						// pMat
							true,						// managePMatrices
							cla_pool);					
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function that deletes the TipData structure allocated by the allocateTipData() function. The intention is for this
|	function to be given to TreeNode objects (in the form of a boost::function object) whenever allocateTipData() is 
|	used to allocate memory for a TipData structure. The TreeNode destructor can then use the boost::function object to
|	delete the memory required by the TipData structure without knowing any details of the TipData structure (i.e.,
|	the file that defines the TreeNode destructor does not need to include a header file that declares the details of
|	the TipData class.
*/
void deallocateTipData(TipData * p)
	{
	delete p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function that deletes the InternalData structure allocated by the allocateInternalData() function. The intention is 
|	for this function to be given to TreeNode objects (in the form of a boost::function object) whenever 
|	allocateInternalData() is used to allocate memory for an InternalData structure. The TreeNode destructor can then 
|	use the boost::function object to delete the memory required by the InternalData structure without knowing any 
|	details of the InternalData structure (i.e., the file that defines the TreeNode destructor does not need to include 
|	a header file that declares the details of the InternalData class.
*/
void deallocateInternalData(InternalData * p)
	{
	delete p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a basic TipData object for all tip nodes and calls allocateInternalData() for all internal nodes in the 
|	supplied Tree. This function does not require a data matrix and does not add observed data to TipData structures.
|	Silently returns if root node already has a TipData structure. This is taken to indicate that prepareForLikelihood
|	or prepareForSimulation was previously called for this tree, in which case all data structures needed for 
|	simulation are already present.
*/
void TreeLikelihood::prepareForSimulation(
  TreeShPtr t)			/**< is the tree to decorate */
	{
	TreeNode::TipDataDeleter		td_deleter	= &deallocateTipData;
	TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;

	preorder_iterator nd = t->begin();

	// Only proceed if root node does not already have a TipData structure
	if (nd->GetTipData() != NULL)
		return;

	for (; nd != t->end(); ++nd)
		{
		if (nd->IsTip())
			{
			TipData * td =	new TipData(num_rates,	num_states, cla_pool);	//@POL should be using shared_ptr here?
			nd->SetTipData(td, td_deleter);
			}
		else
			{
			InternalData * cl = allocateInternalData();
			nd->SetInternalData(cl, cl_deleter);
			}
		}
	}

// BOOKMARK state_list and state_list_pos
/*----------------------------------------------------------------------------------------------------------------------
|	Builds `pattern_map' and `pattern_counts' using uncompressed data stored in `mat'.
*/
void TreeLikelihood::copyDataFromDiscreteMatrix(
  const NxsCXXDiscreteMatrix & mat)		/**< is the data source */
	{
	nTaxa = mat.getNTax();

    // These assignments should be kept before compressDataMatrix because they are used there
	state_list = mat.getStateList(); 
	state_list_pos = mat.getStateListPos();

    // The compressDataMatrix function first erases, then builds, both pattern_map and 
	// pattern_counts using the uncompressed data contained in mat
	num_patterns = compressDataMatrix(mat);

    buildConstantStatesVector();

	// size of likelihood_rate_site vector needs to be revisited if the number of rates subsequently changes 
	recalcRelativeRates();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Builds `pattern_map' and `pattern_counts' using data stored in `sim_data'.
*/
void TreeLikelihood::copyDataFromSimData(
  SimDataShPtr sim_data)	/**< is the data source */
	{
	// Copy simulated data to pattern_map
	pattern_map = sim_data->getSimPatternMap();

	// Build up counts vector
	pattern_counts.clear();
	for (PatternMapType::iterator it = pattern_map.begin(); it != pattern_map.end(); ++it)
		{
		pattern_counts.push_back(it->second);
		}

	nTaxa = sim_data->getPatternLength();
	num_patterns = (unsigned)pattern_map.size();

	model->buildStateList(state_list, state_list_pos);

	// size of likelihood_rate_site vector needs to be revisited if the number of rates subsequently changes 
	recalcRelativeRates();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds data currently stored in `pattern_map' to the patterns already in `other'. Assumes that `pattern_length' for 
|	this SimData object is identical to the `pattern_length' of `other'.
*/
void TreeLikelihood::addDataTo(SimData & other)
	{
	if (pattern_map.empty())
		return;
	if (other.getTotalCount() == 0)
		{
		PHYCAS_ASSERT(nTaxa > 0);
		other.resetPatternLength(nTaxa);
		}
	PHYCAS_ASSERT(nTaxa == other.getPatternLength());
	for (PatternMapType::iterator it = pattern_map.begin(); it != pattern_map.end(); ++it)
		{
		PatternCountType count = it->second;
		VecStateList & other_pattern = other.getCurrPattern();
		std::copy(it->first.begin(), it->first.end(), other_pattern.begin());
		other.insertPattern(count);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls allocateInternalData() to add an InternalData structure to `nd' containing the conditional likelihood arrays
|	needed for likelihood calculations.
*/
void TreeLikelihood::prepareInternalNodeForLikelihood(
  TreeNode * nd)	/**< is the node to decorate */
	{
#if POLPY_NEWWAY
	if (nd)
		{
    	InternalData * ndID = nd->GetInternalData();
        if (ndID == NULL)
            {
    		TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;
	    	InternalData * cl = allocateInternalData();
		    nd->SetInternalData(cl, cl_deleter);
            }
		}
#else
	if (nd)
		{
		TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;
		InternalData * cl = allocateInternalData();
		nd->SetInternalData(cl, cl_deleter);
		}
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Decorates an interior node with conditional likelihood and transition probability structures and add to the tree's
|	StoreInternalNode vector.
*/
void TreeLikelihood::addDecoratedInternalNode(
  TreeShPtr t,		/**< is the tree to decorate */
  unsigned num)		/**< is the number to assign to this node */
	{
	TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;
	InternalData * cl = allocateInternalData();
	TreeNode * nd = t->AllocNewNode();
	nd->SetNodeNum(num);
	nd->SetInternalData(cl, cl_deleter);
	t->StoreInternalNode(nd);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Decorates a tip node with data and add to the tree's tipStorage vector.
*/
void TreeLikelihood::addOrphanTip(
  TreeShPtr t,		/**< is the tree to decorate */
  unsigned row,		/**< is the row in the data matrix to associate with the added tip */
  std::string name) /**< is the taxon name */
	{
	TreeNode::TipDataDeleter td_deleter = &deallocateTipData;
	TipData * td = allocateTipData(row);
	TreeNode * nd = t->AllocNewNode();
	nd->SetTipData(td, td_deleter);
	nd->SetNodeNum(row);
	nd->SetNodeName(name);
	t->StoreLeafNode(nd);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls allocateTipData() for all tip nodes and allocateInternalData() for all internal nodes in the supplied Tree. 
|	Assumes each tip node number in the tree equals the appropriate row in the data matrix. 
*/
void TreeLikelihood::prepareForLikelihood( //POL_BOOKMARK TreeLikelihood::prepareForLikelihood
  TreeShPtr t)									/**< is the tree to decorate */
	{
	// If no_data is true, it means that calcLnL will always return 0.0 immediately and 
	// will thus never need the TipData or InternalData data structures
	if (no_data)
		return;

	TreeNode::TipDataDeleter		td_deleter	= &deallocateTipData;
	TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;

	// Put all existing conditional likelihood arrays already back into storage
	//storeAllCLAs(t);
#if POLPY_NEWWAY	//CLAShPtr
	//cla_pool->clearStack();	//@POL this here only because prepareForLikelihood called in NCatMove::proposeNewState when ncat is increased
#else
	//cla_pool.clearStack();	//@POL this here only because prepareForLikelihood called in NCatMove::proposeNewState when ncat is increased
#endif

	for (preorder_iterator nd = t->begin(); nd != t->end(); ++nd)
		{
		if (nd->IsTip())
			{
			unsigned row = nd->GetNodeNumber();
			TipData * td = allocateTipData(row);
			nd->SetTipData(td, td_deleter);
			}
		else
			{
			InternalData * cl = allocateInternalData();
			nd->SetInternalData(cl, cl_deleter);
			}
		}
	useAsLikelihoodRoot(NULL);

	//@POL	should decorate any nodes in tipStorage and nodeStorage now, but should wait until those 
	// are changed to vectors
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns reference to the vector data member `all_missing', which is a list of indices of sites for which all taxa
|   show completely missing data. The first site has index 0 in the `all_missing' vector.
*/
const std::vector<unsigned> & TreeLikelihood::getListOfAllMissingSites() const
    {
    return all_missing;
    }

#if 0
/*----------------------------------------------------------------------------------------------------------------------
|   Returns a list of all sites with patterns comprising only two primary states and no missing data or ambiguities for
|   any taxon.
*/
std::vector<unsigned> TreeLikelihood::findDataBipartitions() const
    {
    // needs to be written
    }
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Builds up the vector data member `constant_states' based on the patterns in `pattern_map'. The `constant_states' 
|   vector has a similar structure to the `state_list' vector. Here is an example `constant_states' vector for 4 sites:
|>
|   +---+---+---+---+---+---+---+---+
|   | 0 | 1 | 0 | 1 | 3 | 2 | 0 | 2 |
|   +---+---+---+---+---+---+---+---+
|     ^   ^       ^       ^
|     |   |       |       |
|     |   |       |       4th site can potentially be constant for either A (state 0) or G (state 2)
|     |   |       3rd site can be potentially constant for T (state 3)
|     |   2nd site can be potentially constant for 1 state: that state (A, state 0) follows in the next cell
|     1st site is definitely variable (hence the 0 meaning no states follow)
|>
|   The function returns the number of potentially constant patterns found (note the return value is NOT the number
|   of potentially constant sites).
*/
unsigned TreeLikelihood::buildConstantStatesVector()
	{
    PHYCAS_ASSERT(!pattern_map.empty());
    constant_states.clear();

    unsigned num_potentially_constant = 0;
    std::set<int> common_states;    // holds running intersection across taxa for this pattern
    std::set<int> curr_states;      // holds states for taxon under consideration
    std::vector<int> xset;          // temporarily holds intersection of common_states and curr_states

    for (PatternMapType::const_iterator pat = pattern_map.begin(); pat != pattern_map.end(); ++pat)
        {
        bool site_potentially_constant = true;
        for (unsigned taxon = 0; taxon < nTaxa; ++taxon)
            {
            if (!site_potentially_constant)                
                break;

            // Reminder of how data members state_list and state_list_pos are laid out:
            //                                         ?         N      {CGT}  {ACG}    R,(AG)    Y
            // states          | A | C | G | T |  - A C G T | A C G T | C G T | A G T |  A  G  | C T
            // state_list      1 0 1 1 1 2 1 3 5 -1 0 1 2 3 4 0 1 2 3 3 1 2 3 3 0 2 3 2  0  2  4 1 3
            // state_list_pos  0   2   4   6   8           14        19      23      27       30

            int code = (int)(pat->first)[taxon];
            PHYCAS_ASSERT(code >= 0);   // Mark, why don't ? and - states trigger this assert? Does this have to do with the fact that we have abandoned the Cipres version of NCL? We translate - to ? in the NxsCXXDiscreteMatrix constructor
            unsigned pos = (unsigned)state_list_pos[code];
            unsigned n = (unsigned)state_list[pos];
            curr_states.clear();

            // Insert all states for the current taxon into the curr_states set
            for (unsigned x = pos + 1; x < pos + n + 1; ++x)
                {
                int c = (int)state_list[x];
                PHYCAS_ASSERT(c >= 0);
                if (site_potentially_constant)
                    curr_states.insert(c);
                }

            if (taxon == 0)
                {
                // For first taxon, let common_states equal curr_states
                common_states.clear();
                common_states.insert(curr_states.begin(), curr_states.end());
                }
            else
                {
                // For subsequent taxa, common_states should contain only those states possessed
                // by all taxa (i.e. it will be the empty set if the site is variable)
                xset.clear();
                std::set_intersection(common_states.begin(), common_states.end(), curr_states.begin(), curr_states.end(), back_inserter(xset));
                if (xset.empty())
                    site_potentially_constant = false;
                else
                    {
                    common_states.clear();
                    common_states.insert(xset.begin(), xset.end());
                    }
                }
            }

        if (site_potentially_constant)
            {
            ++num_potentially_constant;
            constant_states.push_back((unsigned)common_states.size());
            for (std::set<int>::const_iterator it = common_states.begin(); it != common_states.end(); ++it)
                {
                constant_states.push_back(*it);
                }
            }
        else
            {
            constant_states.push_back(0);
            }
        }
    return num_potentially_constant;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Copies data from `mat' to the map `pattern_map'. The resulting map holds pairs whose key is the pattern for one site
|	and whose value is a count of the number of sites having that pattern. The counts from `pattern_map' are transferred  
|	to the `pattern_counts' vector (vectors are more efficient containers for use during likelihood calculations).
*/
unsigned TreeLikelihood::compressDataMatrix(const NxsCXXDiscreteMatrix & mat) //POL_BOOKMARK TreeLikelihood::compressDataMatrix
	{
	pattern_map.clear();
	unsigned ntax = mat.getNTax();
	unsigned nchar = mat.getNChar();
	charIndexToPatternIndex.assign(nchar, UINT_MAX);

    // Create actingWeights vector and copy the integer weights from mat into it
    // If there are no integer weights in mat, copy the floating point weights instead
    // if floating point weights have been defined
	const std::vector<int> & iwts = mat.getIntWeightsConst();
	std::vector<double> actingWeights(nchar, 1.0);
	if (!iwts.empty())
		{
		PHYCAS_ASSERT(iwts.size() >= nchar);
		for (unsigned j = 0; j < nchar; ++j)
			actingWeights[j] = (double)iwts.at(j);
		}
	else
		{
		const std::vector<double> & dwts = mat.getDblWeightsConst();
		if (!dwts.empty())
			{
			actingWeights = dwts;
			PHYCAS_ASSERT(actingWeights.size() >= nchar);
			}
		}

    // Set corresponding actingWeights elements to zero if any characters have been excluded in mat
	const std::set<unsigned> & excl = mat.getExcludedCharIndices();
	for (std::set<unsigned>::const_iterator eIt = excl.begin(); eIt != excl.end(); ++eIt)
		{
		PHYCAS_ASSERT(*eIt < nchar);
		actingWeights[*eIt] = 0.0;
		}
	const double * wts = &(actingWeights[0]);

	// patternToIndex is a map that associates a list of character indices with each pattern. Thus, if 
	// some pattern is found at sites 0, 15, and 167, then patternToIndex.first is the pattern and 
	// patternToIndex.second is the list<unsigned> [0, 15, 167]
	typedef std::list<unsigned> IndexList;
	typedef std::map<VecStateList, IndexList> PatternToIndex;
	PatternToIndex patternToIndex;

	if (model->isCodonModel()) //@POL could move this test inside the taxon-loop to avoid being so redundant
		{
		// Loop across each triplet of sites in mat
		for (unsigned j = 0; j < nchar; j += 3)
			{
			// Suppose nchar=5 and j=3: not enough sites to make a second codon. In this example,
			// j+2 = 5, which equals nchar, so break out of loop over sites
			if (j+2 >= nchar)
				break;

			// here we (arbitrarily) use the max weight of any char in the codon
			PatternCountType charWt = (wts ? std::max(wts[j], std::max(wts[j+1], wts[j+2])) : 1.0); 
			if (charWt > 0.0)
				{
				// Build up a vector representing the pattern of state codes at this site
				std::vector<int8_t> pattern;
				for (unsigned i = 0; i < ntax; ++i)
					{
					const int8_t * row	= mat.getRow(i);
					const int8_t   code1 = row[j];
					const int8_t   code2 = row[j+1];
					const int8_t   code3 = row[j+2];
					bool code1_ok = (code1 >= 0 && code1 < 4);
					bool code2_ok = (code2 >= 0 && code2 < 4);
					bool code3_ok = (code3 >= 0 && code3 < 4);
					if (code1_ok && code2_ok && code3_ok)
						{
						const int8_t code = codon_state_codes[16*code1 + 4*code2 + code3]; // CGT = 27 = 16*1 + 4*2 + 3
						if (code > 60)
							throw XLikelihood(str(boost::format("Stop codon encountered for taxon %d at sites %d-%d") % (i+1) % (j+1) % (j+4)));
						pattern.push_back(code);
						}
					else
						pattern.push_back((int8_t)61);
					}
				
				//@POL below here same as the else block
	
				// Add the pattern to the map if it has not yet been seen, otherwise increment 
				// the count of this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
				PatternMapType::iterator lowb = pattern_map.lower_bound(pattern);
				if (lowb != pattern_map.end() && !(pattern_map.key_comp()(pattern, lowb->first)))
					{
					// pattern is already in pattern_map, increment count
					lowb->second += charWt;
					}
				else
					{
					// pattern has not yet been stored in pattern_map
					pattern_map.insert(lowb, PatternMapType::value_type(pattern, charWt));
					}
				
				// Add the pattern to the map if it has not yet been seen, otherwise increment 
				// the count of this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
				PatternToIndex::iterator pToILowB = patternToIndex.lower_bound(pattern);
				if (pToILowB != patternToIndex.end() && !(patternToIndex.key_comp()(pattern, pToILowB->first)))
					pToILowB->second.push_back(j);
				else
					{
					IndexList ilist(1, j);	// create a list<unsigned> containing 1 element whose value is j
					patternToIndex.insert(pToILowB, PatternToIndex::value_type(pattern, ilist));
					}
				}
			}
		}
	else // not codon model
		{
		std::vector<int8_t> pattern;
        unsigned nStates = getNStates();

		// Loop across each site in mat
		for (unsigned j = 0; j < nchar; ++j)
			{
			PatternCountType charWt = (wts ? wts[j] : 1.0); 
			if (charWt > 0.0)
				{
				// Build up a vector representing the pattern of state codes at this site
				pattern.clear();
                unsigned num_all_missing = 0;
				for (unsigned i = 0; i < ntax; ++i)
					{
					const int8_t * row	= mat.getRow(i);
					const int8_t   code = row[j];
					pattern.push_back(code);
                    if ((unsigned) code == nStates)	
                        ++num_all_missing;
					}

                // Do not include the pattern if it contains only completely missing data
                // for all taxa
                if (num_all_missing == ntax)
                    {
                    all_missing.push_back(j);
                    continue;
                    }

				// Add the pattern to the map if it has not yet been seen, otherwise increment 
				// the count of this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
				PatternMapType::iterator lowb = pattern_map.lower_bound(pattern);
				if (lowb != pattern_map.end() && !(pattern_map.key_comp()(pattern, lowb->first)))
					{
					// pattern is already in pattern_map, increment count
					lowb->second += charWt;
					}
				else
					{
					// pattern has not yet been stored in pattern_map
					pattern_map.insert(lowb, PatternMapType::value_type(pattern, charWt));
					}
				
				// Add the pattern to the patternToIndex map if not already present, and then
				// append the character index to the list of character indices associated with
				// the pattern
				PatternToIndex::iterator pToILowB = patternToIndex.lower_bound(pattern);
				if (pToILowB != patternToIndex.end() && !(patternToIndex.key_comp()(pattern, pToILowB->first)))
					pToILowB->second.push_back(j);
				else
					{
					IndexList ilist(1, j);	// create a list<unsigned> containing 1 element whose value is j
					patternToIndex.insert(pToILowB, PatternToIndex::value_type(pattern, ilist));
					}
				}
			}
		}   // if model->isCodonModel() ... else ...

	// Copy counts to pattern_counts before returning
	pattern_counts.clear();
	pattern_counts.reserve(pattern_map.size());
	unsigned patternIndex = 0;
	unsigned n_inc_chars = 0;
	for (PatternMapType::iterator mapit = pattern_map.begin(); mapit != pattern_map.end(); ++mapit, ++patternIndex)
		{
		pattern_counts.push_back(mapit->second);

		// Find pattern in patternToIndex, which provides a list of indices of sites having that pattern
		// For each site index in the list, add an element to the map charIndexToPatternIndex
		// Now, charIndexToPatternIndex[i] points to the index in pattern_map for the pattern found at site i
		const IndexList & inds = patternToIndex[mapit->first];
		for (IndexList::const_iterator indIt = inds.begin(); indIt != inds.end(); ++indIt)
			{
			charIndexToPatternIndex[*indIt] = patternIndex;
			++n_inc_chars;
			}
		}

	if (using_unimap)
		return n_inc_chars;
	else
		return (unsigned)pattern_map.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a string representation of the supplied `state'. For example, if `state' equals 1 (where model is a 
|	standard DNA model), the string returned would be "C". If, however, `state' was 7 (again, standard DNA model), then
|	there is some ambiguity present and the string returned might look something like "{AC}" (the string returned 
|	depends of course on the actual meaning of the global state code 7).
*/
std::string TreeLikelihood::getStateStr(
  int8_t state) const /**< is the global state code to be converted to a std::string */
	{
	std::string s;
	const int8_t nsPlusOne = num_states + 1;

	if (state < nsPlusOne)
		{
		// either no ambiguity or complete ambiguity
		s << model->lookupStateRepr((int)state);
		}
	else
		{
		// `state' represents partial ambiguity

		// First, find location of the definition of `state' in the global state list
		unsigned pos = state_list_pos[(unsigned)state];
		VecStateList::const_iterator it = state_list.begin() + pos;

		// Now get the number of basic states composing `state'
		unsigned n = *it++;

		// Walk down global state list converting states into strings
		s << "{";
		for (unsigned i = 0; i < n; ++i)
			{
			s << model->lookupStateRepr((int)*it++);
			}
		s << "}";
		}
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets private data member `nevals' to 0.
*/
void TreeLikelihood::resetNEvals()
	{
	nevals = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor providing access to the current value of the private data member `nevals'.
*/
unsigned TreeLikelihood::getNEvals()
	{
	return nevals;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assuming TreeLikelihood::compressDataMatrix has been called, so that `pattern_map' is up-to-date, returns a string 
|	listing all observed patterns and their frequencies.
*/
std::string TreeLikelihood::listPatterns(
  bool show_coded_states)	/**< if true, output global state codes used internally; otherwise, do the translation back to the original state representations */
	{
	std::vector<std::string> state_repr;
	std::string s;
	unsigned i = 0;
	PatternMapType::iterator it = pattern_map.begin();
	for (; it != pattern_map.end(); ++it, ++i)
		{
		const std::vector<int8_t> & p = it->first;
		PatternCountType c = it->second;
		s << str(boost::format("%6d %6.1f ") % i % c);
		unsigned ntax = (unsigned)p.size();
		if (show_coded_states)
			{
			for (unsigned j = 0; j < ntax; ++j)
				{
				s << str(boost::format("%d ") % (int)(p[j]));
				}
			}
		else
			{
			for (unsigned j = 0; j < ntax; ++j)
				{
				s << str(boost::format("%s") % getStateStr(p[j]));
				}
			}
		s << '\n';
		}
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First erases `pattern_repr', then builds it into a vector of pattern representation strings. When the function 
|	returns, each element of `pattern_repr' is a string containing global state codes separated by whitespace. This 
|	function is intended to be used for debugging purposes, where it is sometimes helpful to be sure of exactly which
|	data pattern is being operated upon. Thus, no particular effort has been made to make this function efficient.
*/
void TreeLikelihood::buildPatternReprVector(std::vector<std::string> & pattern_repr, TreeShPtr t)
	{
	unsigned nTips		= t->GetNTips();
	//unsigned nStates	= getNStates();
	//unsigned nPatterns	= getNPatterns();

	pattern_repr.clear();
		pattern_repr.reserve(num_patterns);

	std::cerr << "\nglobal_pos:" << std::endl;
	for (unsigned z = 0; z < state_list_pos.size(); ++z)
		{
		std::cerr << str(boost::format("%3d") % z) << "	 " << (int)state_list_pos[z] << std::endl;
		}

	// Recreate the data matrix in terms of tip-specific state codes
	int8_t * * m = NewTwoDArray<int8_t>(nTips, num_patterns);

	for (preorder_iterator nd = t->begin(); nd != t->end(); ++nd)
		{
		if (nd->IsTip())
			{
			// get node number i
			unsigned i = nd->GetNodeNumber();
			PHYCAS_ASSERT(i >= 0 && i < nTips);

			// get tip-specific state code array
			TipData * tipData = nd->GetTipData();
			PHYCAS_ASSERT(tipData);
			const int8_t * tipCodes = tipData->getConstStateCodes();

			// get position vector that allows us to translate tip-specific state codes into global state codes
			const std::vector<unsigned int> & global_statelist_pos = tipData->getConstStateListPos();

			std::cerr << "\nglobal_statelist_pos vector for node " << i << ": " << std::endl;
			for (unsigned z = 0; z < global_statelist_pos.size(); ++z)
				{
				std::cerr << "	" << (int)global_statelist_pos[z] << std::endl;
				}

			// fill row i of data matrix
			for (unsigned j = 0; j < num_patterns; ++j)
				{
				int8_t global_code = tipCodes[j];
				int8_t offset = global_code - ((int8_t)num_states + 1);
				if (offset >= 0)
					{
					unsigned pos = global_statelist_pos[offset];
					for (unsigned m = 0; m < (unsigned)state_list_pos.size(); ++m)
						{
						if (state_list_pos[m] == pos)
							global_code = (int8_t)(m);
						}
					}
				m[i][j] = global_code;

				//std::cerr << j << " | " << (int)tipCodes[j] << " | ";
				//if (offset >= 0)
				//	std::cerr << (int)global_statelist_pos[offset];
				//else
				//	std::cerr << "(" << (int)offset << ")";
				//std::cerr << std::endl;
				}
			}
		}

	for (unsigned j = 0; j < num_patterns; ++j)
		{
		std::string s;
		for (unsigned i = 0; i < nTips; ++i)
			{
			s << (int)m[i][j] << " ";
			}
		pattern_repr.push_back(s);
		}

	DeleteTwoDArray<int8_t>(m);
	}

}	// namespace phycas
