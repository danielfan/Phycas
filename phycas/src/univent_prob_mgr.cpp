/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2008 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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
#include <boost/bind.hpp>
#include "phycas/src/univent_prob_mgr.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/basic_tree_node.hpp"
#include "phycas/src/edge_iterators.hpp"
#include "phycas/src/tip_data.hpp"
#include "phycas/src/internal_data.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor
*/
UniventProbMgr::UniventProbMgr(
  ModelShPtr modelArg)  /**< is a shared pointer to the model object */
  : lnUMatMemMgt(modelArg->getNStates(), 0.0)
  , numStates(modelArg->getNStates())
  , model(modelArg)
  , sampleTimes(false)
  , scratchMatOne(modelArg->getNStates(), 0.0)
  , scratchMatTwo(modelArg->getNStates(), 0.0)
  , storeUnivents(true), isMappingValidVar(false)
	{
	lnUMat = lnUMatMemMgt.GetMatrix();
	
	unsigned init_max_m = 5;
	maxm = init_max_m;

	// build up logmfact
	// m   ->  m!          log(m!)
	// 0   ->  1   = 1        0.0 =   0.0
	// 1   ->  1   = 1*1      0.0 =   0.0 +   0.0
	// 2   ->  2   = 1*2    0.693 =   0.0 + 0.693
	// 3   ->  6   = 2*3    1.792 = 0.693 + 1.099
	// 4   -> 24   = 6*4    3.178 = 1.792 + 1.386
	logmfact.resize(init_max_m + 1, 0.0);
	logmfact[0] = 0.0;
	for (unsigned m = 1; m <= init_max_m; ++m)
		logmfact[m] = logmfact[m-1] + log((double)m);
	}

/*
double UniventProbMgr::GetLambda() 
	{
	if (model->timestamp() != this->modelTimestamp)
		recalcUMat();
	return this->lambda;
	}
*/
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UniventProbMgr::recalcUMat()
    {
	this->lambda = model->calcLMat(lnUMat);
	uMatVect.resize(2);
    unsigned prev_maxm = maxm;
	for (unsigned m = 0; m <= 2; ++m)
        {
        if (uMatVect[m].GetDimension() == 0)
            uMatVect[m].CreateMatrix(numStates, 0.0);
        if (m == 0)
            uMatVect[m].Identity();
        else if (m == 1)
            {
            // uniformized transition matrix
            double * * oneUPtr = uMatVect[m].GetMatrix();
            lambda = model->calcUMat(oneUPtr);
            }
		}
	maxm = 1;

	 // the reduceMaxm bit here is a hack to try to reduce maxm as opposed to allowing it to continue to creep up.
	const unsigned reduceMaxm = 2;
	expandUMatVect(prev_maxm > ( reduceMaxm + 2) ? prev_maxm - reduceMaxm : prev_maxm);
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Chooses a value of m, the number of univents on a particular edge for a particular site. 
*/
unsigned UniventProbMgr::sampleM(
  int8_t start_state,               /**< is the state at the beginning of the edge */
  int8_t end_state,                 /**< is the state at the end of the edge */
  double transition_prob,           /**< is the probability of `end_state' given `start_state' (marginalized over all possible numbers of univents) */
  double edgelen,                   /**< is the length of the edge in expected number of substitutions per site */
  Lot & rng) const /**< is the random number generator to use for the mapping */
    {
    unsigned m = UINT_MAX;
    std::vector<double> probm;
    std::vector<double> cumprm;

    double u = rng.Uniform(FILE_AND_LINE);

    // ok will be set to true if an m value can be sampled. The reason m might not be sampled is
    // because maxm might be not set high enough, in which case maxm will be doubled and another
    // attempt to sample m will be made.
    bool ok = false;
    while (!ok)
        {
        probm.clear();

        // First draw m, the number of univents on this edge
        double lambda_t       = edgelen*lambda;
        double log_lambda_t   = log(lambda_t);
        double log_trans_prob = log(transition_prob);
        for (unsigned z = 0; z <= maxm; ++z)
            {
            double logprm = (double)z*log_lambda_t - lambda_t - logmfact[z] - log_trans_prob;
            double pij = uMatVect[z][start_state][end_state]; //@POL should uMatVect hold L matrices rather than U matrices?
            probm.push_back(pij*exp(logprm));
            }
        cumprm.resize(probm.size());

        // partial_sum adds successive elements of prob together and stores in cumprm
        std::partial_sum(probm.begin(), probm.end(), cumprm.begin());

        // lower_bound returns an iterator positioned at the first element in cumprm with a value greater than or equal
        // to a uniform random deviate. 
        std::vector<double>::iterator it = std::lower_bound(cumprm.begin(), cumprm.end(), u);
        if (it != cumprm.end())
            {
            ok = true;
            m = (unsigned)(it - cumprm.begin());
            }
        else
            expandUMatVect(maxm*2);
        }
    // Subtract cumprm.begin() from it to yield the sampled value of m
    PHYCAS_ASSERT(m != UINT_MAX);
    return m;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UniventProbMgr::expandUMatVect(unsigned m) const 
	{
	if (m <= maxm)
		return;
		
	PHYCAS_ASSERT(maxm >= 1);
	uMatVect.resize(m + 2);
	// we could get better cache efficiency by storing the transpose of one uMatVect
	// whenever we calculate it in 
	const double * const * onePtr = const_cast<const double * const * >(uMatVect[1].GetMatrix());
	for (unsigned k = maxm + 1; k <= m + 1; ++k)
		{
		
		if (uMatVect[k].GetDimension() == 0)
			uMatVect[k].CreateMatrix(numStates, 0.0);
		double * * currPtr = uMatVect[k].GetMatrix();
		const double * const * prevPtr = const_cast<const double * const * >(uMatVect[k-1].GetMatrix());
		for (unsigned i = 0; i < numStates; ++i)
			{
			double * currRow = currPtr[i];
			const double * prevRow = prevPtr[i];
			for (unsigned j = 0; j < numStates; ++j)
				{
				double sum = 0.0;
				for (unsigned k = 0; k < numStates; ++k)
					sum += prevRow[k]*onePtr[k][j];
				currRow[j] = sum;
				}
			}
		}

	// The m that was sampled was greater than maxm, indicating that 
	// maxm is not large enough, so double maxm and try again
	const unsigned prevlogmfactsize = logmfact.size();
	if (m > prevlogmfactsize)
		{
		// extend logmfact vector
		logmfact.resize(m + 1, 0.0);
		for (unsigned m = prevlogmfactsize; m <= maxm; ++m)
			logmfact[m] = logmfact[m - 1] + log((double)m);
		std::cerr << "\nIncreasing maxm to " << maxm << std::endl;
		for (unsigned i = 0; i <= maxm; ++i)
			std::cerr << '\t' << i << '\t' << logmfact[i] << std::endl;
		}
	maxm = m;

	}

/*----------------------------------------------------------------------------------------------------------------------
|   Refreshes the uniformized mapping for one site on one particular edge. 
*/
void UniventProbMgr::unimapEdgeOneSite(
  Univents & u,  /**< site-specific vector to hold state-time pairs representing the mapping on this edge */
  unsigned pattern_index,
  int8_t start_state,               /**< is the state at the beginning of the edge */
  int8_t end_state,                 /**< is the state at the end of the edge */
  double transition_prob,           /**< is the probability of `end_state' given `start_state' (marginalized over all possible numbers of univents) */
  double edgelen,                   /**< is the length of the edge in expected number of substitutions per site */
  bool doSampleTimes,
  Lot & rng) const /**< is the random number generator to use for the mapping */
	{
    unsigned m = sampleM(start_state, end_state, transition_prob, edgelen, rng);
    u.mdot += m;
	if (doSampleTimes)
		{
		std::vector<double> & t = u.times.at(pattern_index);
		t.clear();
		t.resize(m);
		for (unsigned k = 0; k < m; ++k)
            t[k] = (float)rng.Uniform(FILE_AND_LINE);
		std::sort(t.begin(), t.end());
		}
	StateMapping & sm = u.univents.at(pattern_index);
	sm.clear();
	sm.resize(m);

    // Now sample the m states
    if (m == 0)
        {
        PHYCAS_ASSERT(start_state == end_state);
        }
    else if (m == 1)
        sm[0] = end_state;
    else
        {
#		if defined(SIMULATE_MAPPINGS_GIVEN_M)
			// Now draw m states until a sample is obtained in which the mth state equals end_state
			double u;
			bool done = false;
			unsigned ntries = 0;
			while (!done)
				{
				ntries++;
				int8_t prev_state = start_state;
				for (k = 0; k < m; ++k)
					{
					const double * uMat_row = uMat[prev_state];
					u = rng.Uniform(FILE_AND_LINE);
					double cump = 0.0;
					int8_t new_state;
					bool found = false;
					for (unsigned j = 0; j < numStates; ++j)
						{
						cump += exp(uMat_row[j])/lambda;
						if (u < cump)
							{
							new_state = j;
							found = true;
							break;
							}
						}
					PHYCAS_ASSERT(found);
					sm[k] = new_state;
					prev_state = new_state;
					}
				if (ntries > 1000 || sm[m-1] == end_state)
					done = true;
				}
			PHYCAS_ASSERT(ntries <= 1000);
			
#		else //SIMULATE_MAPPINGS_GIVEN_M

			const double * const * one_umat = getUMatConst(1);
			const double * const * curr_umat = getUMatConst(m);
			int8_t prev_state = start_state;
			for (unsigned curr_m = 0; curr_m < m - 1; ++curr_m)
				{
				const double * const * rest_umat = getUMatConst(m - curr_m - 1);
				const double prob_all_mappings = curr_umat[prev_state][end_state];
				double u = rng.Uniform(FILE_AND_LINE)*prob_all_mappings;
				unsigned s = 0;
				/*TEMP should just loop over single-mutation neighbors */
				for (; s < numStates; ++s)
					{
					u -= one_umat[prev_state][s]* rest_umat[s][end_state];
					if (u < 0.0)
						break;
					}
                if (s >= numStates && u >= 1.0e-7)
                    {
                    std::cerr << "woah!" << std::endl;
                    }
				PHYCAS_ASSERT(s < numStates || u < 1.0e-7);
				sm[curr_m] = s;
				curr_umat = rest_umat;
                prev_state = s;
				}
			sm[m-1] = end_state;
#		endif //SIMULATE_MAPPINGS_GIVEN_M
		}
    }

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
/*==============================================================================
| uses scratchMatTwo (and indirectly scratchMatOne)
*/
void UniventProbMgr::sampleUniventsKeepEndStates(Univents & u, const double edgelen, const int8_t * par_states, const double * * p_mat_transposed, Lot &rng) const
	{
	const int8_t * des_states = &(u.end_states_vec[0]);
	double ** p_mat = scratchMatTwo.GetMatrix();
	fillTranspose(p_mat, p_mat_transposed, numStates);
	sampleUnivents(u, edgelen, par_states, des_states, const_cast<const double * const *>(p_mat), rng, NULL);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
/*==============================================================================
| uses scratchMatOne
*/
void UniventProbMgr::sampleUnivents(Univents & u, const double edgelen,  const int8_t * par_states, const int8_t * des_states, const double * const * p_mat, Lot & rng, unsigned ** s_mat) const
	{
	const unsigned num_patterns = u.univents.size();
	u.mdot = 0;
	scratchMatOne.Fill(1.0);
	const bool doSampleTimes = this->sampleTimes;
    const double lambda_t       = edgelen*lambda;
    const double log_lambda_t   = log(lambda_t);
	std::vector<double> elogprmVec; // holds factors that do not depend on the start and end states
	for (unsigned pattern_index = 0; pattern_index < num_patterns; ++pattern_index)
		{
		const int8_t start_state = *par_states++;
	    const int8_t end_state = *des_states++;
	    // begin UniventProbMgr::unimapEdgeOneSite inlined manually 
	    const double trans_prob = p_mat[start_state][end_state];
		/* begin sampleM inlined manually.... */
		unsigned m = UINT_MAX;
		double uni_variate = rng.Uniform(FILE_AND_LINE);
			// ok will be set to true if an m value can be sampled. The reason m might not be sampled is
			// because maxm might be not set high enough, in which case maxm will be doubled and another
			// attempt to sample m will be made.
		double total_prob = 0.0;
		unsigned next_z = 0;
		//std::cerr << "Sampling " << (int)start_state << " -> " << (int)end_state << " which has trans_prob = " << trans_prob << " uni_variate=" <<uni_variate;		
		for (;;)
			{
			for (unsigned z = next_z; z <= maxm; ++z)
				{
				if (z >= elogprmVec.size())
					{
					elogprmVec.reserve(maxm);
					elogprmVec.push_back(exp((double)z*log_lambda_t - lambda_t - logmfact[z]));
					PHYCAS_ASSERT(z < elogprmVec.size());
					}
				const double elogprm = elogprmVec[z];
				double pij = uMatVect[z][start_state][end_state]; //@POL should uMatVect hold L matrices rather than U matrices?
				//std::cerr << "(elogprm=" << elogprm << ", pij=" << pij << ") ";
				const double pr = pij*elogprm/trans_prob;
				total_prob += pr;
				if (total_prob > uni_variate)
					{
					m = z;
					break;
					}
				}
			if (m == UINT_MAX)
				{
				next_z = maxm + 1;
				expandUMatVect(maxm*2);
				}
			else
				break;
			}
		/* end  sampleM inlined manually */
    	//std::cerr << " m = " << m << '\n';
		u.mdot += m;
		if (doSampleTimes)
			{
			std::vector<double> & t = u.times.at(pattern_index);
			t.clear();
			t.resize(m);
			for (unsigned k = 0; k < m; ++k)
				t[k] = (float)rng.Uniform(FILE_AND_LINE);
			std::sort(t.begin(), t.end());
			}
		StateMapping * sm = 0L;
		if (storeUnivents)
			{
			sm = &(u.univents.at(pattern_index));
			sm->clear();
			sm->resize(m);
			}
		// Now sample the m states
		if (m == 0)
			{
			PHYCAS_ASSERT(start_state == end_state);
			}
		else if (m == 1)
			{
			if (s_mat)
				s_mat[start_state][end_state] += 1;		
			if (sm)
				(*sm)[0] = end_state;
			}
		else
			{
			const double * const * one_umat = getUMatConst(1);
			const double * const * curr_umat = getUMatConst(m);
			int8_t prev_state = start_state;
			for (unsigned curr_m = 0; curr_m < m - 1; ++curr_m)
				{
				const double * const * rest_umat = getUMatConst(m - curr_m - 1);
				const double prob_all_mappings = curr_umat[prev_state][end_state];
				double u = rng.Uniform(FILE_AND_LINE)*prob_all_mappings;
				unsigned s = 0;
				/*TEMP should just loop over single-mutation neighbors */
				for (; s < numStates; ++s)
					{
					u -= one_umat[prev_state][s]* rest_umat[s][end_state];
					if (u < 0.0)
						break;
					}
				if (s >= numStates && u >= 1.0e-7)
					{
					std::cerr << "woah!" << std::endl;
					}
				PHYCAS_ASSERT(s < numStates || u < 1.0e-7);
				if (s_mat)
					s_mat[start_state][s] += 1;
				if (sm)
					(*sm)[curr_m] = s;
				curr_umat = rest_umat;
				prev_state = s;
				}
			if (s_mat)
				s_mat[prev_state][end_state] += 1;
			if (sm)
				(*sm)[m-1] = end_state;
			}

	  	}
	u.setValid(true);
	u.times_valid = this->sampleTimes;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UniventProbMgr::sampleDescendantStates(
	const unsigned num_patterns, 
	int8_t * nd_states, 
	const double * const * p_mat, 
	const LikeFltType * des_cla, 
	const int8_t * parent_states,
	Lot & rng) const
	{
	std::vector<double> post_prob(numStates);
	
	for (unsigned i = 0 ; i < num_patterns; ++i)
		{
		const int8_t par_state = *parent_states++;
		double total = 0.0;
		for (unsigned j = 0; j < numStates; ++j)
			{
			post_prob[j] =  (*des_cla++)* p_mat[par_state][j];
			total += post_prob[j];
			}
		for (unsigned j = 0; j < numStates; ++j)
			post_prob[j] /= total;
		nd_states[i] = rng.MultinomialDraw(&post_prob[0], numStates);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UniventProbMgr::sampleRootStates(
	const unsigned num_patterns,
	int8_t * nd_states,
	const LikeFltType * rootStatePosterior,
	Lot &rng,
	bool posteriors_normalized, 
	unsigned * obs_state_counts) const
	{
	if (obs_state_counts)
		{
		for (unsigned i = 0 ; i < numStates; ++i)
			obs_state_counts[i] = 0;
		}
	for (unsigned i = 0 ; i < num_patterns; ++i)
		{
		double total = 1.0;
		if (!posteriors_normalized)
			total = std::accumulate(rootStatePosterior, rootStatePosterior + numStates, 0.0); 
		const int8_t st = rng.MultinomialDraw(rootStatePosterior, numStates, total);
		nd_states[i] = st;
		if (obs_state_counts)
			obs_state_counts[st] += 1;
		rootStatePosterior += numStates;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This is the function that is called (from calcLnL) to recompute the log-likelihood if `using_unimap' is true. It 
|	computes the log-likelihood using the stored matrix of uniformized mapping transition counts (`sMat') as well as the
|	stored `mdot' values in the TipData or InternalData structure associated with each node. This method assumes that
|	both `lambda' and the uniformized transition matrix `uMat' are up-to-date, and that the counts of observed states 
|	in the tip root (which always serves as the likelihood root when using_unimap is true) have been stored in
*/
double UniventProbMgr::calcUnimapLnL(
  const Tree &  t,
  const unsigned num_patterns,
  const unsigned * obs_state_counts,
  const unsigned * const * sMat)
	{
	PHYCAS_ASSERT(isMappingValidVar);
	lambda = model->calcLMat(lnUMat);
	double nsites = (double)num_patterns;

	// Compute term that is the product of relative frequencies of each state
	// sampled in the node serving as the likelihood root (the subroot node).
	double log_basal_freqs = 0.0;
	const std::vector<double> & freqs = model->getStateFreqs();
	for (unsigned i = 0; i < numStates; ++i)
		log_basal_freqs += (double)(obs_state_counts[i])*log(freqs[i]);

	// Compute term that is the product of all edge lengths each raised to the power
	// mdot, where mdot is the total number of univents over all sites along that edge. Also use 
	// this opportunity to calculate the tree length, which will be used later.
	double tree_length		   = 0.0;
	double log_edgelen_to_mdot = 0.0;
	const TreeNode * nd = t.GetFirstPreorderConst()->GetNextPreorderConst();
	for (; nd != NULL; nd = nd->GetNextPreorderConst())
		{
		double edgelen = nd->GetEdgeLen();
		tree_length += edgelen;
		const Univents &u =  getUniventsConstRef(*nd);
		log_edgelen_to_mdot += u.getMDot()*log(edgelen);
		}

	// Compute term that is the product of uniformized transition probabilities for all univents
	double log_uij_to_sij = 0.0;
	for (unsigned i = 0; i < numStates; ++i)
		{
		for (unsigned j = 0; j < numStates; ++j)
			{
			// If you are wondering why I haven't taken the log of the uMat entry, I should remind
			// you that it is the log of a transition probability that is stored in uMat[i][j]
			log_uij_to_sij += (double)(sMat[i][j])*lnUMat[i][j];
			//std::cerr << "sMat[i][j] = " << sMat[i][j] << '\n';
			}
		}

	double log_likelihood = log_basal_freqs + log_edgelen_to_mdot + log_uij_to_sij - nsites*lambda*tree_length;
	//std::cerr << "log_likelihood = log_basal_freqs + log_edgelen_to_mdot + log_uij_to_sij - nsites*lambda*tree_length;\n";
	//std::cerr << log_likelihood << " = "<< log_basal_freqs << " + " << log_edgelen_to_mdot << " + " << log_uij_to_sij << " - " << nsites<< " * " << lambda << " *" << tree_length << std::endl;
	return log_likelihood;
	}





} // namespace phycas

