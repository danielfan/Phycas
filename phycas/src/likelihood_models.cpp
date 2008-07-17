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

#include <cmath>
#include <iostream>
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/cipres/AllocateMatrix.hpp"
#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
#	include <boost/python/numeric.hpp>
#	include "phycas/src/thirdparty/num_util/num_util.h"
#endif
using std::cout;
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Recalculates all relative rate means and their probabilities, storing these in the `rates' and `probs' vectors
|	supplied as reference arguments. If using a discrete gamma ("G") model only, the number of elements in both `rates'
|	and `probs' equals the number of gamma rate categories. If using a pinvar model ("I") only, the number of elements 
|	in both vectors will be 2, with rates[0] = 0.0, rates[1] = 1.0/(1.0 - pinvar), probs[0] = pinvar and probs[1] =
|	1.0 - pinvar. If using an "I+G" model, there will be num_gamma_rates + 1 elements in both vectors, with rates[0] =
|	0.0 and probs[0] = pinvar, and the remaining elements containing the gamma rates and probabilities corrected for
|	the value of pinvar.
*/
void Model::recalcRatesAndProbs( //POL_BOOKMARK Model::recalcRatesAndProbs
  std::vector<double> & rates, 
  std::vector<double> & probs) const
	{
	if (is_flex_model)
		{
		normalizeRatesAndProbs(rates, probs);
		}
	else
		{
		std::vector<double> boundaries;

	    if (is_pinvar_model)
		    {
		    PHYCAS_ASSERT(pinvar < 1.0);
		    double one_minus_pinvar = 1.0 - pinvar;

#if POLPY_NEWWAY
		    // First calculate the gamma rates alone
		    std::vector<double> gamma_rates;
		    recalcGammaRatesAndBoundaries(gamma_rates, boundaries);

		    // Copy the rate probabilities directly
		    probs.resize(num_gamma_rates, 0.0);
		    std::copy(gamma_rate_probs.begin(), gamma_rate_probs.end(), probs.begin());

            //temporary!
            //std::cerr << "\nRate probability vector:" << std::endl;
            //for (unsigned i = 0; i < num_gamma_rates; ++i)
            //    {
            //    std::cerr << boost::str(boost::format("%6d %12.5f") % (i+1) % probs[i]) << std::endl;
            //    }

		    // Adjust gamma_rates using pinvar and save to rates vector
		    rates.resize(num_gamma_rates, 0.0);
		    std::vector<double>::iterator       rates_iter       = rates.begin();
		    std::vector<double>::const_iterator gamma_rates_iter = gamma_rates.begin();
		    for (unsigned i = 0; i < num_gamma_rates; ++i)
			    {
			    *rates_iter++ = (*gamma_rates_iter++)/one_minus_pinvar;
			    }

            //temporary!
            //std::cerr << "\nGamma rates vector before correction:" << std::endl;
            //for (unsigned i = 0; i < num_gamma_rates; ++i)
            //    {
            //    std::cerr << boost::str(boost::format("%6d %12.5f") % (i+1) % gamma_rates[i]) << std::endl;
            //    }

            //temporary!
            //std::cerr << "\nGamma rates vector after correction:" << std::endl;
            //for (unsigned i = 0; i < num_gamma_rates; ++i)
            //    {
            //    std::cerr << boost::str(boost::format("%6d %12.5f") % (i+1) % rates[i]) << std::endl;
            //    }
#else
		    // First calculate the gamma rates alone
		    std::vector<double> gamma_rates;
		    recalcGammaRatesAndBoundaries(gamma_rates, boundaries);

		    // Build up the rates and probs vectors using the local vector gamma_rates and the data members 
		    // pinvar and gamma_rate_probs
		    rates.resize(num_gamma_rates + 1, 0.0);
		    probs.resize(num_gamma_rates + 1, 0.0);
		    std::vector<double>::iterator rates_iter = rates.begin();
		    std::vector<double>::iterator probs_iter = probs.begin();
		    std::vector<double>::const_iterator gamma_probs_iter = gamma_rate_probs.begin();
		    std::vector<double>::const_iterator gamma_rates_iter = gamma_rates.begin();
		    (*rates_iter++) = 0.0;
		    (*probs_iter++) = pinvar;
		    for (unsigned i = 0; i < num_gamma_rates; ++i)
			    {
			    *rates_iter++ = (*gamma_rates_iter++)/one_minus_pinvar;
			    *probs_iter++ = (*gamma_probs_iter++)*one_minus_pinvar;
			    }
#endif
		    }
	    else
		    {
		    // The rates and probs computed by recalcGammaRatesAndBoundaries are what we need if is_pinvar_model is false
		    recalcGammaRatesAndBoundaries(rates, boundaries);
		    probs.resize(num_gamma_rates, 0.0);
		    std::copy(gamma_rate_probs.begin(), gamma_rate_probs.end(), probs.begin());
		    }
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Determines the category means of a discrete gamma distribution with `num_gamma_rates' categories. A continuous gamma
|	distribution is divided into `num_gamma_rates' sections, the boundaries of which are determined by the probability 
|	mass assumed to lie in each section. The vector `gamma_rate_probs' holds the probabilities corresponding to each 
|	category, which by default are all 1/`num_gamma_rates'. Upon return, rates[i] holds the mean rate for the category 
|	starting at boundaries[i]. For i < `num_gamma_rates' - 1, the upper limit of the rate category is at boundaries[i+1],
|	but for the last rate category (i.e. i = `num_gamma_rates' - 1), the upper boundary is infinity. Assumes 
|	`num_gamma_rates' is greater than zero.
|
|	Assume there are four rate categories, and let the means of these categories be labeled r_0, r_1, r_2 and r_3. Let 
|	the boundaries of the (e.g. four) rate categories be b_0, b_1, b_2, b_3 and b_4, with b_0 = 0.0 and b_4 = infinity.
|	The mean of rate category i is r_i and can be computed as the expected value of r in the interval [b_i, b_{i+1})
|	divided by the integral of the continuous gamma distribution over the same interval. For a gamma distribution having
|	shape alpha and scale beta, we have (in pseudo-LaTeX):
|>	
|	       / b_{i+1}   r^{alpha - 1} exp{-r/beta}
|	       |         r -------------------------- dr
|	       / b_i        beta^alpha  Gamma(alpha)      
|	r_i = ------------------------------------------
|	       / b_{i+1}   r^{alpha - 1} exp{-r/beta}   
|	       |           -------------------------- dr
|	       / b_i        beta^alpha  Gamma(alpha)
|
|	      (alpha)(beta)[CumGamma(b_{i+1}, alpha + 1, beta) - CumGamma(b_i, alpha + 1, beta)]
|	    = ----------------------------------------------------------------------------------
|	                 CumGamma(b_{i+1}, alpha, beta) - CumGamma(b_i, alpha, beta)
|>
*/
void Model::recalcGammaRatesAndBoundaries(std::vector<double> & rates, std::vector<double> & boundaries) const
	{
	PHYCAS_ASSERT(num_gamma_rates > 0);
	rates.resize(num_gamma_rates, 0.0);
	boundaries.resize(num_gamma_rates, 0.0);

	if (num_gamma_rates == 1)
		{
		rates[0] = 1.0;
		boundaries[0] = 0.0;
		return;
		}

	PHYCAS_ASSERT(gamma_shape > 0.0);
	double alpha = gamma_shape;
	double beta = 1.0/gamma_shape;
	
	//std::cerr << "\nGamma rate category information for shape = " << gamma_shape << std::endl;

	double cum_upper		= 0.0;
	double cum_upper_plus	= 0.0;
	double upper			= 0.0;
	double cum_prob			= 0.0;
	for (unsigned i = 1; i <= num_gamma_rates; ++i)
		{
		double lower				= upper;
		double cum_lower_plus		= cum_upper_plus;
		double cum_lower			= cum_upper;
		cum_prob					+= gamma_rate_probs[i-1];

		if (i < num_gamma_rates)
			{
			upper					= cdf.SampleGamma(cum_prob, alpha, beta);
			cum_upper_plus			= cdf.CumGamma(upper, alpha + 1.0, beta);
			cum_upper				= cdf.CumGamma(upper, alpha, beta);
			}
		else
			{
			cum_upper_plus			= 1.0;
			cum_upper				= 1.0;
			}

		double numer				= cum_upper_plus - cum_lower_plus;
		double denom				= cum_upper - cum_lower;
		double r_mean				= (denom > 0.0 ? (alpha*beta*numer/denom) : 0.0);
		rates[i-1]					= r_mean;
		boundaries[i-1]				= lower;

		//std::cerr << str(boost::format("%6d %12.6f %12.6f") % (i-1) % boundaries[i-1] % rates[i-1]) << std::endl;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `gamma_rates_unnorm[param_index]'.
*/
double Model::getFlexRateUnnorm(
  unsigned param_index)		/**< the 0-based index into the `gamma_rates_unnorm' vector of the element to return */
	{
	return gamma_rates_unnorm[param_index];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `gamma_rate_probs[param_index]'.
*/
double Model::getFlexProbUnnorm(
  unsigned param_index)		/**< the 0-based index into the `gamma_rate_probs' vector of the element to return */
	{
	return gamma_rate_probs[param_index];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `rel_rates[param_index]'.
*/
double GTR::getRelRateUnnorm(
  unsigned param_index)		/**< the 0-based index into the `gamma_rate_probs' vector of the element to return */
	{
	return rel_rates[param_index];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `state_freq_unnorm[param_index]'.
*/
double Model::getStateFreqUnnorm(
  unsigned param_index)		/**< the 0-based index into the `state_freq_unnorm' vector of the element to return */
	{
	return state_freq_unnorm[param_index];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The base class version of this function should be called by all derived classes because it is where the edge 
|	length parameters, the edge length hyperparameters and rate heterogeneity (gamma shape and pinvar) parameters are 
|	created.
*/
void Model::createParameters(
  TreeShPtr t,									    /**< is the tree (the nodes of which are needed for creating edge length parameters) */
  MCMCUpdaterVect & edgelens_vect_ref,			    /**< is the vector of edge length parameters to fill */
  MCMCUpdaterVect & edgelen_hyperparams_vect_ref,	/**< is the vector of edge length hyperparameters to fill */
  MCMCUpdaterVect & parameters_vect_ref) const	    /**< is the vector of model-specific parameters to fill */
    {
	PHYCAS_ASSERT(t);
	PHYCAS_ASSERT(edgelens_vect_ref.empty());
	PHYCAS_ASSERT(edgelen_hyperparams_vect_ref.empty());
	PHYCAS_ASSERT(parameters_vect_ref.empty());

	// Add the edge length parameter(s)
    if (separate_int_ext_edgelen_priors)
        {
	    // Add two edge length parameters to manage the priors for all edge lengths in the tree;
        // One of these two parameters will be in charge of the prior on internal edge lengths,
        // whereas the other will be in charge of external edge lengths. These edge length parameters 
        // do not actually update any edge lengths: a Metropolis proposal such as the LargetSimonMove 
        // is responsible for updating edge lengths.

        // First the parameter governing external edge length priors
        MCMCUpdaterShPtr p_external = MCMCUpdaterShPtr(new EdgeLenMasterParam(EdgeLenMasterParam::external));
	    std::string nm_external = str(boost::format("external edge length parameter"));
	    p_external->setName(nm_external);
	    p_external->setTree(t);
	    p_external->setPrior(externalEdgeLenPrior);
	    if (edge_lengths_fixed)
		    p_external->fixParameter();
	    edgelens_vect_ref.push_back(p_external);

        // Now the parameter governing internal edge length priors
	    MCMCUpdaterShPtr p_internal = MCMCUpdaterShPtr(new EdgeLenMasterParam(EdgeLenMasterParam::internal));
	    std::string nm_internal = str(boost::format("internal edge length parameter"));
	    p_internal->setName(nm_internal);
	    p_internal->setTree(t);
	    p_internal->setPrior(internalEdgeLenPrior);
	    if (edge_lengths_fixed)
		    p_internal->fixParameter();
	    edgelens_vect_ref.push_back(p_internal);
        }
    else
        {
	    // Add one edge length parameter to manage the priors for all edge lengths in the tree. 
        // This edge length parameter does not actually update any edge lengths: a Metropolis proposal 
        // such as the LargetSimonMove is responsible for updating edge lengths.
        MCMCUpdaterShPtr p = MCMCUpdaterShPtr(new EdgeLenMasterParam(EdgeLenMasterParam::both));
	    std::string nm_external = str(boost::format("master edge length parameter"));
	    p->setName(nm_external);
	    p->setTree(t);
	    p->setPrior(externalEdgeLenPrior);
	    if (edge_lengths_fixed)
		    p->fixParameter();
	    edgelens_vect_ref.push_back(p);
        }

	// Save a vector of shared pointers so that we can modify their fixed/free status if we need to
	edgelen_params.resize(edgelens_vect_ref.size());
	std::copy(edgelens_vect_ref.begin(), edgelens_vect_ref.end(), edgelen_params.begin());

	// Add edge length hyperparameters if requested
	if (edgeLenHyperPrior)
		{
        // For each edge length parameter, create an edge length hyperparameter, passing it a shared
        // pointer to the corresponding edge length parameter 
        for (MCMCUpdaterVect::iterator it = edgelen_params.begin(); it != edgelen_params.end(); ++it)
            {
            EdgeLenMasterParamShPtr pit = boost::shared_dynamic_cast<EdgeLenMasterParam>(*it);
		    MCMCUpdaterShPtr p = MCMCUpdaterShPtr(new HyperPriorParam(pit));
		    p->setName(std::string("edge length hyperparameter"));
		    p->setTree(t);
		    p->setPrior(edgeLenHyperPrior);
		    if (edgelen_hyperprior_fixed)
			    p->fixParameter();
		    edgelen_hyperparams_vect_ref.push_back(p);
            }

        // Save a vector of shared pointers so that we can modify their fixed/free status if we need to
	    edgelen_hyper_params.resize(edgelen_hyperparams_vect_ref.size());
	    std::copy(edgelen_hyperparams_vect_ref.begin(), edgelen_hyperparams_vect_ref.end(), edgelen_hyper_params.begin());
		}

	// Create any model-specific parameters and add to the parameters vector
	if (is_flex_model)
		{
		gamma_rates_unnorm.resize(num_gamma_rates, 0.0);
		PHYCAS_ASSERT(flex_rate_params.empty());
		PHYCAS_ASSERT(flex_prob_params.empty());
		for (unsigned i = 0; i < num_gamma_rates; ++i)
			{ //POL_BOOKMARK creating flex model parameters
			// start with rates drawn from Uniform(0.0, flex_upper_rate_bound)
			//@POL to do this right, need to draw from prior, but if number of spacers is small, drawing from
            // Uniform(0,flex_upper_rate_bound) will give almost the same results
			double u = flex_prob_param_prior->GetLot()->Uniform(FILE_AND_LINE);
			gamma_rates_unnorm[i] = flex_upper_rate_bound*u;
			//old way: gamma_rates_unnorm[i] = flex_upper_rate_bound*(double)(i + 1)/(double)(num_gamma_rates + 1);

			// start with probabilities drawn from the prior
			PHYCAS_ASSERT(flex_prob_param_prior);
			gamma_rate_probs[i] = flex_prob_param_prior->Sample();
			}

		// Rates must be sorted from lowest to highest to begin with
		std::sort(gamma_rates_unnorm.begin(), gamma_rates_unnorm.end());

		MCMCUpdaterShPtr rate_param = MCMCUpdaterShPtr(new FlexRateParam(num_flex_spacers, flex_upper_rate_bound, gamma_rates_unnorm));
		rate_param->setName("FLEX rates"); //@POL shouldn't this be done in the constructor?
		rate_param->setTree(t);
		rate_param->setPrior(flex_rate_param_prior);
		if (flex_rates_fixed)
			rate_param->fixParameter();
		parameters_vect_ref.push_back(rate_param);
		flex_rate_params.push_back(rate_param);

		MCMCUpdaterShPtr prob_param = MCMCUpdaterShPtr(new FlexProbParam(gamma_rate_probs));
		prob_param->setName("FLEX probs"); //@POL shouldn't this be done in the constructor?
		prob_param->setTree(t);
		prob_param->setPrior(flex_prob_param_prior);
		if (flex_probs_fixed)
			prob_param->fixParameter();
		parameters_vect_ref.push_back(prob_param);
		flex_prob_params.push_back(prob_param);
		}
	else if (num_gamma_rates > 1)
		{
		PHYCAS_ASSERT(num_gamma_rates > 1);
		PHYCAS_ASSERT(!gamma_shape_param);
		gamma_shape_param = MCMCUpdaterShPtr(new DiscreteGammaShapeParam(invert_shape));
        if (invert_shape)
    		gamma_shape_param->setName("Discrete gamma variance"); //@POL shouldn't this be done in the constructor?
        else
    		gamma_shape_param->setName("Discrete gamma shape"); //@POL shouldn't this be done in the constructor?
        gamma_shape_param->setStartingValue(gamma_shape);
		gamma_shape_param->setTree(t);
		gamma_shape_param->setPrior(gamma_shape_prior);
		if (gamma_shape_fixed)
			gamma_shape_param->fixParameter();
		parameters_vect_ref.push_back(gamma_shape_param);
		}

	if (is_pinvar_model)
		{
		PHYCAS_ASSERT(!pinvar_param);
		pinvar_param = MCMCUpdaterShPtr(new PinvarParam());
		pinvar_param->setName("Proportion of invariable sites");
        pinvar_param->setStartingValue(pinvar);
		pinvar_param->setTree(t);
		pinvar_param->setPrior(pinvar_prior);
		if (pinvar_fixed)
			pinvar_param->fixParameter();
		parameters_vect_ref.push_back(pinvar_param);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Stores a flattened version of the supplied 2-dimensional array `twoDarr', storing the result in the supplied VecDbl
|	reference variable `p'. The supplied `twoDarr' should be laid out so that rows occupy contiguous memory.
*/
inline void Model::flattenTwoDMatrix(VecDbl & p, double * * twoDarr, unsigned dim) const
	{
	unsigned flat_length = dim*dim;
	p.reserve(flat_length);
	double * twoD_begin = &twoDarr[0][0];
	double * twoD_end   = twoD_begin + flat_length;
	std::copy(twoD_begin, twoD_end, p.begin());
	}
	
#if defined(PYTHON_ONLY)
#if defined(USING_NUMARRAY)
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
boost::python::numeric::array Model::getPMatrix(double edgeLength) const
	{
	//@POL this function, along with calcPMat and calcPMatrices, should be provided by a policy class
	// (which would be QMatrix for GTR model). This function is full of ad hoc nonsense at the moment!
	double * * pMat = NewTwoDArray<double>(num_states, num_states);
	calcPMat(pMat, edgeLength);

	VecDbl p;
	flattenTwoDMatrix(p, pMat, num_states);

	DeleteTwoDArray<double>(pMat);

	// create vector of dimensions
	std::vector<int> dim_vect;
	dim_vect.push_back((int)num_states);
	dim_vect.push_back((int)num_states);

	return num_util::makeNum(&p[0], dim_vect);	//PELIGROSO
	}
#else
/*----------------------------------------------------------------------------------------------------------------------
|	This function was written for the purpose of displaying a transition probability matrix on the Python side. It is 
|	not intended to be used for anything other than display purposes (just use calcPMat directly if you want to 
|	efficiently calculate transition probabilities. It first creates a temporary 2-dimensional array, calls calcPMat 
|	to compute the transition probabilities and fill this array, then flattens out the matrix, returning a 1-dimensional
|	vector containing the first row, then second row, etc.
*/
VecDbl Model::getPMatrix(double edgeLength) const
	{
	//@POL this function, along with calcPMat and calcPMatrices, should be provided by a policy class
	// (which would be QMatrix for GTR model).
	double * * pMat = NewTwoDArray<double>(num_states, num_states);
	calcPMat(pMat, edgeLength);

	VecDbl p;
	flattenTwoDMatrix(p, pMat, num_states);

	DeleteTwoDArray<double>(pMat);

	return p;
	}
#endif
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the transition probability matrix given an edge length. Overrides the pure virtual function inherited from 
|	the base class Model. For the JC69 model, the transition probabilities are:
|>
|	Pii = 0.25 + 0.75*exp{-4*beta*t}
|	Pij = 0.25 - 0.25*exp{-4*beta*t}
|>
|	The evolutionary distance `edgeLength' = 3*beta*t, so as a function of edge length, the transition probabilities
|	are:
|>
|	Pii = 0.25 + 0.75*exp{-4*edgeLength/3}
|	Pij = 0.25 - 0.25*exp{-4*edgeLength/3}
|>
*/
void JC::calcPMat(double * * pMat, double edgeLength) const
	{
	const double exp_term = exp(-(4.0/3.0)*edgeLength);
	const double prob_change = 0.25 - 0.25*exp_term;
	const double prob_nochange = 0.25 + 0.75*exp_term;
	pMat[0][0] = prob_nochange;
	pMat[0][1] = prob_change;
	pMat[0][2] = prob_change;
	pMat[0][3] = prob_change;
	pMat[1][0] = prob_change;
	pMat[1][1] = prob_nochange;
	pMat[1][2] = prob_change;
	pMat[1][3] = prob_change;
	pMat[2][0] = prob_change;
	pMat[2][1] = prob_change;
	pMat[2][2] = prob_nochange;
	pMat[2][3] = prob_change;
	pMat[3][0] = prob_change;
	pMat[3][1] = prob_change;
	pMat[3][2] = prob_change;
	pMat[3][3] = prob_nochange;
	//for (unsigned i = 0 ; i < 4; ++i)
	//	cout << pMat[i][0] << ' '<< pMat[i][1] << ' '<< pMat[i][2] << ' '<< pMat[i][3] << '\n';
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function 
|   inherited from the base class Model. For the JC69 model, the uniformized transition probability matrix is:
|>
|   1 - 3*beta/lambda      beta/lambda         beta/lambda        beta/lambda
|     beta/lambda       1 - 3*beta/lambda      beta/lambda        beta/lambda
|     beta/lambda          beta/lambda      1 - 3*beta/lambda     beta/lambda
|     beta/lambda          beta/lambda         beta/lambda      1 - 3*beta/lambda
|>
|	where lambda is slightly larger than the largest element in the rate matrix (i.e. slightly larger than beta).
|   See Mateiu, L., and B. Rannala (2006. Systematic Biology 55:259-269) for details.
|   For efficiency, what is actually computed is the product of lambda and the above matrix, and what is stored are
|   the logarithms of the elements (because these values are always used on the log scale): i.e.
|>
|   log(lambda - 3*beta)      log(beta)             log(beta)             log(beta)
|       log(beta)         log(lambda - 3*beta)      log(beta)             log(beta)
|       log(beta)             log(beta)         log(lambda - 3*beta)      log(beta)
|       log(beta)             log(beta)             log(beta)         log(lambda - 3*beta)
|>
|   Because rate and time are confounded, beta is arbitrarily scaled so that time is measured in 
|   expected number of substitutions (i.e. set beta such that the edge length v = t):
|>  
|   v = 3*beta*t
|   v = 3*(1/3)*t
|>
|   Setting beta = 1/3 results in v = t.
*/
double JC::calcLMat(double * * lMat) const
	{
    // Mateiu and Rannala:   Me:
    //  B = 4/3              beta = 1/3
    //  Qii = -1             Qii = -3*beta = -1
    //  Qij = 1/3            Qij = beta    = 1/3
    //  Pii = 3/4            Pii = 1 - 3*beta/lambda = 1 - 1/4 = 3/4
    //  Pij = 1/12           Pij = beta/lambda = 1/12
    //  lambda = 4
    double beta = 1.0/3.0;
	const double lambda = 1.1;
	const double diag_term = log(lambda - 3.0*beta);
	const double off_diag_term = log(beta);
	lMat[0][0] = diag_term;
	lMat[0][1] = off_diag_term;
	lMat[0][2] = off_diag_term;
	lMat[0][3] = off_diag_term;
	lMat[1][0] = off_diag_term;
	lMat[1][1] = diag_term;
	lMat[1][2] = off_diag_term;
	lMat[1][3] = off_diag_term;
	lMat[2][0] = off_diag_term;
	lMat[2][1] = off_diag_term;
	lMat[2][2] = diag_term;
	lMat[2][3] = off_diag_term;
	lMat[3][0] = off_diag_term;
	lMat[3][1] = off_diag_term;
	lMat[3][2] = off_diag_term;
	lMat[3][3] = diag_term;
    return lambda;
	}

double JC::calcUMat(double * * uMat) const
	{
    double beta = 1.0/3.0;
	const double lambda = 1.1;
	const double diag_term = 1.0 - 3.0*beta/lambda;
	const double off_diag_term = beta/lambda;
	uMat[0][0] = diag_term;
	uMat[0][1] = off_diag_term;
	uMat[0][2] = off_diag_term;
	uMat[0][3] = off_diag_term;
	uMat[1][0] = off_diag_term;
	uMat[1][1] = diag_term;
	uMat[1][2] = off_diag_term;
	uMat[1][3] = off_diag_term;
	uMat[2][0] = off_diag_term;
	uMat[2][1] = off_diag_term;
	uMat[2][2] = diag_term;
	uMat[2][3] = off_diag_term;
	uMat[3][0] = off_diag_term;
	uMat[3][1] = off_diag_term;
	uMat[3][2] = off_diag_term;
	uMat[3][3] = diag_term;
    return lambda;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function 
|   inherited from the base class Model. 
*/
double HKY::calcLMat(double * * lMat) const
	{
	double piA = state_freqs[0];
	double piC = state_freqs[1];
	double piG = state_freqs[2];
	double piT = state_freqs[3];
	double piR = piA + piG;
	double piY = piC + piT;

    // set beta such that edgelen = t
    double beta = 1.0/(2.0*piR*piY + 2.0*kappa*(piA*piG + piC*piT));
        
    // set lambda to maximim substitution rate + epsilon
    double lambda_epsilon = 0.1;
    double max_rate = beta*(piY + piG*kappa);
    double next_rate = beta*(piR + piT*kappa);
    if (next_rate > max_rate)
        max_rate = next_rate;
    next_rate = beta*(piY + piA*kappa);
    if (next_rate > max_rate)
        max_rate = next_rate;
    next_rate = beta*(piR + piC*kappa);
    if (next_rate > max_rate)
        max_rate = next_rate;
    double lambda = max_rate + lambda_epsilon;

    lMat[0][0] = log(lambda - beta*(piY + piG*kappa));
    lMat[0][1] = log(piC*beta);
    lMat[0][2] = log(piG*kappa*beta);
    lMat[0][3] = log(piT*beta);
    lMat[1][0] = log(piA*beta);
    lMat[1][1] = log(lambda - beta*(piR + piT*kappa));
    lMat[1][2] = log(piG*beta);
    lMat[1][3] = log(piT*kappa*beta);
    lMat[2][0] = log(piA*kappa*beta);
    lMat[2][1] = log(piC*beta);
    lMat[2][2] = log(lambda - beta*(piY + piA*kappa));
    lMat[2][3] = log(piT*beta);
    lMat[3][0] = log(piA*beta);
    lMat[3][1] = log(piC*kappa*beta);
    lMat[3][2] = log(piG*beta);
    lMat[3][3] = log(lambda - beta*(piR + piC*kappa));

    // here are the values on the normal scale in case they are ever needed
    //lMat[0][0] = 1.0 - beta*(piY + piG*kappa)/lambda
    //lMat[0][1] = piC*beta/lambda
    //lMat[0][2] = piG*kappa*beta/lambda
    //lMat[0][3] = piT*beta/lambda
    //lMat[1][0] = piA*beta/lambda
    //lMat[1][1] = 1.0 - beta*(piR + piT*kappa)/lambda
    //lMat[1][2] = piG*beta/lambda
    //lMat[1][3] = piT*kappa*beta/lambda
    //lMat[2][0] = piA*kappa*beta/lambda
    //lMat[2][1] = piC*beta/lambda
    //lMat[2][2] = 1.0 - beta*(piY + piA*kappa)/lambda
    //lMat[2][3] = piT*beta/lambda
    //lMat[3][0] = piA*beta/lambda
    //lMat[3][1] = piC*kappa*beta/lambda
    //lMat[3][2] = piG*beta/lambda
    //lMat[3][3] = 1.0 - beta*(piR + piC*kappa)/lambda

    return lambda;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function 
|   inherited from the base class Model. 
*/
double HKY::calcUMat(double * * uMat) const
	{
	double piA = state_freqs[0];
	double piC = state_freqs[1];
	double piG = state_freqs[2];
	double piT = state_freqs[3];
	double piR = piA + piG;
	double piY = piC + piT;

    // set beta such that edgelen = t
    double beta = 1.0/(2.0*piR*piY + 2.0*kappa*(piA*piG + piC*piT));
        
    // set lambda to maximim substitution rate + epsilon
    double lambda_epsilon = 0.1;
    double max_rate = beta*(piY + piG*kappa);
    double next_rate = beta*(piR + piT*kappa);
    if (next_rate > max_rate)
        max_rate = next_rate;
    next_rate = beta*(piY + piA*kappa);
    if (next_rate > max_rate)
        max_rate = next_rate;
    next_rate = beta*(piR + piC*kappa);
    if (next_rate > max_rate)
        max_rate = next_rate;
    double lambda = max_rate + lambda_epsilon;

    uMat[0][0] = 1.0 - beta*(piY + piG*kappa)/lambda;
    uMat[0][1] = piC*beta/lambda;
    uMat[0][2] = piG*kappa*beta/lambda;
    uMat[0][3] = piT*beta/lambda;
    uMat[1][0] = piA*beta/lambda;
    uMat[1][1] = 1.0 - beta*(piR + piT*kappa)/lambda;
    uMat[1][2] = piG*beta/lambda;
    uMat[1][3] = piT*kappa*beta/lambda;
    uMat[2][0] = piA*kappa*beta/lambda;
    uMat[2][1] = piC*beta/lambda;
    uMat[2][2] = 1.0 - beta*(piY + piA*kappa)/lambda;
    uMat[2][3] = piT*beta/lambda;
    uMat[3][0] = piA*beta/lambda;
    uMat[3][1] = piC*kappa*beta/lambda;
    uMat[3][2] = piG*beta/lambda;
    uMat[3][3] = 1.0 - beta*(piR + piC*kappa)/lambda;

    return lambda;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the transition probability matrix given an edge length. Overrides the pure virtual function inherited from 
|	the base class Model.
*/
void HKY::calcPMat(double * * pMat, double edgeLength) const
	{
	PHYCAS_ASSERT(state_freqs.size() == 4);
	double piA = state_freqs[0];
	double piC = state_freqs[1];
	double piG = state_freqs[2];
	double piT = state_freqs[3];

	double PiA = piA + piG;
	double PiC = piC + piT;
	double PiG = piA + piG;
	double PiT = piC + piT;

	double bigPiInvA = 1.0/PiA;
	double bigPiInvC = 1.0/PiC;
	double bigPiInvG = 1.0/PiG;
	double bigPiInvT = 1.0/PiT;

	double t = edgeLength;
    // The next two lines fix the "Rota" bug; see BUGS file for details
    if (t < 1.e-8) 
        t = 1.e-8; //TreeNode::edgeLenEpsilon;
	double ta, tb, tc, td, y;
	double denom = ((piA + piG)*(piC + piT) + kappa*((piA*piG) + (piC*piT)));
	double beta = 0.5/denom;
	double x = exp(-beta*t);

	// changes to base A
	td			= -beta*(1 + PiA*(kappa - 1.0));
	y			= exp(t*td);
	ta			= piA*(bigPiInvA - 1.0);
	tb			= (PiA - piA)*bigPiInvA;
	tc			= piA*bigPiInvA;
	pMat[0][0]	= piA + (x*ta) + (y*tb);
	pMat[1][0]	= piA*(1.0 - x);
	pMat[2][0]	= piA + (x*ta) - (y*tc);
	pMat[3][0]	= pMat[1][0];

	// changes to base C
	td = -beta*(1 + PiC*(kappa - 1.0));
	y = exp(t*td);
	ta = piC*(bigPiInvC - 1.0);
	tb = (PiC - piC)*bigPiInvC;
	tc = piC*bigPiInvC;
	pMat[0][1] = piC*(1.0 - x);
	pMat[1][1] = piC + (x*ta) + (y*tb);
	pMat[2][1] = pMat[0][1];
	pMat[3][1] = piC + (x*ta) - (y*tc);

	// changes to base G
	td = -beta*(1 + PiG*(kappa - 1.0));
	y = exp(t*td);
	ta = piG*(bigPiInvG - 1.0);
	tb = (PiG - piG)*bigPiInvG;
	tc = piG*bigPiInvG;
	pMat[0][2] = piG + (x*ta) - (y*tc);
	pMat[1][2] = piG*(1.0 - x);
	pMat[2][2] = piG + (x*ta) + (y*tb);
	pMat[3][2] = pMat[1][2];

	// changes to base T
	td = -beta*(1 + PiT*(kappa - 1.0));
	y = exp(t*td);
	ta = piT*(bigPiInvT - 1.0);
	tb = (PiT - piT)*bigPiInvT;
	tc = piT*bigPiInvT;
	pMat[0][3] = piT*(1.0 - x);
	pMat[1][3] = piT + (x*ta) - (y*tc);
	pMat[2][3] = pMat[0][3];
	pMat[3][3] = piT + (x*ta) + (y*tb);
	}
