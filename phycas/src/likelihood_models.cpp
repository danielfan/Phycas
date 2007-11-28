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
			{ //POL_BOOKMARK
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
