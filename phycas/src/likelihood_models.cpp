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
#include "phycas/src/basic_tree.hpp"
#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
#	include <boost/python/numeric.hpp>
#	include "phycas/src/thirdparty/num_util/num_util.h"
#endif
using std::cout;
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The non-public constructor has a single parameter representing the number of basic states this model supports.
|	State frequencies are all set to 1/`num_states'. Assumes that the number of states provided is greater than 0.
*/
Model::Model(
  unsigned numStates)	/**< is the number of basic states (e.g. 4 for DNA) */
	:
	subset_index(-1),
	internal_edgelen_hyperparam(0.1),
	external_edgelen_hyperparam(0.1),
#if POLPY_NEWWAY
   	separate_edgelen_params(false),
#endif
   	separate_int_ext_edgelen_priors(false),
	state_freq_fixed(false),
	edge_lengths_fixed(false),	
	edgelen_hyperprior_fixed(false),
	pinvar_fixed(false),
	gamma_shape_fixed(false),
    invert_shape(false),
    time_stamp(0),
	num_states(numStates), //@POL 31-Oct-2005 what about morphological models, is numStates in this case the maximum number of states allowed?
	num_gamma_rates(1), 
	gamma_rates_unnorm(1, 1.0),
	gamma_rate_probs(1, 1.0), 
	gamma_shape(0.5),
	pinvar(0.0), 
	is_codon_model(false),
	is_pinvar_model(false)
	{
	PHYCAS_ASSERT(num_states > 0);
	setAllFreqsEqual();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The virtual destructor .
*/
Model::~Model()
	{
	//std::cerr << "\n>>>>> Model dying..." << std::endl;
	Clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears vectors of shared pointers.
*/
void Model::Clear()
	{
	//std::cerr << "\n>>>>> Model::Clear..." << std::endl;
	freq_params.clear();
	edgelen_hyper_params.clear();
	edgelen_params.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets shared pointers to MCMCUpdater objects to break cyclic dependency that prevents shared pointers from ever 
|	reaching a use count of 1 (and hence preventing deleting of the pointed-to objects). The cycle is 
|	TreeLikelihood->Model->MCMCUpdater->TreeLikelihood.
*/
void Model::releaseUpdaters()
    {
	//std::cerr << "\n>>>>> Model::releaseUpdaters()..." << std::endl;    
	Clear();
	edgeLenHyperPrior.reset();
	internalEdgeLenPrior.reset();
	externalEdgeLenPrior.reset();
	pinvar_param.reset();	
	pinvar_prior.reset();
	gamma_shape_param.reset();
	gamma_shape_prior.reset();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `is_codon_model'.
*/
bool Model::isCodonModel() const
	{
	return is_codon_model;
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of data member `separate_edgelen_params' to supplied value `param_for_each_edgelen'.
*/
void Model::setEdgeSpecificParams(
  bool param_for_each_edgelen)
    {
    separate_edgelen_params = param_for_each_edgelen;
    }
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of data member `separate_int_ext_edgelen_priors' to supplied value `separate'.
*/
void Model::separateInternalExternalEdgeLenPriors(
  bool separate)
    {
    separate_int_ext_edgelen_priors = separate;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `separate_int_ext_edgelen_priors'.
*/
bool Model::isSeparateInternalExternalEdgeLenPriors() const
    {
    return separate_int_ext_edgelen_priors;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `edgeLenHyperPrior' to the supplied probability distribution `d'. This prior distribution will 
|	be used when the model is asked to create parameters.
*/
void	Model::setEdgeLenHyperPrior(ProbDistShPtr d)		
	{
	edgeLenHyperPrior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that can be used to copy the `internalEdgeLenPrior' shared pointer.
*/
ProbDistShPtr Model::getInternalEdgeLenPrior()
	{
	return internalEdgeLenPrior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `internalEdgeLenPrior' to the supplied probability distribution `d'. This prior distribution
|	will be used when the model is asked to create parameters.
*/
void	Model::setInternalEdgeLenPrior(ProbDistShPtr d)
	{
	internalEdgeLenPrior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that can be used to copy the `externalEdgeLenPrior' shared pointer.
*/
ProbDistShPtr Model::getExternalEdgeLenPrior()
	{
	return externalEdgeLenPrior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `externalEdgeLenPrior' to the supplied probability distribution `d'. This prior distribution
|	will be used when the model is asked to create parameters.
*/
void Model::setExternalEdgeLenPrior(ProbDistShPtr d)
	{
	externalEdgeLenPrior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that can be used to copy the `edgeLenHyperPrior' shared pointer.
*/
ProbDistShPtr Model::getEdgeLenHyperPrior()
	{
	return edgeLenHyperPrior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Predicate that returns true if `edgeLenHyperPrior' actually points to something, and returns false if no edge length
|	hyperprior probability distribution was assigned.
*/
bool Model::hasEdgeLenHyperPrior()
	{
	if (edgeLenHyperPrior)
		return true;
	else 
		return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes all `numRates' transition probability matrices. Assumes edge lengths in `edgeLength' array have already
|	been computed.
*/
void Model::calcPMatrices(
  double * * *		pMat,			/**< is the array of 2-dimensional transition probability matrices (one transition matrix for each relative rate category) */
  const double *	edgeLength,		/**< is the vector of rate-adjusted edge lengths (length is `numRates') */
  unsigned			numRates		/**< is the number of relative rate categories */
  ) const
	{
	for (unsigned i = 0; i < numRates; ++i)
		calcPMat(pMat[i], edgeLength[i]);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns value of data member `num_states'.
*/
unsigned Model::getNStates() const
	{
	return num_states;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns a reference to `state_freqs', the vector of relative state frequencies. The values in
|	the array returned by this function have been normalized such that they sum to 1.0.
*/
const std::vector<double> & Model::getStateFreqs() const
	{
	return state_freqs;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets all four nucleotide frequency parameters. The member function normalizeFreqs is called 
|	automatically to recalculate `state_freqs', so it is not necessary to supply four values that sum to 1.0. Assumes
|	that `num_states' equals 4 and that each of the four values supplied are greater than or equal to zero.
*/
void Model::setNucleotideFreqs(
  double freqA,				/**< the new value of `state_freq_unnorm'[0] (i.e. frequency of base A) */
  double freqC,				/**< the new value of `state_freq_unnorm'[1] (i.e. frequency of base C) */
  double freqG,				/**< the new value of `state_freq_unnorm'[2] (i.e. frequency of base G) */
  double freqT)				/**< the new value of `state_freq_unnorm'[3] (i.e. frequency of base T/U) */
	{
	++time_stamp;
	PHYCAS_ASSERT(num_states == 4);
	PHYCAS_ASSERT(freqA >= 0.0);
	PHYCAS_ASSERT(freqC >= 0.0);
	PHYCAS_ASSERT(freqG >= 0.0);
	PHYCAS_ASSERT(freqT >= 0.0);
	state_freq_unnorm[0] = freqA;
	state_freq_unnorm[1] = freqC;
	state_freq_unnorm[2] = freqG;
	state_freq_unnorm[3] = freqT;
	normalizeFreqs();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets all of the state frequency parameters (the unnormalized values that determine the values
|	in the `state_freqs' vector when normalized) in the data member vector `state_freq_unnorm'. The member function
|	normalizeFreqs is called automatically to recalculate `state_freqs'. 
*/
void Model::setStateFreqsUnnorm(
  const std::vector<double> & values)	/**< the new values */
	{
	++time_stamp;
    PHYCAS_ASSERT(values.size() == state_freq_unnorm.size());
    std::copy(values.begin(), values.end(), state_freq_unnorm.begin());
	normalizeFreqs();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets one of the state frequency parameters (the unnormalized values that determine the values
|	in the `state_freqs' vector when normalized) in the data member vector `state_freq_unnorm'. The member function
|	normalizeFreqs is called automatically to recalculate `state_freqs'. Assumes `param_index' is less than `num_states'
|	(i.e. a valid index into `state_freq_unnorm) and `value' is non-negative.
*/
void Model::setStateFreqUnnorm(
  unsigned param_index,		/**< the 0-based index into the `state_freq_unnorm' vector of the element to modify */
  double value)				/**< the new value of `state_freq_unnorm'[`param_index'] */
	{
	++time_stamp;
	PHYCAS_ASSERT(param_index < num_states);
	PHYCAS_ASSERT(value >= 0.0);
	state_freq_unnorm[param_index] = value;
	normalizeFreqs();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets all state frequency parameters in the data member `stateFreqParams' to 1.0, and all
|	state frequencies in the data member `state_freqs' to the value 1/`num_states'. Assumes `num_states' is greater than zero.
*/
void Model::setAllFreqsEqual()
	{
	++time_stamp;
	PHYCAS_ASSERT(num_states > 0);
	state_freq_unnorm.clear();
	state_freq_unnorm.assign(num_states, 1.0);
	state_freqs.clear();
	state_freqs.assign(num_states, 1.0/num_states);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets all state frequencies in data member `state_freqs' to the normalized values in the 
|	vector `stateFreqParams'. Assumes `num_states' is greater than zero.
*/
void Model::normalizeFreqs()
	{
	PHYCAS_ASSERT(num_states > 0);
	double sum = std::accumulate(state_freq_unnorm.begin(), state_freq_unnorm.end(), 0.0);
	PHYCAS_ASSERT(sum != 0.0);
	PHYCAS_ASSERT(state_freq_unnorm.size() == state_freqs.size());
	std::transform(state_freq_unnorm.begin(), state_freq_unnorm.end(), state_freqs.begin(), boost::lambda::_1/sum);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the total number of rate categories. The total number of rate categories equals the number of discrete 
|	gamma rate categories (i.e. the value returned by the getNGammaRates member function) if `is_pinvar_model' is false.
|	If `is_pinvar_model' is true, however, the total number of rate categories equals the number of discrete gamma rate
|	categories plus 1 (to accommodate the 0-rate category of which the probability is `pinvar').
*/
unsigned Model::getNRatesTotal() const
	{
	return num_gamma_rates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function returning the value of the data member `num_gamma_rates'. This is the number of discrete gamma
|	rate categories. The total number of relative rate categories exceeds this value by 1 if `pinvar' is greater than
|	0.0. See the documentation for the getNRatesTotal member function.
*/
unsigned Model::getNGammaRates() const
	{
	return num_gamma_rates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the value of data member `num_gamma_rates' to the supplied number of rate categories 
|	`nGammaRates'. Recomputes the vector `gamma_rate_probs' if `num_gamma_rates' is changed. Each element of 
|	`gamma_rate_probs' is assigned the value 1/`nGammaRates'. Assumes `nGammaRates' > 0.
*/
void Model::setNGammaRates(
  unsigned nGammaRates)				/**< is the new number of discrete gamma rate categories */
	{
	if (nGammaRates != num_gamma_rates)
		{
		++time_stamp;
		PHYCAS_ASSERT(nGammaRates > 0);
		gamma_rates_unnorm.assign(nGammaRates, 1.0);
		gamma_rate_probs.resize(nGammaRates); //@POL this line not necessary (?) because assign also resizes
		gamma_rate_probs.assign(nGammaRates, 1.0/(double)nGammaRates);
		num_gamma_rates = nGammaRates;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns a reference to `gamma_rate_probs', the vector of probabilities of a site being in 
|	any given rate category.
*/
const std::vector<double> & Model::getGammaRateProbs() const
	{
	return gamma_rate_probs;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets all values in the vector `gamma_rate_probs' to the value 1/`num_gamma_rates'. Assumes 
|	`num_gamma_rates' is greater than zero.
*/
void Model::setAllGammaRateProbsEqual()
	{
	++time_stamp;
	PHYCAS_ASSERT(num_gamma_rates > 0);
	gamma_rate_probs.clear();
	gamma_rate_probs.assign(num_gamma_rates, 1.0/num_gamma_rates);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the current value of the data member `gamma_shape'.
*/
double Model::getShape()
	{
	return gamma_shape;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the data member `gamma_shape' to the supplied value `alpha'. Assumes `alpha' is greater
|	than zero.
*/
void Model::setShape(
  double alpha)		/**< is the new value for the `gamma_shape' data member */
	{
	++time_stamp;
	PHYCAS_ASSERT(alpha > 0.0);
	gamma_shape = alpha;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the data member `invert_shape' to the supplied value `invert'. If true is specified for
|	`invert', then `gamma_shape_param' will manage the inverse of the gamma shape rather than the shape itself.
*/
void Model::setPriorOnShapeInverse(
  bool invert)		/**< is the new value for the `gamma_shape' data member */
	{
	invert_shape = invert;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `gamma_shape_fixed' to true. The fixParameter member function of the DiscreteGammaParam object
|	is either called immediately (if `gamma_shape_param' is a valid pointer) or is called in createParameters (when
|	`gamma_shape_param' is first assigned).
*/
void Model::fixShape()
	{
	gamma_shape_fixed = true;
	if (gamma_shape_param)
		gamma_shape_param->fixParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `gamma_shape_fixed' to false. The freeParameter member function of the DiscreteGammaParam 
|	object is called immediately if `gamma_shape_param' is a valid pointer.
*/
void Model::freeShape()
	{
	gamma_shape_fixed = false;
	if (gamma_shape_param)
		gamma_shape_param->freeParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the current value of the data member `pinvar'.
*/
double Model::getPinvar()
	{
	return pinvar;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the data member `pinvar' to the supplied value `pinv'. Assumes `pinv' is greater
|	than or equal to zero but strictly less than 1.0.
*/
void Model::setPinvar(
  double pinv)		/**< is the new value for the `pinvar' data member */
	{
	++time_stamp;
	PHYCAS_ASSERT(pinv >= 0.0);
	PHYCAS_ASSERT(pinv < 1.0);
	pinvar = pinv;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `gamma_shape_fixed'.
*/
bool Model::shapeFixed() const
	{
	return gamma_shape_fixed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `pinvar_fixed'.
*/
bool Model::pinvarFixed() const
	{
	return pinvar_fixed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `state_freq_fixed'.
*/
bool Model::stateFreqsFixed() const
	{
	return state_freq_fixed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `edgelen_hyperprior_fixed'.
*/
bool Model::edgeLenHyperParamFixed() const
	{
	return edgelen_hyperprior_fixed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `edge_lengths_fixed'.
*/
bool Model::edgeLengthsFixed() const
	{
	return edge_lengths_fixed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `edge_lengths_fixed' to true.
*/
void Model::fixEdgeLengths()
	{
	edge_lengths_fixed = true;
	if (!edgelen_params.empty())
		{
		//@POL good place to use std::for_each with a boost::lambda functor?
		for (MCMCUpdaterVect::iterator it = edgelen_params.begin(); it != edgelen_params.end(); ++it)
			(*it)->fixParameter();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `edge_lengths_fixed' to false.
*/
void Model::freeEdgeLengths()
	{
	edge_lengths_fixed = false;
	if (!edgelen_params.empty())	//TODO create Model::edgelen_params and fill it in createParameters
		{
		//@POL good place to use std::for_each with a boost::lambda functor?
		for (MCMCUpdaterVect::iterator it = edgelen_params.begin(); it != edgelen_params.end(); ++it)
			(*it)->freeParameter();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `edgelen_hyperprior_fixed' to true and calls fixParameter() for all edge length 
|   hyperparameters.
*/
void Model::fixEdgeLenHyperprior()
	{
	edgelen_hyperprior_fixed = true;
	if (!edgelen_hyper_params.empty())
		{
		for (MCMCUpdaterVect::iterator it = edgelen_hyper_params.begin(); it != edgelen_hyper_params.end(); ++it)
			(*it)->fixParameter();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `edgelen_hyperprior_fixed' to false and calls freeParameter() for all edge length 
|   hyperparameters.
*/
void Model::freeEdgeLenHyperprior()
	{
	edgelen_hyperprior_fixed = false;
	if (!edgelen_hyper_params.empty())
		{
		for (MCMCUpdaterVect::iterator it = edgelen_hyper_params.begin(); it != edgelen_hyper_params.end(); ++it)
			(*it)->freeParameter();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `state_freq_fixed' to true. The fixParameter member function of all StateFreqParam objects is 
|	either called immediately (if the `freq_params' vector is not empty) or is called in createParameters (when 
|	`freq_params' is built).
*/
void Model::fixStateFreqs()
	{
	state_freq_fixed = true;
	if (!freq_params.empty())
		{
		//@POL good place to use std::for_each with a boost::lambda functor?
		for (MCMCUpdaterVect::iterator it = freq_params.begin(); it != freq_params.end(); ++it)
			(*it)->fixParameter();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `state_freq_fixed' to false. The freeParameter member function of all StateFreqParam objects is 
|	called immediately if the `freq_params' vector is not empty.
*/
void Model::freeStateFreqs()
	{
	state_freq_fixed = false;
	if (!freq_params.empty())
		{
		//@POL good place to use std::for_each with a boost::lambda functor?
		for (MCMCUpdaterVect::iterator it = freq_params.begin(); it != freq_params.end(); ++it)
			(*it)->freeParameter();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `is_pinvar_model'.
*/
bool Model::isPinvarModel()
	{
	return is_pinvar_model;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `is_pinvar_model' to true. A subsequent call to the createParameters member function will
|	result in a PinvarParam being added to the list of updaters for this model.
*/
void Model::setPinvarModel()
	{
	++time_stamp;
	is_pinvar_model = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `is_pinvar_model' to false. No PinvarParam will be added in a subsequent call to the 
|	createParameters member function.
*/
void Model::setNotPinvarModel()
	{
	++time_stamp;
	is_pinvar_model = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `pinvar_fixed' to true. The fixParameter member function of the PinvarParam object is either 
|	called immediately (if `pinvar_param' is a valid pointer) or is called in createParameters (when `pinvar_param' is 
|	first assigned).
*/
void Model::fixPinvar()
	{
	pinvar_fixed = true;
	if (pinvar_param)
		pinvar_param->fixParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `pinvar_fixed' to false. The freeParameter member function of the PinvarParam object is called 
|	immediately if `pinvar_param' is a valid pointer.
*/
void Model::freePinvar()
	{
	pinvar_fixed = false;
	if (pinvar_param)
		pinvar_param->freeParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `gamma_shape_prior'.
*/
ProbDistShPtr Model::getDiscreteGammaShapePrior()
 	{
	return gamma_shape_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of data member `gamma_shape_prior'.
*/
void Model::setDiscreteGammaShapePrior(ProbDistShPtr d)
 	{
	gamma_shape_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `pinvar_prior'.
*/
ProbDistShPtr Model::getPinvarPrior()
 	{
	return pinvar_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `pinvar_prior'.
*/
void Model::setPinvarPrior(ProbDistShPtr d)
 	{
	pinvar_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears `state_list' and `state_list_pos', then builds up a 4-state version of these vectors suitable for simulated
|	nucleotide data (i.e. no ambiguities).
*/
void Model::buildStateList(
  state_list_t & state_list,				/**< is the vector of state codes (see TipData documentation for explanation) */
  state_list_pos_t & state_list_pos)	/**< is the vector of positions of states within `state_list' */
  const
	{
	state_list.clear();
	state_list_pos.clear();

	// state A 
	state_list.push_back(1);		// next 1 element is state code for base A
	state_list.push_back(0);		// state code for base A is 0
	state_list_pos.push_back(0);	// information for first state starts at position 0 in state_list

	// state C
	state_list.push_back(1);		// next 1 element is state code for base C
	state_list.push_back(1);		// state code for base C is 1
	state_list_pos.push_back(2);	// information for second state starts at position 2 in state_list

	// state G 
	state_list.push_back(1);		// next 1 element is state code for base G
	state_list.push_back(2);		// state code for base G is 2
	state_list_pos.push_back(4);	// information for third state starts at position 4 in state_list

	// state T 
	state_list.push_back(1);		// next 1 element is state code for base T
	state_list.push_back(3);		// state code for base T is 3
	state_list_pos.push_back(6);	// information for fourth state starts at position 6 in state_list
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string in the `state_repr' vector corresponding to `state'. If supplied `state' is negative, returns
|	the gap state ("-"). If `state' equals or exceeds `num_states', returns "?".
*/
std::string Model::lookupStateRepr(
  int state)	const /**< is the state to look up */
	{
	if (state < 0)
		return "-";
	else if ((unsigned)state < num_states)
		return state_repr[state];
	else
		return "?";
	}

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
	std::vector<double> boundaries;

	if (is_pinvar_model)
		{
		PHYCAS_ASSERT(pinvar < 1.0);
		double one_minus_pinvar = 1.0 - pinvar;

		// First calculate the gamma rates alone
		std::vector<double> gamma_rates;
		recalcGammaRatesAndBoundaries(gamma_rates, boundaries);

		// Copy the rate probabilities directly
		probs.resize(num_gamma_rates, 0.0);
		std::copy(gamma_rate_probs.begin(), gamma_rate_probs.end(), probs.begin());

		// Adjust gamma_rates using pinvar and save to rates vector
		rates.resize(num_gamma_rates, 0.0);
		std::vector<double>::iterator       rates_iter       = rates.begin();
		std::vector<double>::const_iterator gamma_rates_iter = gamma_rates.begin();
		for (unsigned i = 0; i < num_gamma_rates; ++i)
			{
			*rates_iter++ = (*gamma_rates_iter++)/one_minus_pinvar;
			}
		}
	else
		{
		// The rates and probs computed by recalcGammaRatesAndBoundaries are what we need if is_pinvar_model is false
		recalcGammaRatesAndBoundaries(rates, boundaries);
		probs.resize(num_gamma_rates, 0.0);
		std::copy(gamma_rate_probs.begin(), gamma_rate_probs.end(), probs.begin());
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
  MCMCUpdaterVect & parameters_vect_ref,		    /**< is the vector of model-specific parameters to fill */
  int subset_pos) 									/**< if 0 (first subset) or -1 (only subset), edge length parameters and hyperparams will be added; otherwise, the `edgelens' and `edgelen_hyperparams' vectors returned will be empty */
    {
	PHYCAS_ASSERT(t);
	PHYCAS_ASSERT(edgelens_vect_ref.empty());
	PHYCAS_ASSERT(edgelen_hyperparams_vect_ref.empty());
	PHYCAS_ASSERT(parameters_vect_ref.empty());
	
	subset_index = subset_pos;
	if (subset_index <= 0)
		{
		// Add the edge length parameter(s)
#if POLPY_NEWWAY
		if (separate_edgelen_params)
			{
			// Add an EdgeLengthParam for every node that has an associated edge length
			for (preorder_iterator nd = t->begin(); nd != t->end(); ++nd)
				{
				if (nd->IsAnyRoot())
					continue;
				EdgeLenParam * p = new EdgeLenParam();
				p->setTreeNode(*nd);
				MCMCUpdaterShPtr edge_len_param = MCMCUpdaterShPtr(p);
				edge_len_param->setName(boost::str(boost::format("edgelen_%d") % nd->GetNodeNumber())); 
				edge_len_param->setStartingValue(nd->GetEdgeLen());
				edge_len_param->setTree(t);
				if (separate_int_ext_edgelen_priors && nd->IsInternal())
					edge_len_param->setPrior(internalEdgeLenPrior);
				else 
					edge_len_param->setPrior(externalEdgeLenPrior);
				if (edge_lengths_fixed)
					edge_len_param->fixParameter();	
				edgelens_vect_ref.push_back(edge_len_param);
				}
			}
		else 	// add one or two EdgeLenMasterParams who have the job of computing the prior
#endif

		if (separate_int_ext_edgelen_priors)
			{
			// Add two edge length parameters to manage the priors for all edge lengths in the tree;
			// The first of these two parameters will be in charge of the prior on external edge lengths,
			// whereas the second will be in charge of internal edge lengths. These edge length parameters 
			// do not actually update any edge lengths: a Metropolis proposal such as the LargetSimonMove 
			// is responsible for updating edge lengths.

			// First the parameter governing external edge length priors
			MCMCUpdaterShPtr p_external = MCMCUpdaterShPtr(new EdgeLenMasterParam(EdgeLenMasterParam::external));
			p_external->setName(std::string("external_edgelen"));
			p_external->setTree(t);
			p_external->setPrior(externalEdgeLenPrior);
			if (edge_lengths_fixed)
				p_external->fixParameter();
			edgelens_vect_ref.push_back(p_external);
			
			if (edgeLenHyperPrior)
				{
				EdgeLenMasterParamShPtr pit = boost::shared_dynamic_cast<EdgeLenMasterParam>(p_external);
				MCMCUpdaterShPtr p = MCMCUpdaterShPtr(new HyperPriorParam(pit, true));
				p->setName(std::string("external_hyper"));
				p->setTree(t);
				p->setPrior(edgeLenHyperPrior);
				if (edgelen_hyperprior_fixed)
					p->fixParameter();
				edgelen_hyperparams_vect_ref.push_back(p);
				}

			// Now the parameter governing internal edge length priors
			MCMCUpdaterShPtr p_internal = MCMCUpdaterShPtr(new EdgeLenMasterParam(EdgeLenMasterParam::internal));
			p_internal->setName(std::string("internal_edgelen"));
			p_internal->setTree(t);
			p_internal->setPrior(internalEdgeLenPrior);
			if (edge_lengths_fixed)
				p_internal->fixParameter();
			edgelens_vect_ref.push_back(p_internal);

			if (edgeLenHyperPrior)
				{
				EdgeLenMasterParamShPtr pit = boost::shared_dynamic_cast<EdgeLenMasterParam>(p_internal);
				MCMCUpdaterShPtr p = MCMCUpdaterShPtr(new HyperPriorParam(pit, false));
				p->setName(std::string("internal_hyper"));
				p->setTree(t);
				p->setPrior(edgeLenHyperPrior);
				if (edgelen_hyperprior_fixed)
					p->fixParameter();
				edgelen_hyperparams_vect_ref.push_back(p);
				}
			}
		else
			{
			// Add one edge length parameter to manage the priors for all edge lengths in the tree. 
			// This edge length parameter does not actually update any edge lengths: a Metropolis proposal 
			// such as the LargetSimonMove is responsible for updating edge lengths.
			MCMCUpdaterShPtr p = MCMCUpdaterShPtr(new EdgeLenMasterParam(EdgeLenMasterParam::both));
			p->setName(std::string("master_edgelen"));
			p->setTree(t);
			p->setPrior(externalEdgeLenPrior);
			if (edge_lengths_fixed)
				p->fixParameter();
			edgelens_vect_ref.push_back(p);
			//std::cerr << "\n>>>>> MASTER EDGE LENGH PARAMETER use count = 2? " << p.use_count() << std::endl;

			if (edgeLenHyperPrior)
				{
				EdgeLenMasterParamShPtr pit = boost::shared_dynamic_cast<EdgeLenMasterParam>(p);
				MCMCUpdaterShPtr p_hyper = MCMCUpdaterShPtr(new HyperPriorParam(pit, true));
				p_hyper->setName(std::string("edgelen_hyper"));
				p_hyper->setTree(t);
				p_hyper->setPrior(edgeLenHyperPrior);
				if (edgelen_hyperprior_fixed)
					p_hyper->fixParameter();
				edgelen_hyperparams_vect_ref.push_back(p_hyper);
				}
			}
		
		// Save a vector of shared pointers to edge length parameters so that we can modify their fixed/free status if we need to
		edgelen_params.resize(edgelens_vect_ref.size());
		std::copy(edgelens_vect_ref.begin(), edgelens_vect_ref.end(), edgelen_params.begin());

		// Save a vector of shared pointers to edge length hyperparameters so that we can modify their fixed/free status if we need to
		edgelen_hyper_params.resize(edgelen_hyperparams_vect_ref.size());
		std::copy(edgelen_hyperparams_vect_ref.begin(), edgelen_hyperparams_vect_ref.end(), edgelen_hyper_params.begin());
	}

	if (num_gamma_rates > 1)
		{
		PHYCAS_ASSERT(num_gamma_rates > 1);
		PHYCAS_ASSERT(!gamma_shape_param);
		gamma_shape_param = MCMCUpdaterShPtr(new DiscreteGammaShapeParam(invert_shape));
        if (invert_shape)
			{
			if (subset_pos < 0)
				gamma_shape_param->setName(std::string("gamma_var")); 
			else
				gamma_shape_param->setName(boost::str(boost::format("gamma_var_%d") % (subset_pos + 1)));
			}
        else
			{
			if (subset_pos < 0)
				gamma_shape_param->setName(std::string("gamma_shape")); 
			else
				gamma_shape_param->setName(boost::str(boost::format("gamma_shape_%d") % (subset_pos + 1)));
			}
        gamma_shape_param->setStartingValue(gamma_shape);
		gamma_shape_param->setTree(t);
		gamma_shape_param->setPrior(gamma_shape_prior);
		if (gamma_shape_fixed)
			gamma_shape_param->fixParameter();
		parameters_vect_ref.push_back(gamma_shape_param);
		//std::cerr << "\n>>>>> gamma_shape_param.use_count() = " << gamma_shape_param.use_count() << std::endl;
		}

	if (is_pinvar_model)
		{
		PHYCAS_ASSERT(!pinvar_param);
		pinvar_param = MCMCUpdaterShPtr(new PinvarParam());
		if (subset_pos < 0)
			pinvar_param->setName(std::string("pinvar")); 
		else
			pinvar_param->setName(boost::str(boost::format("pinvar_%d") % (subset_pos + 1)));
        pinvar_param->setStartingValue(pinvar);
		pinvar_param->setTree(t);
		pinvar_param->setPrior(pinvar_prior);
		if (pinvar_fixed)
			pinvar_param->fixParameter();
		parameters_vect_ref.push_back(pinvar_param);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by 
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for 
|	compatability with MrBayes). This base class version provides columns for rate heterogeneity and edge length
|	hyperparameters.
*/
std::string Model::paramHeader() const
	{
	std::string s;
	for (MCMCUpdaterVect::const_iterator it = edgelen_hyper_params.begin(); it != edgelen_hyper_params.end(); ++it)
		{
		s += boost::str(boost::format("\t%s") % (*it)->getName());
		}
		
	if (is_pinvar_model)
		s += boost::str(boost::format("\t%s") % pinvar_param->getName());

	if (num_gamma_rates > 1)
		{
		s += boost::str(boost::format("\t%s") % gamma_shape_param->getName());
		}
	return s;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by 
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for 
|	compatability with MrBayes). This base class version provides columns for rate heterogeneity and edge length
|	hyperparameters.
*/
std::string Model::paramReport(
  unsigned ndecimals) const	/**< is the number of decimal places precision to use when reporting parameter values */
	{
    std::string fmt = boost::str(boost::format("%%.%df\t") % ndecimals);
	std::string s;
	if (edgeLenHyperPrior)
		{
		if (separate_int_ext_edgelen_priors)
			{
			s += boost::str(boost::format(fmt) % external_edgelen_hyperparam);
			s += boost::str(boost::format(fmt) % internal_edgelen_hyperparam);
			}
		else 
			{
			s += boost::str(boost::format(fmt) % external_edgelen_hyperparam);
			}
		}
		
	if (is_pinvar_model)
        {
		s += boost::str(boost::format(fmt) % pinvar);
        }
		
	if (num_gamma_rates > 1)
		{
		s += boost::str(boost::format(fmt) % gamma_shape);
		}
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Stores a flattened version of the supplied 2-dimensional array `twoDarr', storing the result in the supplied VecDbl
|	reference variable `p'. The supplied `twoDarr' should be laid out so that rows occupy contiguous memory.
*/
void Model::flattenTwoDMatrix(VecDbl & p, double * * twoDarr, unsigned dim) const
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

