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
    	separate_int_ext_edgelen_priors(false),
	num_states(numStates), //@POL 31-Oct-2005 what about morphological models, is numStates in this case the maximum number of states allowed?
	num_gamma_rates(1), 
	gamma_rates_unnorm(1, 1.0),
	gamma_rate_probs(1, 1.0), 
	state_freq_fixed(false),
	edge_lengths_fixed(false),	
	edgelen_hyperprior_fixed(false),
	is_codon_model(false),
	is_pinvar_model(false),
	is_flex_model(false),
	flex_upper_rate_bound(1.0),
	num_flex_spacers(1),
	flex_probs_fixed(false),
	flex_rates_fixed(false),
	pinvar_fixed(false),
	pinvar(0.0), 
	gamma_shape_fixed(false),
	gamma_shape(0.5),
    invert_shape(false)
	{
	PHYCAS_ASSERT(num_states > 0);
	setAllFreqsEqual();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The virtual destructor .
*/
Model::~Model()
	{
	std::cerr << "\n>>>>> Model dying..." << std::endl;
	Clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears vectors of shared pointers.
*/
void Model::Clear()
	{
	std::cerr << "\n>>>>> Model::Clear..." << std::endl;
	freq_params.clear();
	edgelen_hyper_params.clear();
	edgelen_params.clear();
	flex_rate_params.clear();
	flex_prob_params.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `is_codon_model'.
*/
bool Model::isCodonModel() const
	{
	return is_codon_model;
	}

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
void	Model::setExternalEdgeLenPrior(ProbDistShPtr d)
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
|	Sets the data member `is_flex_model' to true. A subsequent call to the createParameters member function will
|	result in `num_gamma_rates' FlexRateParam and FlexProbParam objects being added to the list of updaters for this 
|	model.
*/
void Model::setFlexModel()
	{
	is_flex_model = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `flex_upper_rate_bound' to `new_upper_bound'. Rescales all existing unnormalized relative rate
|	parameter values in the vector `gamma_rates_unnorm' to accommodate the new upper bound. Assumes that 
|	`new_upper_bound' is greater than zero. Returns immediately if `new_upper_bound' is identical to current value of
|	`flex_upper_rate_bound'.
*/
void Model::setFlexRateUpperBound(double new_upper_bound)
	{
	double old_upper_bound = flex_upper_rate_bound;
	if (new_upper_bound == old_upper_bound)
		return;
	PHYCAS_ASSERT(new_upper_bound > 0.0);
	flex_upper_rate_bound = new_upper_bound;
	double mult_factor = new_upper_bound/old_upper_bound;
	std::for_each(gamma_rates_unnorm.begin(), gamma_rates_unnorm.end(), boost::lambda::_1 *= mult_factor);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `num_flex_spacers' to `s'. See the data member `num_flex_spacers' and the function 
|	FlexRateParam::recalcPrior for more information on FLEX model spacers.
*/
void Model::setNumFlexSpacers(unsigned s)
	{
	num_flex_spacers = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `is_flex_model' to false. No FlexRateParam or FlexProbParam objects will be added in a 
|	subsequent call to the createParameters member function.
*/
void Model::setNotFlexModel()
	{
	is_flex_model = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `flex_probs_fixed' to true. The fixParameter member function of all FlexProbParam objects is 
|	either called immediately (if the `flex_prob_params' vector is not empty) or is called in createParameters (when 
|	`flex_prob_params' is built).
*/
void Model::fixFlexProbs()
	{
	flex_probs_fixed = true;
	if (!flex_prob_params.empty())
		{
		//@POL good place to use std::for_each with a boost::lambda functor?
		for (MCMCUpdaterVect::iterator it = flex_prob_params.begin(); it != flex_prob_params.end(); ++it)
			(*it)->fixParameter();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `flex_probs_fixed' to false. The freeParameter member function of all FlexProbParam objects is 
|	called immediately if the `flex_prob_params' vector is not empty.
*/
void Model::freeFlexProbs()
	{
	flex_probs_fixed = false;
	if (!flex_prob_params.empty())
		{
		//@POL good place to use std::for_each with a boost::lambda functor?
		for (MCMCUpdaterVect::iterator it = flex_prob_params.begin(); it != flex_prob_params.end(); ++it)
			(*it)->freeParameter();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `flex_rates_fixed' to true. The fixParameter member function of all FlexRateParam objects is 
|	either called immediately (if the `flex_rate_params' vector is not empty) or is called in createParameters (when 
|	`flex_rate_params' is built).
*/
void Model::fixFlexRates()
	{
	flex_rates_fixed = true;
	if (!flex_rate_params.empty())
		{
		//@POL good place to use std::for_each with a boost::lambda functor?
		for (MCMCUpdaterVect::iterator it = flex_rate_params.begin(); it != flex_rate_params.end(); ++it)
			(*it)->fixParameter();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `flex_rates_fixed' to false. The freeParameter member function of all FlexRateParam objects is 
|	called immediately if the `flex_rate_params' vector is not empty.
*/
void Model::freeFlexRates()
	{
	flex_rates_fixed = false;
	if (!flex_rate_params.empty())
		{
		//@POL good place to use std::for_each with a boost::lambda functor?
		for (MCMCUpdaterVect::iterator it = flex_rate_params.begin(); it != flex_rate_params.end(); ++it)
			(*it)->freeParameter();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of data member `flex_prob_param_prior'.
*/
void Model::setFLEXProbParamPrior(ProbDistShPtr d)
 	{
	flex_prob_param_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `flex_prob_param_prior'.
*/
ProbDistShPtr Model::getFLEXProbParamPrior()
 	{
	return flex_prob_param_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets one of the FLEX rate parameters (the unnormalized values that determine the values in 
|	the `TreeLikelihood::rate_means' vector when normalized) in the data member vector `gamma_rates_unnorm'. The 
|	function Model::recalcRatesAndProbs should be called to recalculate `TreeLikelihood::rate_means' before the 
|	likelihood is computed (ideally immediately after setFlexRateUnnorm is called). Assumes `param_index' is less than 
|	`num_gamma_rates' (i.e. a valid index into `gamma_rates_unnorm') and `value' is non-negative.
*/
void Model::setFlexRateUnnorm(
  unsigned param_index,		/**< the 0-based index into the `gamma_rates_unnorm' vector of the element to modify */
  double value)				/**< the new value of `gamma_rates_unnorm'[`param_index'] */
	{
	PHYCAS_ASSERT(gamma_rates_unnorm.size() == num_gamma_rates);
	PHYCAS_ASSERT(param_index < num_gamma_rates);
	PHYCAS_ASSERT(value >= 0.0);
	gamma_rates_unnorm[param_index] = value;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets one of the FLEX rate probability parameters (the unnormalized values that determine the 
|	values in the `TreeLikelihood::rate_probs' vector when normalized) in the data member vector `gamma_rate_probs'. 
|	The function Model::recalcRatesAndProbs should be called to recalculate `TreeLikelihood::rate_probs' before the 
|	likelihood is computed (ideally immediately after setFlexProbUnnorm is called). Assumes `param_index' is less than 
|	`num_gamma_rates' (i.e. a valid index into `gamma_rate_probs') and `value' is non-negative.
*/
void Model::setFlexProbUnnorm(
  unsigned param_index,		/**< the 0-based index into the `gamma_rate_probs' vector of the element to modify */
  double value)				/**< the new value of `gamma_rate_probs'[`param_index'] */
	{
	PHYCAS_ASSERT(gamma_rate_probs.size() == num_gamma_rates);
	PHYCAS_ASSERT(param_index < num_gamma_rates);
	PHYCAS_ASSERT(value >= 0.0);
	gamma_rate_probs[param_index] = value;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Normalizes the rates stored in `gamma_rates_unnorm' so that their mean (using probabilities in `gamma_rate_probs')
|	equals 1.0. Normalizes the probabilities stored in `gamma_rate_probs' so that their sum is 1.0. Copies normalized
|	rates to supplied vector `rates' and normalized probabilities to the supplied vector `probs'. Resizes `rates' and
|	`probs' if necessary.
*/
void Model::normalizeRatesAndProbs(
  std::vector<double> & rates,		/**< is the vector to receive the normalized rates */
  std::vector<double> & probs)		/**< is the vector to receive the normalized probabilities */
  const
	{
	// c*r0*p0 + c*r1*p1 + c*r2*p2 + c*r3*p3 = 1.0
	// c*(r0*p0 + r1*p1 + r2*p2 + r3*p3) = 1.0
	// c = 1/(r0*p0 + r1*p1 + r2*p2 + r3*p3)
	//
	PHYCAS_ASSERT(gamma_rates_unnorm.size() == num_gamma_rates);
	PHYCAS_ASSERT(gamma_rate_probs.size() == num_gamma_rates);
	rates.resize(num_gamma_rates, 0.0);
	probs.resize(num_gamma_rates, 0.0);
	std::vector<double>::const_iterator rates_it = gamma_rates_unnorm.begin();
	std::vector<double>::const_iterator probs_it = gamma_rate_probs.begin();
	unsigned i;

	// normalize the probabilities first (we will need them to normalize the rates)
	double prob_sum = 0.0;
	for (i = 0; i < num_gamma_rates; ++i)
		{
		PHYCAS_ASSERT(gamma_rate_probs[i] > 0.0);
		prob_sum += gamma_rate_probs[i];
		}
	double prob_normalizer = 1.0/prob_sum;

	// normalize the rates and copy the probs
	double rate_prob_sum = 0.0;
	for (i = 0; i < num_gamma_rates; ++i)
		{
		probs[i] = gamma_rate_probs[i]*prob_normalizer;
		PHYCAS_ASSERT(gamma_rates_unnorm[i] > 0.0);
		rate_prob_sum += gamma_rates_unnorm[i]*probs[i];
		}
	double rate_normalizer = 1.0/rate_prob_sum;

	// copy the rates
	for (i = 0; i < num_gamma_rates; ++i)
		{
		rates[i] = gamma_rates_unnorm[i]*rate_normalizer;
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
	is_pinvar_model = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `is_pinvar_model' to false. No PinvarParam will be added in a subsequent call to the 
|	createParameters member function.
*/
void Model::setNotPinvarModel()
	{
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
  VecStateList & state_list,		/**< is the vector of state codes (see TipData documentation for explanation) */
  VecStateListPos & state_list_pos)	/**< is the vector of positions of states within `state_list' */
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
		std::cerr << "\n>>>>> MASTER EDGE LENGH PARAMETER use count = 2? " << p.use_count() << std::endl;
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
		std::cerr << "\n>>>>> gamma_shape_param.use_count() = " << gamma_shape_param.use_count() << std::endl;
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

