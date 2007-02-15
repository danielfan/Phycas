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

#if ! defined(LIKELIHOOD_MODELS_INL)
#define LIKELIHOOD_MODELS_INL

#include "boost/format.hpp"	// for boost::format

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	The non-public constructor has a single parameter representing the number of basic states this model supports.
|	State frequencies are all set to 1/`num_states'. Assumes that the number of states provided is greater than 0.
*/
inline Model::Model(
  unsigned numStates)	/**< is the number of basic states (e.g. 4 for DNA) */
	:
	num_states(numStates), //@POL 31-Oct-2005 what about morphological models, is numStates in this case the maximum number of states allowed?
	num_gamma_rates(1), 
	gamma_rates_unnorm(1, 1.0),
	gamma_rate_probs(1, 1.0), 
	state_freq_fixed(false),
	edge_lengths_fixed(false),	
	edgelen_hyperprior_fixed(false),
	is_pinvar_model(false),
	is_flex_model(false),
	flex_upper_rate_bound(1.0),
	num_flex_spacers(1),
	flex_probs_fixed(false),
	flex_rates_fixed(false),
	pinvar_fixed(false),
	pinvar(0.0), 
	gamma_shape_fixed(false),
	gamma_shape(0.5)
	{
	PHYCAS_ASSERT(num_states > 0);
	setAllFreqsEqual();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The virtual destructor does nothing.
*/
inline Model::~Model()
	{
	//std::cerr << "Model dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `edgeLenHyperPrior' to the supplied probability distribution `d'. This prior distribution will 
|	be used when the model is asked to create parameters.
*/
inline void	Model::setEdgeLenHyperPrior(ProbDistShPtr d)		
	{
	edgeLenHyperPrior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `edgeLenPrior' to the supplied probability distribution `d'. This prior distribution will 
|	be used when the model is asked to create parameters.
*/
inline void	Model::setEdgeLenPrior(ProbDistShPtr d)
	{
	edgeLenPrior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that can be used to copy the `edgeLenHyperPrior' shared pointer.
*/
inline ProbDistShPtr Model::getEdgeLenHyperPrior()
	{
	return edgeLenHyperPrior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that can be used to copy the `edgeLenPrior' shared pointer.
*/
inline ProbDistShPtr Model::getEdgeLenPrior()
	{
	return edgeLenPrior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes all `numRates' transition probability matrices. Assumes edge lengths in `edgeLength' array have already
|	been computed.
*/
inline void Model::calcPMatrices(
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
inline unsigned Model::getNStates() const
	{
	return num_states;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns a reference to `state_freqs', the vector of relative state frequencies. The values in
|	the array returned by this function have been normalized such that they sum to 1.0.
*/
inline const std::vector<double> & Model::getStateFreqs() const
	{
	return state_freqs;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets all four nucleotide frequency parameters. The member function normalizeFreqs is called 
|	automatically to recalculate `state_freqs', so it is not necessary to supply four values that sum to 1.0. Assumes
|	that `num_states' equals 4 and that each of the four values supplied are greater than or equal to zero.
*/
inline void Model::setNucleotideFreqs(
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
|	Modifier function that sets one of the state frequency parameters (the unnormalized values that determine the values
|	in the `state_freqs' vector when normalized) in the data member vector `state_freq_unnorm'. The member function
|	normalizeFreqs is called automatically to recalculate `state_freqs'. Assumes `param_index' is less than `num_states'
|	(i.e. a valid index into `state_freq_unnorm) and `value' is non-negative.
*/
inline void Model::setStateFreqUnnorm(
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
inline void Model::setAllFreqsEqual()
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
inline void Model::normalizeFreqs()
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
inline unsigned Model::getNRatesTotal() const
	{
	return num_gamma_rates + (is_pinvar_model ? 1 : 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function returning the value of the data member `num_gamma_rates'. This is the number of discrete gamma
|	rate categories. The total number of relative rate categories exceeds this value by 1 if `pinvar' is greater than
|	0.0. See the documentation for the getNRatesTotal member function.
*/
inline unsigned Model::getNGammaRates() const
	{
	return num_gamma_rates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the value of data member `num_gamma_rates' to the supplied number of rate categories 
|	`nGammaRates'. Recomputes the vector `gamma_rate_probs' if `num_gamma_rates' is changed. Each element of 
|	`gamma_rate_probs' is assigned the value 1/`nGammaRates'. Assumes `nGammaRates' > 0.
*/
inline void Model::setNGammaRates(
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
inline const std::vector<double> & Model::getGammaRateProbs() const
	{
	return gamma_rate_probs;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets all values in the vector `gamma_rate_probs' to the value 1/`num_gamma_rates'. Assumes 
|	`num_gamma_rates' is greater than zero.
*/
inline void Model::setAllGammaRateProbsEqual()
	{
	PHYCAS_ASSERT(num_gamma_rates > 0);
	gamma_rate_probs.clear();
	gamma_rate_probs.assign(num_gamma_rates, 1.0/num_gamma_rates);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the current value of the data member `gamma_shape'.
*/
inline double Model::getShape()
	{
	return gamma_shape;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the data member `gamma_shape' to the supplied value `alpha'. Assumes `alpha' is greater
|	than zero.
*/
inline void Model::setShape(
  double alpha)		/**< is the new value for the `gamma_shape' data member */
	{
	PHYCAS_ASSERT(alpha > 0.0);
	gamma_shape = alpha;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the data member `invert_shape' to the supplied value `invert'. If true is specified for
|	`invert', then `gamma_shape_param' will manage the inverse of the gamma shape rather than the shape itself.
*/
inline void Model::setPriorOnShapeInverse(
  bool invert)		/**< is the new value for the `gamma_shape' data member */
	{
	invert_shape = invert;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `gamma_shape_fixed' to true. The fixParameter member function of the DiscreteGammaParam object
|	is either called immediately (if `gamma_shape_param' is a valid pointer) or is called in createParameters (when
|	`gamma_shape_param' is first assigned).
*/
inline void Model::fixShape()
	{
	gamma_shape_fixed = true;
	if (gamma_shape_param)
		gamma_shape_param->fixParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `gamma_shape_fixed' to false. The freeParameter member function of the DiscreteGammaParam 
|	object is called immediately if `gamma_shape_param' is a valid pointer.
*/
inline void Model::freeShape()
	{
	gamma_shape_fixed = false;
	if (gamma_shape_param)
		gamma_shape_param->freeParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the current value of the data member `pinvar'.
*/
inline double Model::getPinvar()
	{
	return pinvar;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the data member `pinvar' to the supplied value `pinv'. Assumes `pinv' is greater
|	than or equal to zero but strictly less than 1.0.
*/
inline void Model::setPinvar(
  double pinv)		/**< is the new value for the `pinvar' data member */
	{
	PHYCAS_ASSERT(pinv >= 0.0);
	PHYCAS_ASSERT(pinv < 1.0);
	pinvar = pinv;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `gamma_shape_fixed'.
*/
inline bool Model::shapeFixed() const
	{
	return gamma_shape_fixed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `pinvar_fixed'.
*/
inline bool Model::pinvarFixed() const
	{
	return pinvar_fixed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `state_freq_fixed'.
*/
inline bool Model::stateFreqsFixed() const
	{
	return state_freq_fixed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `edgelen_hyperprior_fixed'.
*/
inline bool Model::edgeLenHyperParamFixed() const
	{
	return edgelen_hyperprior_fixed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `edge_lengths_fixed'.
*/
inline bool Model::edgeLengthsFixed() const
	{
	return edge_lengths_fixed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `edge_lengths_fixed' to true.
*/
inline void Model::fixEdgeLengths()
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
inline void Model::freeEdgeLengths()
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
|	Sets the data member `edgelen_hyperprior_fixed' to true.
*/
inline void Model::fixEdgeLenHyperprior()
	{
	edgelen_hyperprior_fixed = true;
	if (edgelen_hyper_param)	//TODO create Model::edgelen_hyper_param and assign it in createParameters
		edgelen_hyper_param->fixParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `edgelen_hyperprior_fixed' to false.
*/
inline void Model::freeEdgeLenHyperprior()
	{
	edgelen_hyperprior_fixed = false;
	if (edgelen_hyper_param)
		edgelen_hyper_param->freeParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `state_freq_fixed' to true. The fixParameter member function of all BaseFreqParam objects is 
|	either called immediately (if the `freq_params' vector is not empty) or is called in createParameters (when 
|	`freq_params' is built).
*/
inline void Model::fixStateFreqs()
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
|	Sets the data member `state_freq_fixed' to false. The freeParameter member function of all BaseFreqParam objects is 
|	called immediately if the `freq_params' vector is not empty.
*/
inline void Model::freeStateFreqs()
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
inline void Model::setFlexModel()
	{
	is_flex_model = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `flex_upper_rate_bound' to `new_upper_bound'. Rescales all existing unnormalized relative rate
|	parameter values in the vector `gamma_rates_unnorm' to accommodate the new upper bound. Assumes that 
|	`new_upper_bound' is greater than zero. Returns immediately if `new_upper_bound' is identical to current value of
|	`flex_upper_rate_bound'.
*/
inline void Model::setFlexRateUpperBound(double new_upper_bound)
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
inline void Model::setNumFlexSpacers(unsigned s)
	{
	num_flex_spacers = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `is_flex_model' to false. No FlexRateParam or FlexProbParam objects will be added in a 
|	subsequent call to the createParameters member function.
*/
inline void Model::setNotFlexModel()
	{
	is_flex_model = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `flex_probs_fixed' to true. The fixParameter member function of all FlexProbParam objects is 
|	either called immediately (if the `flex_prob_params' vector is not empty) or is called in createParameters (when 
|	`flex_prob_params' is built).
*/
inline void Model::fixFlexProbs()
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
inline void Model::freeFlexProbs()
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
inline void Model::fixFlexRates()
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
inline void Model::freeFlexRates()
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
inline void Model::setFLEXProbParamPrior(ProbDistShPtr d)
 	{
	flex_prob_param_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `flex_prob_param_prior'.
*/
inline ProbDistShPtr Model::getFLEXProbParamPrior()
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
inline void Model::setFlexRateUnnorm(
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
inline void Model::setFlexProbUnnorm(
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
inline void Model::normalizeRatesAndProbs(
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
|	Sets the data member `is_pinvar_model' to true. A subsequent call to the createParameters member function will
|	result in a PinvarParam being added to the list of updaters for this model.
*/
inline void Model::setPinvarModel()
	{
	is_pinvar_model = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `is_pinvar_model' to false. No PinvarParam will be added in a subsequent call to the 
|	createParameters member function.
*/
inline void Model::setNotPinvarModel()
	{
	is_pinvar_model = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `pinvar_fixed' to true. The fixParameter member function of the PinvarParam object is either 
|	called immediately (if `pinvar_param' is a valid pointer) or is called in createParameters (when `pinvar_param' is 
|	first assigned).
*/
inline void Model::fixPinvar()
	{
	pinvar_fixed = true;
	if (pinvar_param)
		pinvar_param->fixParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `pinvar_fixed' to false. The freeParameter member function of the PinvarParam object is called 
|	immediately if `pinvar_param' is a valid pointer.
*/
inline void Model::freePinvar()
	{
	pinvar_fixed = false;
	if (pinvar_param)
		pinvar_param->freeParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `gamma_shape_prior'.
*/
inline ProbDistShPtr Model::getDiscreteGammaShapePrior()
 	{
	return gamma_shape_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of data member `gamma_shape_prior'.
*/
inline void Model::setDiscreteGammaShapePrior(ProbDistShPtr d)
 	{
	gamma_shape_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `pinvar_prior'.
*/
inline ProbDistShPtr Model::getPinvarPrior()
 	{
	return pinvar_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `pinvar_prior'.
*/
inline void Model::setPinvarPrior(ProbDistShPtr d)
 	{
	pinvar_prior = d;
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
inline void Model::recalcRatesAndProbs(
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
			//@POL use std::transform algorithm with a boost::lambda here?
			for (unsigned i = 0; i < num_gamma_rates; ++i)
				{
				*rates_iter++ = (*gamma_rates_iter++)/one_minus_pinvar;
				*probs_iter++ = (*gamma_probs_iter++)*one_minus_pinvar;
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
inline void Model::recalcGammaRatesAndBoundaries(std::vector<double> & rates, std::vector<double> & boundaries) const
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
|	Constructor sets `num_states' data member to 4.
*/
inline JC::JC()
  : Model(4)
	{
	state_repr.reserve(4);
	state_repr.push_back("A");
	state_repr.push_back("C");
	state_repr.push_back("G");
	state_repr.push_back("T");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string indicating the name of this model, namely "JC69", "JC69+G, "JC69+I" or "JC69+G+I".
*/
inline std::string JC::getModelName() const
	{
	std::string s = "JC69";
	if (is_pinvar_model)
		s += "+I";
	if (num_gamma_rates > 1)
		s += "+G";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears `state_list' and `state_list_pos', then builds up a 4-state version of these vectors suitable for simulated
|	nucleotide data (i.e. no ambiguities).
*/
inline void Model::buildStateList(
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
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by 
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for 
|	compatability with MrBayes). The JC model provides additional column labels for the gamma shape parameter (if the 
|	number of rates is greater than 1) and the pinvar parameter (if an invariable sites model is being used).
*/
inline std::string JC::paramHeader() const
	{
	std::string s = std::string("Gen\tLnL\tTL");
	if (is_flex_model)
		{
		s += "\tncat";
		//unsigned i;
		//for ( i = 0; i < num_gamma_rates; ++i)
		//	{
		//	s += "\trate";
		//	s += i;
		//	}
		//for ( i = 0; i < num_gamma_rates; ++i)
		//	{
		//	s += "\trateprob";
		//	s += i;
		//	}
		}
	else if (num_gamma_rates > 1)
		{
		s += "\tshape";
		}
	if (is_pinvar_model)
		s += "\tpinvar";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The only parameters in the JC model are edge lengths, which are not reported by this function, but we must provide
|	an override of the pure virtual base class version.
*/
inline std::string JC::paramReport() const
	{
	std::string s;
	if (is_flex_model)
		{
		s += str(boost::format("\t%d") % num_gamma_rates);
		//std::vector<double> rates(num_gamma_rates, 0.0);
		//std::vector<double> probs(num_gamma_rates, 0.0);
		//normalizeRatesAndProbs(rates, probs);
		//unsigned i;
		//for ( i = 0; i < num_gamma_rates; ++i)
		//	{
		//	s += str(boost::format("\t%.5f") % rates[i]);
		//	}
		//for ( i = 0; i < num_gamma_rates; ++i)
		//	{
		//	s += str(boost::format("\t%.5f") % probs[i]);
		//	}
		}
	else if (num_gamma_rates > 1)
		{
		s += str(boost::format("\t%.5f") % gamma_shape);
		}
	if (is_pinvar_model)
		s += str(boost::format("\t%.5f") % pinvar);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string in the `state_repr' vector corresponding to `state'. If supplied `state' is negative, returns
|	the gap state ("-"). If `state' equals or exceeds `num_states', returns "?".
*/
inline std::string Model::lookupStateRepr(
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
|	Constructor sets `num_states' data member to 4, the base frequencies `piA', `piC', `piG' and `piT' to 0.25, and `kappa'
|	to 1.0.
*/
inline HKY::HKY()
  : Model(4), kappa(1.0), kappa_fixed(false)
	{
	state_repr.reserve(4);
	state_repr.push_back("A");
	state_repr.push_back("C");
	state_repr.push_back("G");
	state_repr.push_back("T");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string indicating the name of this model, namely "HKY85", "HKY85+G", "HKY85+I" or "HKY85+G+I".
*/
inline std::string HKY::getModelName() const
	{
	std::string s = "HKY85";
	if (num_gamma_rates > 1)
		s += "+G";
	if (is_pinvar_model)
		s += "+I";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Model::createParameters to create the edge length parameters, the edge length hyperparameter, and any 
|	parameters related to rate heterogeneity. This function then adds additional HKY85-specific parameters to the 
|	supplied `parameters' vector. This incudes the four base frequencies as well as the transition/transversion rate 
|	ratio kappa.
*/
inline void	HKY::createParameters(
  TreeShPtr t,								/**< is the tree (the nodes of which are needed for creating edge length parameters) */
  MCMCUpdaterVect & edgelens,				/**< is the vector of edge length parameters to fill */
  MCMCUpdaterShPtr & edgelen_hyperparam,	/**< is the edge length hyperparameter */
  MCMCUpdaterVect & parameters,				/**< is the vector of model-specific parameters to fill */
  bool separate_edgelens) const				/**< specifies (if true) that each edge should have its own parameter or (if false) that one edge length master parameter should be created */
	{
	Model::createParameters(t, edgelens, edgelen_hyperparam, parameters, separate_edgelens);

	PHYCAS_ASSERT(!kappa_param);
	kappa_param = MCMCUpdaterShPtr(new KappaParam());
	kappa_param->setName("trs/trv rate ratio");
	kappa_param->setStartingValue(4.0);
	kappa_param->setTree(t);
	kappa_param->setPrior(kappa_prior);
	if (kappa_fixed)
		kappa_param->fixParameter();
	parameters.push_back(kappa_param);

	PHYCAS_ASSERT(freq_params.empty());

	MCMCUpdaterShPtr freqA_param = MCMCUpdaterShPtr(new BaseFreqParam(0));
	freqA_param->setName("base freq. A");
	freqA_param->setTree(t);
	freqA_param->setStartingValue(1.0);
	freqA_param->setPrior(freq_param_prior);
	if (state_freq_fixed)
		freqA_param->fixParameter();
	parameters.push_back(freqA_param);
	freq_params.push_back(freqA_param);

	MCMCUpdaterShPtr freqC_param = MCMCUpdaterShPtr(new BaseFreqParam(1));
	freqC_param->setName("base freq. C");
	freqC_param->setTree(t);
	freqC_param->setStartingValue(1.0);
	freqC_param->setPrior(freq_param_prior);
	if (state_freq_fixed)
		freqC_param->fixParameter();
	parameters.push_back(freqC_param);
	freq_params.push_back(freqC_param);

	MCMCUpdaterShPtr freqG_param = MCMCUpdaterShPtr(new BaseFreqParam(2));
	freqG_param->setName("base freq. G");
	freqG_param->setTree(t);
	freqG_param->setStartingValue(1.0);
	freqG_param->setPrior(freq_param_prior);
	if (state_freq_fixed)
		freqG_param->fixParameter();
	parameters.push_back(freqG_param);
	freq_params.push_back(freqG_param);

	MCMCUpdaterShPtr freqT_param = MCMCUpdaterShPtr(new BaseFreqParam(3));
	freqT_param->setName("base freq. T");
	freqT_param->setTree(t);
	freqT_param->setStartingValue(1.0);
	freqT_param->setPrior(freq_param_prior);
	if (state_freq_fixed)
		freqT_param->fixParameter();
	parameters.push_back(freqT_param);
	freq_params.push_back(freqT_param);
	}


/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by 
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for 
|	compatability with MrBayes). The HKY model provide additional columns for kappa, the base frequencies, the gamma 
|	shape parameter (if the number of rates is greater than 1) and the pinvar parameter (if an invariable sites model 
|	is being used)
*/
inline std::string HKY::paramHeader() const
	{
	std::string s = std::string("Gen\tLnL\tTL\tkappa\tfreqA\tfreqC\tfreqG\tfreqT");
	if (is_flex_model)
		{
		s += "\tncat";

		//std::ofstream ratef("flex_rates.txt");
		//ratef << str(boost::format("%12s\t%12s\t%12s\t%12s\n") % "i" % "n" % "rate" % "prob");
		//ratef.close();

		//unsigned i;
		//for ( i = 0; i < num_gamma_rates; ++i)
		//	{
		//	s += str(boost::format("\trate%d") % i);
		//	}
		//for ( i = 0; i < num_gamma_rates; ++i)
		//	{
		//	s += str(boost::format("\trateprob%d") % i);
		//	}
		}
	else if (num_gamma_rates > 1)
		{
		s += "\tshape";
		}
	if (is_pinvar_model)
		s += "\tpinvar";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides the pure virtual base class version to generate a string of tab-separated values of model-specific 
|	parameters suitable for saving in a sampled parameter file (e.g. like the .p files saved by MrBayes).
*/
inline std::string HKY::paramReport() const
	{
	std::string s = str(boost::format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f") % kappa % state_freqs[0] % state_freqs[1] % state_freqs[2] % state_freqs[3]);;
	if (is_flex_model)
		{
		s += str(boost::format("\t%d") % num_gamma_rates);

		//std::ofstream ratef("flex_rates.txt", std::ios::out | std::ios::app);
		//for (unsigned i = 0; i < num_gamma_rates; ++i)
		//	ratef << str(boost::format("%12d\t%12d\t%12.5f\t%12.5f\n") % i % num_gamma_rates % gamma_rates_unnorm[i] % gamma_rate_probs[i]);
		//ratef.close();

		//std::vector<double> rates(num_gamma_rates, 0.0);
		//std::vector<double> probs(num_gamma_rates, 0.0);
		//normalizeRatesAndProbs(rates, probs);
		//unsigned i;
		//for ( i = 0; i < num_gamma_rates; ++i)
		//	{
		//	s += str(boost::format("\t%.5f") % rates[i]);
		//	}
		//for ( i = 0; i < num_gamma_rates; ++i)
		//	{
		//	s += str(boost::format("\t%.5f") % probs[i]);
		//	}
		}
	else if (num_gamma_rates > 1)
		{
		s += str(boost::format("\t%.5f") % gamma_shape);
		}
	if (is_pinvar_model)
		s += str(boost::format("\t%.5f") % pinvar);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `kappa_fixed' to true. The fixParameter member function of the KappaParam object is either 
|	called immediately (if `kappa_param' is a valid pointer) or is called in createParameters (when `kappa_param' is 
|	first assigned).
*/
inline void HKY::fixKappa()
	{
	kappa_fixed = true;
	if (kappa_param)
		kappa_param->fixParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `kappa_fixed' to false. The freeParameter member function of the KappaParam object is called 
|	immediately if `kappa_param' is a valid pointer.
*/
inline void HKY::freeKappa()
	{
	kappa_fixed = false;
	if (kappa_param)
		kappa_param->freeParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `kappa'.
*/
inline double HKY::getKappa()
 	{
	return kappa;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa' to supplied value `k'. Throws XLikelihood exception if `k' is less than or equal to 0.0.
*/
inline void HKY::setKappa(double k)
 	{
	if (k <= 0.0)
		throw XLikelihood();
	kappa = k;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `kappa_prior'.
*/
inline ProbDistShPtr HKY::getKappaPrior()
 	{
	return kappa_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
inline void HKY::setKappaPrior(ProbDistShPtr d)
 	{
	kappa_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_param_prior'.
*/
inline ProbDistShPtr HKY::getBaseFreqParamPrior()
 	{
	return freq_param_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_param_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
inline void HKY::setBaseFreqParamPrior(ProbDistShPtr d)
 	{
	freq_param_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa' using supplied value `tratio', which is the transition/transversion ratio, not the ts/tv rate ratio. 
|	The quantity `kappa' is related to `tratio' as follows (where piA means `state_freqs'[0], piC means 
|	`state_freqs'[1], piC means	`state_freqs'[2] and piT means `state_freqs'[3]):
|>
|	         tratio (piA + piG) (piC + piT)

|	kappa = -------------------------------

|	              (piA piG + piC piT)

|>
|	Throws XLikelihood exception if `tratio' is less than or equal to 0.0.
*/
inline void HKY::setKappaFromTRatio(double tratio)
 	{
	if (tratio <= 0.0)
		throw XLikelihood();
	double numerator   = tratio*(state_freqs[0] + state_freqs[2])*(state_freqs[1] + state_freqs[3]);

	double denominator = ((state_freqs[0]*state_freqs[2]) + (state_freqs[1]*state_freqs[3]));

	PHYCAS_ASSERT(denominator != 0.0);

	kappa = numerator/denominator;
	}

/*----------------------------------------------------------------------------------------------------------------------

|	Calculates the transition/transversion ratio given the transition/transversion rate ratio `kappa' and the relative

|	base frequencies, which are stored in the `state_freqs' vector. Here are the details of the calculation (for brevity,

|	`state_freqs'[0] has been symbolized piA, `state_freqs'[1] by piC, `state_freqs'[2] by piG and `state_freqs'[3] by piT):

|>

|	Parameters: b = transversion rate, k = kappa

|	Pr(any transition | dt)   = Pr(AG) + Pr(CT) + Pr(GA) + Pr(TC) 

|	                          = (piA piG k b dt) + (piC piT k b dt) + (piG piA k b dt) + (piT piC k b dt)

|	                          = 2 k b dt (piA piG + piC piT)

|	

|	Pr(any transversion | dt) = Pr(AC) + Pr(AT) + Pr(CA) + Pr(CG) + Pr(GC) + Pr(GT) + Pr(TA) + Pr(TG)

|	                          = (piA piC b dt) + (piA piT b dt) + (piC piA b dt) + (piC piG b dt)

|	                            + (piG piC b dt) + (piG piT b dt) + (piT piA b dt) + (piT piG b dt)

|	                          = 2 b dt (piA + piG) (piC + piT)

|	

|	          2 k b dt (piA piG + piC piT)     k (piA piG + piC piT)

|	TRatio = ------------------------------ = -----------------------

|	         2 b dt (piA + piG) (piC + piT)   (piA + piG) (piC + piT)

|<

*/

inline double HKY::calcTRatio()

	{

	double numerator   = kappa*((state_freqs[0]*state_freqs[2]) + (state_freqs[1]*state_freqs[3]));

	double denominator = (state_freqs[0] + state_freqs[2])*(state_freqs[1] + state_freqs[3]);

	PHYCAS_ASSERT(denominator != 0.0);

	double tratio = numerator/denominator;

	return tratio;

	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `num_states' data member to 4, the base frequencies `piA', `piC', `piG' and `piT' to 0.25, and all
|	six relative rates to 1.0.
*/
inline GTR::GTR()
  : Model(4), rel_rates_fixed(false)
	{
	state_repr.reserve(4);
	state_repr.push_back("A");
	state_repr.push_back("C");
	state_repr.push_back("G");
	state_repr.push_back("T");

	rel_rates.assign(6, 1.0);

	q_matrix.setRelativeRates(rel_rates);
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string indicating the name of this model, namely "GTR", "GTR+G", "GTR+I" or "GTR+G+I".
*/
inline std::string GTR::getModelName() const
	{
	std::string s = "GTR";
	if (num_gamma_rates > 1)
		s += "+G";
	if (is_pinvar_model)
		s += "+I";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Model::createParameters to create the edge length parameters, the edge length hyperparameter, and any 
|	parameters related to rate heterogeneity. This function then adds additional GTR-specific parameters to the 
|	supplied `parameters' vector. This incudes the four base frequencies and six relative rate parameters.
*/
inline void	GTR::createParameters(
  TreeShPtr t,								/**< is the tree (the nodes of which are needed for creating edge length parameters) */
  MCMCUpdaterVect & edgelens,				/**< is the vector of edge length parameters to fill */
  MCMCUpdaterShPtr & edgelen_hyperparam,	/**< is the edge length hyperparameter */
  MCMCUpdaterVect & parameters,				/**< is the vector of model-specific parameters to fill */
  bool separate_edgelens) const				/**< specifies (if true) that each edge should have its own parameter or (if false) that one edge length master parameter should be created */
	{
	Model::createParameters(t, edgelens, edgelen_hyperparam, parameters, separate_edgelens);

	PHYCAS_ASSERT(rel_rate_params.empty());

	MCMCUpdaterShPtr rAC_param = MCMCUpdaterShPtr(new GTRRateParam(0));
	rAC_param->setName("rAC");
	rAC_param->setTree(t);
	rAC_param->setStartingValue(1.0);
	rAC_param->setPrior(rel_rate_prior);
	if (rel_rates_fixed)
		rAC_param->fixParameter();
	parameters.push_back(rAC_param);
	rel_rate_params.push_back(rAC_param);

	MCMCUpdaterShPtr rAG_param = MCMCUpdaterShPtr(new GTRRateParam(1));
	rAG_param->setName("rAG");
	rAG_param->setTree(t);
	rAG_param->setStartingValue(4.0);
	rAG_param->setPrior(rel_rate_prior);
	if (rel_rates_fixed)
		rAG_param->fixParameter();
	parameters.push_back(rAG_param);
	rel_rate_params.push_back(rAG_param);

	MCMCUpdaterShPtr rAT_param = MCMCUpdaterShPtr(new GTRRateParam(2));
	rAT_param->setName("rAT");
	rAT_param->setTree(t);
	rAT_param->setStartingValue(1.0);
	rAT_param->setPrior(rel_rate_prior);
	if (rel_rates_fixed)
		rAT_param->fixParameter();
	parameters.push_back(rAT_param);
	rel_rate_params.push_back(rAT_param);

	MCMCUpdaterShPtr rCG_param = MCMCUpdaterShPtr(new GTRRateParam(3));
	rCG_param->setName("rCG");
	rCG_param->setTree(t);
	rCG_param->setStartingValue(1.0);
	rCG_param->setPrior(rel_rate_prior);
	if (rel_rates_fixed)
		rCG_param->fixParameter();
	parameters.push_back(rCG_param);
	rel_rate_params.push_back(rCG_param);

	MCMCUpdaterShPtr rCT_param = MCMCUpdaterShPtr(new GTRRateParam(4));
	rCT_param->setName("rCT");
	rCT_param->setTree(t);
	rCT_param->setStartingValue(4.0);
	rCT_param->setPrior(rel_rate_prior);
	if (rel_rates_fixed)
		rCT_param->fixParameter();
	parameters.push_back(rCT_param);
	rel_rate_params.push_back(rCT_param);

	MCMCUpdaterShPtr rGT_param = MCMCUpdaterShPtr(new GTRRateParam(5));
	rGT_param->setName("rGT");
	rGT_param->setStartingValue(1.0);
	rGT_param->setTree(t);
	rGT_param->setPrior(rel_rate_prior);
	if (rel_rates_fixed)
		rGT_param->fixParameter();
	parameters.push_back(rGT_param);
	rel_rate_params.push_back(rGT_param);

	PHYCAS_ASSERT(freq_params.empty());

	MCMCUpdaterShPtr freqA_param = MCMCUpdaterShPtr(new BaseFreqParam(0));
	freqA_param->setName("base freq. A");
	freqA_param->setTree(t);
	freqA_param->setStartingValue(1.0);
	freqA_param->setPrior(freq_param_prior);
	if (state_freq_fixed)
		freqA_param->fixParameter();
	parameters.push_back(freqA_param);
	freq_params.push_back(freqA_param);

	MCMCUpdaterShPtr freqC_param = MCMCUpdaterShPtr(new BaseFreqParam(1));
	freqC_param->setName("base freq. C");
	freqC_param->setTree(t);
	freqC_param->setStartingValue(1.0);
	freqC_param->setPrior(freq_param_prior);
	if (state_freq_fixed)
		freqC_param->fixParameter();
	parameters.push_back(freqC_param);
	freq_params.push_back(freqC_param);

	MCMCUpdaterShPtr freqG_param = MCMCUpdaterShPtr(new BaseFreqParam(2));
	freqG_param->setName("base freq. G");
	freqG_param->setTree(t);
	freqG_param->setStartingValue(1.0);
	freqG_param->setPrior(freq_param_prior);
	if (state_freq_fixed)
		freqG_param->fixParameter();
	parameters.push_back(freqG_param);
	freq_params.push_back(freqG_param);

	MCMCUpdaterShPtr freqT_param = MCMCUpdaterShPtr(new BaseFreqParam(3));
	freqT_param->setName("base freq. T");
	freqT_param->setTree(t);
	freqT_param->setStartingValue(1.0);
	freqT_param->setPrior(freq_param_prior);
	if (state_freq_fixed)
		freqT_param->fixParameter();
	parameters.push_back(freqT_param);
	freq_params.push_back(freqT_param);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the transition probability matrix given an edge length. Overrides the pure virtual function inherited from 
|	the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
inline void GTR::calcPMat(double * * pMat, double edgeLength) const
	{
	q_matrix.recalcPMat(pMat, edgeLength);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `rel_rates_fixed' to true. The fixParameter member function of all GTRRateParam objects is 
|	either called immediately (if the `rel_rate_params' vector is not empty) or is called in createParameters (when 
|	`rel_rate_params' is built).
*/
inline void GTR::fixRelRates()
	{
	rel_rates_fixed = true;
	if (!rel_rate_params.empty())
		{
#if 0
		//@POL good place to use std::for_each with a boost::lambda functor?
		for (MCMCUpdaterVect::iterator it = rel_rate_params.begin(); it != rel_rate_params.end(); ++it)
			(*it)->fixParameter();
#else
		std::for_each(rel_rate_params.begin(), rel_rate_params.end(), boost::lambda::bind(&MCMCUpdater::fixParameter, *boost::lambda::_1));
#endif
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `rel_rates_fixed' to false. The freeParameter member function of all GTRRateParam objects is 
|	called immediately if the `rel_rate_params' vector is not empty.
*/
inline void GTR::freeRelRates()
	{
	rel_rates_fixed = false;
	if (!rel_rate_params.empty())
		{
#if 0
		//@POL good place to use std::for_each with a boost::lambda functor?
		for (MCMCUpdaterVect::iterator it = rel_rate_params.begin(); it != rel_rate_params.end(); ++it)
			(*it)->freeParameter();
#else
		std::for_each(rel_rate_params.begin(), rel_rate_params.end(), boost::lambda::bind(&MCMCUpdater::freeParameter, *boost::lambda::_1));
#endif
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a copy of the `rel_rates' vector. 
*/
inline std::vector<double> GTR::getRelRates()
 	{
	return rel_rates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Copies the supplied vector `rates' into the data member vector `rel_rates'. Throws XLikelihood exception if any of
|	the values of `rates' is less than 0.0, or if the number of elements in `rates' is not equal to 6. Assumes that
|	`rel_rates' has 6 elements.
*/
inline void GTR::setRelRates(const std::vector<double> & rates)
 	{
	PHYCAS_ASSERT(rel_rates.size() == 6);

	// Ensure that we are not about to set any negative relative rates
	if (std::find_if(rates.begin(), rates.end(), boost::lambda::_1 < 0.0) != rates.end())
		throw XLikelihood("GTR relative rates cannot be less than 0.0");

	// Ensure that user supplied the correct number of relative rates
	if (rates.size() != 6)
		throw XLikelihood(str(boost::format("The number of GTR relative rates is 6, but %d were supplied") % rates.size()));

	// Ok, seems safe to copy now
	std::copy(rates.begin(), rates.end(), rel_rates.begin());

	// Make sure QMatrix knows about this change in relative rates
	q_matrix.setRelativeRates(rel_rates);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets one of the GTR relative rate parameters. Throws an XLikelihood exception if 
|	`param_index' is not less than 6 or if `value' is negative. Assumes that `rel_rates' has length 6.
*/
inline void GTR::setRelRateUnnorm(
  unsigned param_index,		/**< the 0-based index into the `rel_rates' vector of the element to modify */
  double value)				/**< the new value of `rel_rates'[`param_index'] */
	{
	PHYCAS_ASSERT(rel_rates.size() == 6);

	// Ensure that user supplied a correct index
	if (param_index >= 6)
		throw XLikelihood(str(boost::format("Index of a GTR rate parameter must be less than 6, but %d was supplied") % param_index));

	// Ensure that we are not about to set a negative relative rate
	if (value < 0.0)
		throw XLikelihood(str(boost::format("GTR rate parameters must be greater than or equal to 0.0, but the value %f was supplied") % value));

	// Ok, seems safe to overwrite value now
	rel_rates[param_index] = value;	//@POL do we need rel_rates if q_matrix is maintaining a copy?

	// Make sure QMatrix knows about this change in relative rates
	q_matrix.setRelativeRates(rel_rates);	//@POL could speed up by just changing one rate rather than copying all 6 each time
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `rel_rate_prior'.
*/
inline ProbDistShPtr GTR::getRelRatePrior()
 	{
	return rel_rate_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `rel_rate_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
inline void GTR::setRelRatePrior(ProbDistShPtr d)
 	{
	rel_rate_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of Model base class function that sets all four nucleotide frequency parameters. The base class version is
|	called to do most of the work, and this function is responsible only for ensuring that the `q_matrix' data member
|	knows about the change in base frequencies.
*/
inline void GTR::setNucleotideFreqs(
  double freqA,				/**< the new frequency of base A */
  double freqC,				/**< the new frequency of base C */
  double freqG,				/**< the new frequency of base G */
  double freqT)				/**< the new frequency of base T/U */
	{
	Model::setNucleotideFreqs(freqA, freqC, freqG, freqT);
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of Model base class function that sets all four nucleotide frequency parameters to 0.25. The base class 
|	version is called to do most of the work, and this function is responsible only for ensuring that the `q_matrix' 
|	data member knows about the change in base frequencies.
*/
inline void GTR::setAllFreqsEqual()
	{
	Model::setAllFreqsEqual();
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of Model base class function that sets one of the state frequency parameters (the unnormalized values that
|	determine the values in the `state_freqs' vector when normalized) in the data member vector `state_freq_unnorm'. 
|	The base class version is called to do most of the work, and this function is responsible only for ensuring that the
|	`q_matrix' data member knows about the change in base frequencies.
*/
inline void GTR::setStateFreqUnnorm(
  unsigned param_index,		/**< the 0-based index into the `state_freq_unnorm' vector of the element to modify */
  double value)				/**< the new value of `state_freq_unnorm'[`param_index'] */
	{
	Model::setStateFreqUnnorm(param_index, value);
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_param_prior'.
*/
inline ProbDistShPtr GTR::getBaseFreqParamPrior()
 	{
	return freq_param_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_param_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
inline void GTR::setBaseFreqParamPrior(ProbDistShPtr d)
 	{
	freq_param_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by 
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for 
|	compatability with MrBayes). The GTR model provide additional columns for the six relative rates, the base 
|	frequencies, the gamma shape parameter (if the number of rates is greater than 1) and the pinvar parameter (if
|	an invariable sites model is being used).
*/
inline std::string GTR::paramHeader() const
	{
	std::string s = std::string("Gen\tLnL\tTL\trAC\trAG\trAT\trCG\trCT\trGT\tfreqA\tfreqC\tfreqG\tfreqT");
	if (is_flex_model)
		{
		s += "\tncat";
		}
	else if (num_gamma_rates > 1)
		{
		s += "\tshape";
		}
	if (is_pinvar_model)
		s += "\tpinvar";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides the pure virtual base class version to generate a string of tab-separated values of model-specific 
|	parameters suitable for saving in a sampled parameter file (e.g. like the .p files saved by MrBayes). Assumes that
|	the `state_freqs' vector has length 4 and the `rel_rates' vector has length 6.
*/
inline std::string GTR::paramReport() const
	{
	PHYCAS_ASSERT(rel_rates.size() == 6);
	PHYCAS_ASSERT(state_freqs.size() == 4);
	std::string s = str(boost::format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f") % rel_rates[0] % rel_rates[1] % rel_rates[2] % rel_rates[3] % rel_rates[4] % rel_rates[5]);
	s += str(boost::format("\t%.5f\t%.5f\t%.5f\t%.5f") % state_freqs[0] % state_freqs[1] % state_freqs[2] % state_freqs[3]);
	if (is_flex_model)
		{
		s += str(boost::format("\t%d") % num_gamma_rates);
		}
	else if (num_gamma_rates > 1)
		{
		s += str(boost::format("\t%.5f") % gamma_shape);
		}
	if (is_pinvar_model)
		s += str(boost::format("\t%.5f") % pinvar);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------

|	Calculates the transition/transversion ratio given the six relative rates stored in the `rel_rates' vector, and the 

|	relative base frequencies, which are stored in the `state_freqs' vector. Here are the details of the calculation 

|	(for brevity, `state_freqs'[0] has been symbolized piA, `state_freqs'[1] by piC, `state_freqs'[2] by piG and 

|	`state_freqs'[3] by piT, `rel_rates'[0] by rAC, `rel_rates'[1] by rAG `rel_rates'[2] by rAT, `rel_rates'[3] by rCG, 

|	`rel_rates'[4] by rCT and `rel_rates'[5] by rGT):

|>

|	Pr(any transition | dt)   = Pr(AG) + Pr(CT) + Pr(GA) + Pr(TC) 

|	                          = (piA piG rAG dt) + (piC piT rCT dt) + (piG piA rAG dt) + (piT piC rCT dt)

|	                          = 2 dt (piA piG rAG + piT piC rCT dt)

|	

|	Pr(any transversion | dt) = Pr(AC) + Pr(AT) + Pr(CA) + Pr(CG) + Pr(GC) + Pr(GT) + Pr(TA) + Pr(TG)

|	                          = (piA piC rAC dt) + (piA piT rAT dt) + (piC piA rAC dt) + (piC piG rCG dt)

|	                            + (piG piC rCG dt) + (piG piT rGT dt) + (piT piA rAT dt) + (piT piG rGT dt)

|	                          = 2 dt (piA piC rAC + piA piT rAT + piC piG rCG +  + piG piT rGT)

|	

|	                      2 dt (piA piG rAG + piT piC rCT)

|	TRatio = ------------------------------------------------------------

|	         2 dt (piA piC rAC + piA piT rAT + piC piG rCG + piG piT rGT)

|	

|	                     piA piG rAG + piT piC rCT

|	       = -----------------------------------------------------

|	         piA piC rAC + piA piT rAT + piC piG rCG + piG piT rGT

|<

*/

inline double GTR::calcTRatio()

	{

	double numerator = state_freqs[0]*state_freqs[2]*rel_rates[1] + state_freqs[3]*state_freqs[1]*rel_rates[4];

	double denominator = state_freqs[0]*state_freqs[1]*rel_rates[0] + state_freqs[0]*state_freqs[3]*rel_rates[2] + state_freqs[1]*state_freqs[2]*rel_rates[3] + state_freqs[2]*state_freqs[3]*rel_rates[5];

	PHYCAS_ASSERT(denominator != 0.0);

	double tratio = numerator/denominator;

	return tratio;

	}

}	// namespace phycas

#endif

