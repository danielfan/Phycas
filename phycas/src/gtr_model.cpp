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

//#include <cmath>
//#include <iostream>
#include "phycas/src/gtr_model.hpp"
//#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
//#	include <boost/python/numeric.hpp>
//#	include "phycas/src/thirdparty/num_util/num_util.h"
//#endif
//using std::cout;
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `num_states' data member to 4, the base frequencies `piA', `piC', `piG' and `piT' to 0.25, and all
|	six relative rates to 1.0.
*/
GTR::GTR()
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
std::string GTR::getModelName() const
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
void	GTR::createParameters(
  TreeShPtr t,								/**< is the tree (the nodes of which are needed for creating edge length parameters) */
  MCMCUpdaterVect & edgelens,				/**< is the vector of edge length parameters to fill */
  MCMCUpdaterVect & edgelen_hyperparams,	/**< is the edge length hyperparameter */
  MCMCUpdaterVect & parameters,				/**< is the vector of model-specific parameters to fill */
  int subset_pos) 							/**< if 0 (first subset) or -1 (only subset), edge length parameters and hyperparams will be added; otherwise, the `edgelens' and `edgelen_hyperparams' vectors returned will be empty */
	{
	Model::createParameters(t, edgelens, edgelen_hyperparams, parameters, subset_pos);

	PHYCAS_ASSERT(rel_rate_params.empty());
    PHYCAS_ASSERT(rel_rate_param_prior || rel_rate_prior);
	relrate_name.clear();
	freq_name.clear();
	if (subset_pos < 0)
		{
		relrate_name.push_back("rAC");
		relrate_name.push_back("rAG");
		relrate_name.push_back("rAT");
		relrate_name.push_back("rCG");
		relrate_name.push_back("rCT");
		relrate_name.push_back("rGT");
		freq_name.push_back("freqA");
		freq_name.push_back("freqC");
		freq_name.push_back("freqG");
		freq_name.push_back("freqT");
		}
	else 
		{
		unsigned m = subset_pos + 1;
		relrate_name.push_back(boost::str(boost::format("rAC_%d") % m));
		relrate_name.push_back(boost::str(boost::format("rAG_%d") % m));
		relrate_name.push_back(boost::str(boost::format("rAT_%d") % m));
		relrate_name.push_back(boost::str(boost::format("rCG_%d") % m));
		relrate_name.push_back(boost::str(boost::format("rCT_%d") % m));
		relrate_name.push_back(boost::str(boost::format("rGT_%d") % m));
		freq_name.push_back(boost::str(boost::format("freqA_%d") % m));
		freq_name.push_back(boost::str(boost::format("freqC_%d") % m));
		freq_name.push_back(boost::str(boost::format("freqG_%d") % m));
		freq_name.push_back(boost::str(boost::format("freqT_%d") % m));
		}
    if (rel_rate_param_prior)
        {
        // The fact that rel_rate_param_prior points to something means that the user
        // chose to update the relative rates separately using slice sampling. If the 
        // user had chosen to update the relative rates jointly (using the RelRateMove
        // Metropolis-Hastings move), then relrate_prior would be set and rel_rate_param_prior
        // would be empty.
    
        MCMCUpdaterShPtr rAC_param = MCMCUpdaterShPtr(new GTRRateParam(0));
		rAC_param->setName(relrate_name[0]);
        rAC_param->setTree(t);
        rAC_param->setStartingValue(1.0);
        rAC_param->setPrior(rel_rate_param_prior);
        if (rel_rates_fixed)
            rAC_param->fixParameter();
        parameters.push_back(rAC_param);
        rel_rate_params.push_back(rAC_param);
    
        MCMCUpdaterShPtr rAG_param = MCMCUpdaterShPtr(new GTRRateParam(1));
		rAG_param->setName(relrate_name[1]);
        rAG_param->setTree(t);
        rAG_param->setStartingValue(4.0);
        rAG_param->setPrior(rel_rate_param_prior);
        if (rel_rates_fixed)
            rAG_param->fixParameter();
        parameters.push_back(rAG_param);
        rel_rate_params.push_back(rAG_param);
    
        MCMCUpdaterShPtr rAT_param = MCMCUpdaterShPtr(new GTRRateParam(2));
		rAT_param->setName(relrate_name[2]);
        rAT_param->setTree(t);
        rAT_param->setStartingValue(1.0);
        rAT_param->setPrior(rel_rate_param_prior);
        if (rel_rates_fixed)
            rAT_param->fixParameter();
        parameters.push_back(rAT_param);
        rel_rate_params.push_back(rAT_param);
    
        MCMCUpdaterShPtr rCG_param = MCMCUpdaterShPtr(new GTRRateParam(3));
		rCG_param->setName(relrate_name[3]);
        rCG_param->setTree(t);
        rCG_param->setStartingValue(1.0);
        rCG_param->setPrior(rel_rate_param_prior);
        if (rel_rates_fixed)
            rCG_param->fixParameter();
        parameters.push_back(rCG_param);
        rel_rate_params.push_back(rCG_param);
    
        MCMCUpdaterShPtr rCT_param = MCMCUpdaterShPtr(new GTRRateParam(4));
		rCT_param->setName(relrate_name[4]);
        rCT_param->setTree(t);
        rCT_param->setStartingValue(4.0);
        rCT_param->setPrior(rel_rate_param_prior);
        if (rel_rates_fixed)
            rCT_param->fixParameter();
        parameters.push_back(rCT_param);
        rel_rate_params.push_back(rCT_param);
    
        MCMCUpdaterShPtr rGT_param = MCMCUpdaterShPtr(new GTRRateParam(5));
		rGT_param->setName(relrate_name[5]);
        rGT_param->setStartingValue(1.0);
        rGT_param->setTree(t);
        rGT_param->setPrior(rel_rate_param_prior);
        if (rel_rates_fixed)
            rGT_param->fixParameter();
        parameters.push_back(rGT_param);
        rel_rate_params.push_back(rGT_param);
        }

	PHYCAS_ASSERT(freq_params.empty());
    PHYCAS_ASSERT(freq_param_prior || freq_prior);
    if (freq_param_prior)
        {
        // Only add frequency parameters if freqs will be updated separately
        // The other option is to update the frequencies jointly using the 
        // StateFreqMove Metropolis-Hastings move (in which case freq_prior
        // will be set and freq_param_prior will be empty)

	    MCMCUpdaterShPtr freqA_param = MCMCUpdaterShPtr(new StateFreqParam(0));
		freqA_param->setName(freq_name[0]);
	    freqA_param->setTree(t);
	    freqA_param->setStartingValue(1.0);
	    freqA_param->setPrior(freq_param_prior);
	    if (state_freq_fixed)
		    freqA_param->fixParameter();
	    parameters.push_back(freqA_param);
	    freq_params.push_back(freqA_param);

	    MCMCUpdaterShPtr freqC_param = MCMCUpdaterShPtr(new StateFreqParam(1));
		freqC_param->setName(freq_name[1]);
	    freqC_param->setTree(t);
	    freqC_param->setStartingValue(1.0);
	    freqC_param->setPrior(freq_param_prior);
	    if (state_freq_fixed)
		    freqC_param->fixParameter();
	    parameters.push_back(freqC_param);
	    freq_params.push_back(freqC_param);

	    MCMCUpdaterShPtr freqG_param = MCMCUpdaterShPtr(new StateFreqParam(2));
		freqG_param->setName(freq_name[2]);
	    freqG_param->setTree(t);
	    freqG_param->setStartingValue(1.0);
	    freqG_param->setPrior(freq_param_prior);
	    if (state_freq_fixed)
		    freqG_param->fixParameter();
	    parameters.push_back(freqG_param);
	    freq_params.push_back(freqG_param);

	    MCMCUpdaterShPtr freqT_param = MCMCUpdaterShPtr(new StateFreqParam(3));
		freqT_param->setName(freq_name[3]);
	    freqT_param->setTree(t);
	    freqT_param->setStartingValue(1.0);
	    freqT_param->setPrior(freq_param_prior);
	    if (state_freq_fixed)
		    freqT_param->fixParameter();
	    parameters.push_back(freqT_param);
	    freq_params.push_back(freqT_param);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the transition probability matrix given an edge length. Overrides the pure virtual function inherited from 
|	the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
void GTR::calcPMat(double * * pMat, double edgeLength) const
	{
	q_matrix.recalcPMat(pMat, edgeLength);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Needs work.
*/
double GTR::calcUniformizationLambda() const
	{
    PHYCAS_ASSERT(0);
    return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function 
|   inherited from the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
double GTR::calcLMat(double * * lMat) const
	{
    std::cerr << "Error in GTR::calcUMat: q_matrix does not yet have the required recalcLMat function" << std::endl;
    PHYCAS_ASSERT(0);
	//q_matrix.recalcLMat(lMat, edgeLength);
    return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function 
|   inherited from the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
double GTR::calcUMat(double * * uMat) const
	{
    std::cerr << "Error in GTR::calcLMat: q_matrix does not yet have the required recalcUMat function" << std::endl;
    PHYCAS_ASSERT(0);
	//q_matrix.recalcUMat(uMat, edgeLength);
    return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `rel_rates_fixed' to true. The fixParameter member function of all GTRRateParam objects is 
|	either called immediately (if the `rel_rate_params' vector is not empty) or is called in createParameters (when 
|	`rel_rate_params' is built).
*/
void GTR::fixRelRates()
	{
	rel_rates_fixed = true;
	if (!rel_rate_params.empty())
		{
#if 1
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
void GTR::freeRelRates()
	{
	rel_rates_fixed = false;
	if (!rel_rate_params.empty())
		{
#if 1
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
std::vector<double> GTR::getRelRates()
 	{
	return rel_rates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Copies the supplied vector `rates' into the data member vector `rel_rates'. Throws XLikelihood exception if any of
|	the values of `rates' is less than 0.0, or if the number of elements in `rates' is not equal to 6. Assumes that
|	`rel_rates' has 6 elements.
*/
void GTR::setRelRates(const std::vector<double> & rates)
 	{
	++time_stamp;
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
void GTR::setRelRateUnnorm(
  unsigned param_index,		/**< the 0-based index into the `rel_rates' vector of the element to modify */
  double value)				/**< the new value of `rel_rates'[`param_index'] */
	{
	++time_stamp;
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
MultivarProbDistShPtr GTR::getRelRatePrior()
 	{
	return rel_rate_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `rel_rate_param_prior'.
*/
ProbDistShPtr GTR::getRelRateParamPrior()
 	{
	return rel_rate_param_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `rel_rate_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void GTR::setRelRatePrior(MultivarProbDistShPtr d)
 	{
	rel_rate_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `rel_rate_param_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void GTR::setRelRateParamPrior(ProbDistShPtr d)
 	{
	rel_rate_param_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of Model base class function that sets all four nucleotide frequency parameters. The base class version is
|	called to do most of the work, and this function is responsible only for ensuring that the `q_matrix' data member
|	knows about the change in base frequencies.
*/
void GTR::setNucleotideFreqs(
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
void GTR::setAllFreqsEqual()
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
void GTR::setStateFreqUnnorm(
  unsigned param_index,		/**< the 0-based index into the `state_freq_unnorm' vector of the element to modify */
  double value)				/**< the new value of `state_freq_unnorm'[`param_index'] */
	{
	Model::setStateFreqUnnorm(param_index, value);
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets all of the state frequency parameters (the unnormalized values that determine the values
|	in the `state_freqs' vector when normalized) in the data member vector `state_freq_unnorm'. The base class version
|	is called to do most of the work, and this function is responsible only for ensuring that the `q_matrix' data member
|	knows about the change in base frequencies.
*/
void GTR::setStateFreqsUnnorm(
  const std::vector<double> & values)	/**< the new unnormalized state frequencies */
	{
	Model::setStateFreqsUnnorm(values);
		//std::copy(state_freqs.begin(), state_freqs.end(), std::ostream_iterator<double>(std::cerr, "|"));
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_param_prior'.
*/
ProbDistShPtr GTR::getStateFreqParamPrior()
 	{
	return freq_param_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_param_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void GTR::setStateFreqParamPrior(ProbDistShPtr d)
 	{
	freq_param_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_prior'.
*/
MultivarProbDistShPtr GTR::getStateFreqPrior()
 	{
	return freq_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_prior' data member to the supplied MultivariateProbabilityDistribution shared pointer `d'.
*/
void GTR::setStateFreqPrior(MultivarProbDistShPtr d)
 	{
	freq_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by 
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for 
|	compatability with MrBayes). The GTR model provide additional columns for the six relative rates, the base 
|	frequencies, the gamma shape parameter (if the number of rates is greater than 1) and the pinvar parameter (if
|	an invariable sites model is being used).
*/
std::string GTR::paramHeader() const	/**< is the suffix to tack onto the parameter names for this model (useful for partitioned models to show to which partition subset the parameter belongs) */
	{
	std::string s;
	string_vect_t::const_iterator it;
	for (it = relrate_name.begin(); it != relrate_name.end(); ++it)
		{
		s += boost::str(boost::format("\t%s") % (*it)); 
		}
	for (it = freq_name.begin(); it != freq_name.end(); ++it)
		{
		s += boost::str(boost::format("\t%s") % (*it)); 
		}
	s += Model::paramHeader();
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides the pure virtual base class version to generate a string of tab-separated values of model-specific 
|	parameters suitable for saving in a sampled parameter file (e.g. like the .p files saved by MrBayes). Assumes that
|	the `state_freqs' vector has length 4 and the `rel_rates' vector has length 6.
*/
std::string GTR::paramReport(
  unsigned ndecimals, 						/**< floating point precision to use */
  bool include_edgelen_hyperparams) const	/**< if true, include values of edge length hyperparameters */
	{
	PHYCAS_ASSERT(rel_rates.size() == 6);
	PHYCAS_ASSERT(state_freqs.size() == 4);
    std::string fmt = boost::str(boost::format("%%.%df\t") % ndecimals);

	std::string s = str(boost::format(fmt) % rel_rates[0]);
	s += str(boost::format(fmt) % rel_rates[1]);
	s += str(boost::format(fmt) % rel_rates[2]);
	s += str(boost::format(fmt) % rel_rates[3]);
	s += str(boost::format(fmt) % rel_rates[4]);
	s += str(boost::format(fmt) % rel_rates[5]);
	//std::string s = str(boost::format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f") % rel_rates[0] % rel_rates[1] % rel_rates[2] % rel_rates[3] % rel_rates[4] % rel_rates[5]);
	
	s += str(boost::format(fmt) % state_freqs[0]);
	s += str(boost::format(fmt) % state_freqs[1]);
	s += str(boost::format(fmt) % state_freqs[2]);
	s += str(boost::format(fmt) % state_freqs[3]);
	//s += str(boost::format("\t%.5f\t%.5f\t%.5f\t%.5f") % state_freqs[0] % state_freqs[1] % state_freqs[2] % state_freqs[3]);
	
	s += Model::paramReport(ndecimals, include_edgelen_hyperparams);
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
double GTR::calcTRatio()
	{
	double numerator = state_freqs[0]*state_freqs[2]*rel_rates[1] + state_freqs[3]*state_freqs[1]*rel_rates[4];
    double denominator = state_freqs[0]*state_freqs[1]*rel_rates[0] + state_freqs[0]*state_freqs[3]*rel_rates[2] + state_freqs[1]*state_freqs[2]*rel_rates[3] + state_freqs[2]*state_freqs[3]*rel_rates[5];
	PHYCAS_ASSERT(denominator != 0.0);
	double tratio = numerator/denominator;
	return tratio;
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
 */
void GTR::beagleGetStateFreqs(std::vector<double> & freqs)
{
	//std::cerr << "in GTR::beagleGetStateFreqs\n";
	q_matrix.beagleGetStateFreqs(freqs);
	//std::copy(freqs.begin(), freqs.end(), std::ostream_iterator<double>(std::cerr, "|"));
}

void GTR::beagleGetEigenValues(std::vector<double> & eigenValues)
{
	q_matrix.beagleGetEigenValues(eigenValues);
}

void GTR::beagleGetEigenVectors(std::vector<double> & eigenVectors)
{
	q_matrix.beagleGetEigenVectors(eigenVectors);
}

void GTR::beagleGetInverseEigenVectors(std::vector<double> & inverseEigenVectors)
{
	q_matrix.beagleGetInverseEigenVectors(inverseEigenVectors);
}

double GTR::beagleGetEdgelenScaler()
{
	return q_matrix.beagleGetEdgelenScaler();
}

