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
|	Constructor sets `num_states' data member to 4, the base frequencies `piA', `piC', `piG' and `piT' to 0.25, and `kappa'
|	to 1.0.
*/
HKY::HKY()
  : Model(4), kappa(1.0), kappa_fixed(false)
	{
	state_repr.reserve(4);
	state_repr.push_back("A");
	state_repr.push_back("C");
	state_repr.push_back("G");
	state_repr.push_back("T");
	}

HKY::~HKY()
	{
	//std::cerr << "HKY dying..." << std::endl;
	Model::Clear();
	//std::cerr << "  kappa_prior use_count() = " << kappa_prior.use_count() << std::endl;
	}
										
/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string indicating the name of this model, namely "HKY85", "HKY85+G", "HKY85+I" or "HKY85+G+I".
*/
std::string HKY::getModelName() const
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
void	HKY::createParameters(
  TreeShPtr t,								/**< is the tree (the nodes of which are needed for creating edge length parameters) */
  MCMCUpdaterVect & edgelens,				/**< is the vector of edge length parameters to fill */
  MCMCUpdaterVect & edgelen_hyperparams,	/**< is the edge length hyperparameter */
  MCMCUpdaterVect & parameters) const		/**< is the vector of model-specific parameters to fill */
	{
	Model::createParameters(t, edgelens, edgelen_hyperparams, parameters);

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

    PHYCAS_ASSERT(freq_param_prior || freq_prior);
    if (freq_param_prior)
        {
        // Only add frequency parameters if freqs will be updated separately
        // The other option is to update the frequencies jointly using the 
        // StateFreqMove Metropolis-Hastings move (in which case freq_param
        // will be set and freq_param_prior will be empty)

	    MCMCUpdaterShPtr freqA_param = MCMCUpdaterShPtr(new StateFreqParam(0));
	    freqA_param->setName("base freq. A");
	    freqA_param->setTree(t);
	    freqA_param->setStartingValue(1.0);
	    freqA_param->setPrior(freq_param_prior);
	    if (state_freq_fixed)
		    freqA_param->fixParameter();
	    parameters.push_back(freqA_param);
	    freq_params.push_back(freqA_param);

	    MCMCUpdaterShPtr freqC_param = MCMCUpdaterShPtr(new StateFreqParam(1));
	    freqC_param->setName("base freq. C");
	    freqC_param->setTree(t);
	    freqC_param->setStartingValue(1.0);
	    freqC_param->setPrior(freq_param_prior);
	    if (state_freq_fixed)
		    freqC_param->fixParameter();
	    parameters.push_back(freqC_param);
	    freq_params.push_back(freqC_param);

	    MCMCUpdaterShPtr freqG_param = MCMCUpdaterShPtr(new StateFreqParam(2));
	    freqG_param->setName("base freq. G");
	    freqG_param->setTree(t);
	    freqG_param->setStartingValue(1.0);
	    freqG_param->setPrior(freq_param_prior);
	    if (state_freq_fixed)
		    freqG_param->fixParameter();
	    parameters.push_back(freqG_param);
	    freq_params.push_back(freqG_param);

	    MCMCUpdaterShPtr freqT_param = MCMCUpdaterShPtr(new StateFreqParam(3));
	    freqT_param->setName("base freq. T");
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
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by 
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for 
|	compatability with MrBayes). The HKY model provide additional columns for kappa, the base frequencies, the gamma 
|	shape parameter (if the number of rates is greater than 1) and the pinvar parameter (if an invariable sites model 
|	is being used)
*/
std::string HKY::paramHeader() const
	{
	std::string s = std::string("\tkappa\tfreqA\tfreqC\tfreqG\tfreqT");
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
std::string HKY::paramReport(
  unsigned ndecimals) const /**< floating point precision to use */
	{
    std::string fmt = boost::str(boost::format("%%.%df\t") % ndecimals);
	std::string s = str(boost::format(fmt) % kappa);
	s += str(boost::format(fmt) % state_freqs[0]);
	s += str(boost::format(fmt) % state_freqs[1]);
	s += str(boost::format(fmt) % state_freqs[2]);
	s += str(boost::format(fmt) % state_freqs[3]);
	//std::string s = str(boost::format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f") % kappa % state_freqs[0] % state_freqs[1] % state_freqs[2] % state_freqs[3]);
	if (is_flex_model)
		{
		s += str(boost::format("%d\t") % num_gamma_rates);
		//s += str(boost::format("\t%d") % num_gamma_rates);

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
		s += str(boost::format(fmt) % gamma_shape);
		//s += str(boost::format("\t%.5f") % gamma_shape);
		}
	if (is_pinvar_model)
        {
		s += str(boost::format(fmt) % pinvar);
		//s += str(boost::format("\t%.5f") % pinvar);
        }
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `kappa_fixed' to true. The fixParameter member function of the KappaParam object is either 
|	called immediately (if `kappa_param' is a valid pointer) or is called in createParameters (when `kappa_param' is 
|	first assigned).
*/
void HKY::fixKappa()
	{
	kappa_fixed = true;
	if (kappa_param)
		kappa_param->fixParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `kappa_fixed' to false. The freeParameter member function of the KappaParam object is called 
|	immediately if `kappa_param' is a valid pointer.
*/
void HKY::freeKappa()
	{
	kappa_fixed = false;
	if (kappa_param)
		kappa_param->freeParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `kappa'.
*/
double HKY::getKappa()
 	{
	return kappa;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa' to supplied value `k'. Throws XLikelihood exception if `k' is less than or equal to 0.0.
*/
void HKY::setKappa(double k)
 	{
	if (k <= 0.0)
		throw XLikelihood();
	kappa = k;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `kappa_prior'.
*/
ProbDistShPtr HKY::getKappaPrior()
 	{
	return kappa_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void HKY::setKappaPrior(ProbDistShPtr d)
 	{
	kappa_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_param_prior'.
*/
ProbDistShPtr HKY::getStateFreqParamPrior()
 	{
	return freq_param_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_param_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void HKY::setStateFreqParamPrior(ProbDistShPtr d)
 	{
	freq_param_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_prior'.
*/
MultivarProbDistShPtr HKY::getStateFreqPrior()
 	{
	return freq_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_prior' data member to the supplied MultivariateProbabilityDistribution shared pointer `d'.
*/
void HKY::setStateFreqPrior(MultivarProbDistShPtr d)
 	{
	freq_prior = d;
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
void HKY::setKappaFromTRatio(double tratio)
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

double HKY::calcTRatio()
	{
	double numerator   = kappa*((state_freqs[0]*state_freqs[2]) + (state_freqs[1]*state_freqs[3]));
	double denominator = (state_freqs[0] + state_freqs[2])*(state_freqs[1] + state_freqs[3]);
	PHYCAS_ASSERT(denominator != 0.0);
	double tratio = numerator/denominator;
	return tratio;
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
