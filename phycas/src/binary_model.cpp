/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2011 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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
|	Constructor sets `num_states' data member to 2 and 'root_state' data member to 1.
*/
Binary::Binary()
  : Model(2),phi(1.0),rho(1.0)
	{
    state_freq_fixed = false;
        
    state_repr.reserve(2);
	state_repr.push_back("0.5");
	state_repr.push_back("0.5");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `phi_fixed' data member to false.
*/
void Binary::freeScalingFactor() 
    {
    phi_fixed = false;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `phi_fixed' data member to true.
*/
void Binary::fixScalingFactor() 
    {
    phi_fixed = true;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `phi' data member.
*/
double Binary::getScalingFactor()
    {
    return phi;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `phi' data member to the supplied scaling factor `sf'. Assumes `sf' is greater than or equal to zero.
*/
void Binary::setScalingFactor(
  double sf)    /**< is the new value of phi */
    {
    PHYCAS_ASSERT(sf >= 0.0);
    phi = sf;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `rho_fixed' data member to false.
*/
void Binary::freeRevForRateRatio() 
    {
    rho_fixed = false;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `rho_fixed' data member to true.
*/
void Binary::fixRevForRateRatio() 
    {
    rho_fixed = true;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `rho' data member.
*/
double Binary::getRevForRateRatio()
    {
    return rho;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `rho' data member to the supplied reverse/forward rate ratio `r'. Assumes `r' is greater than or equal to zero.
*/
void Binary::setRevForRateRatio(
  double r)    /**< is the new value of phi */
    {
    PHYCAS_ASSERT(r >= 0.0);
    rho = r;
}   

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string indicating the name of this model, namely "Binary", "Binary+I", "Binary+G", or "Binary+I+G".
*/
std::string Binary::getModelName() const
	{
	std::string s = "Binary";
	if (is_pinvar_model)
		s += "+I";
	if (num_gamma_rates > 1)
		s += "+G";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by 
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for 
|	compatability with MrBayes). The Binary model provides additional column labels for the gamma shape parameter (if
|   the number of rates is greater than 1) and the pinvar parameter (if an invariable sites model is being used).
*/
std::string Binary::paramHeader() const	/**< is the suffix to tack onto the parameter names for this model (useful for partitioned models to show to which partition subset the parameter belongs) */
	{
	return Model::paramHeader();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides the pure virtual base class version to generate a string of tab-separated values of model-specific 
|	parameters suitable for saving in a sampled parameter file (e.g. like the .p files saved by MrBayes).
*/
std::string Binary::paramReport(
  unsigned ndecimals,						/**< floating point precision to use */
  bool include_edgelen_hyperparams) const	/**< if true, include values of edge length hyperparameters */
    {
    std::string fmt = boost::str(boost::format("%%.%df\t") % ndecimals);
	std::string s = "";
    s += str(boost::format(fmt) % phi);
    s += str(boost::format(fmt) % rho);
    s += str(boost::format(fmt) % state_freqs[0]);
	s += str(boost::format(fmt) % state_freqs[1]);
	s += Model::paramReport(ndecimals, include_edgelen_hyperparams);
	return s;
    }

/*----------------------------------------------------------------------------------------------------------------------
 |	Computes the transition probability matrix given an edge length. Overrides the pure virtual function inherited from 
 |	the base class Model. For the Binary model, the transition probabilities are:
 |>
 |	P00 = pi[0] + (1 - pi[0])*exp{-beta*t}
 |	P01 = pi[1] - pi[1]*exp{-beta*t}
 |	P10 = pi[0] - pi[0]*exp{-beta*t}
 |	P11 = pi[1] + (1 - pi[1])*exp{-beta*t}
 |>
 |	The evolutionary distance `edgeLength' = beta*t, so as a function of edge length, the transition probabilities
 |	are:
 |>
 |	P00 = pi[0] + (1 - pi[0])*exp{-edgeLength}
 |	P01 = pi[1] - pi[1]*exp{-edgeLength}
 |	P10 = pi[0] - pi[0]*exp{-edgeLength}
 |	P11 = pi[1] + (1 - pi[1])*exp{-edgeLength}
 |>
 */
void Binary::calcPMat(double * * pMat, double edgeLength) const
{
	// The next two lines fix the "Jockusch" bug (see bug 8 in the BUGS file for details)
    if (edgeLength < 1.e-8) 
        edgeLength = 1.e-8; //TreeNode::edgeLenEpsilon;
	const double exp_term = exp(-edgeLength);
	double pi0 = state_freqs[0];
	double pi1 = state_freqs[1];
	pMat[0][0] = pi0 + (1 - pi0)*exp_term;
	pMat[0][1] = pi1 - pi1*exp_term;
	pMat[1][0] = pi0 - pi0*exp_term;
	pMat[1][1] = pi1 + (1 - pi1)*exp_term;
}

/*----------------------------------------------------------------------------------------------------------------------
|   Needs work.
*/
double Binary::calcUniformizationLambda() const
	{
    return 1.1;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Needs work.
*/
double Binary::calcLMat(double * * lMat) const
	{
    // The body of this function is totally bogus right now and this model should not be used
    // with uniformization until something better is put in place
    const double lambda = 1.1;
	lMat[0][0] = 0.0;
	lMat[0][1] = 0.0;
	lMat[1][0] = 0.0;
	lMat[1][1] = 0.0;
    return lambda;
	}

double Binary::calcUMat(double * * uMat) const
	{
    // The body of this function is totally bogus right now and this model should not be used
    // with uniformization until something better is put in place
    const double lambda = 1.1;
    uMat[0][0] = 0.0;
    uMat[0][1] = 0.0;
    uMat[1][0] = 0.0;
    uMat[1][1] = 0.0;
    return lambda;
	}

