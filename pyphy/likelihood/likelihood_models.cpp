#include <cmath>
#include <iostream>
#include "likelihood_models.hpp"
#include "CipresCommlib/AllocateMatrix.hpp"
#include <boost/python/numeric.hpp>
#include "pyphy/include/num_util.h"
using std::cout;
using namespace phycas;

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Normalizes the rates in `unnorm_pat_spec_rates' so that the mean rate across all sites is 1.0, and copies the 
|	normalized values into the supplied vector.
*/
void Model::normalizePatternSpecificRates(
  std::vector<double> & rates,  /**< is the vector that should receive the normalized rates */
  CountVectorType & counts)		/**< is the vector of pattern counts needed to normalize the rates to have mean 1.0 */
	{
	assert(unnorm_pat_spec_rates.size() == rates.size());
	assert(unnorm_pat_spec_rates.size() == counts.size());
	//@POL when there is time, convert this to use algorithms instead of loops
	CountVectorType::iterator cit = counts.begin();
	std::vector<double>::iterator rit = unnorm_pat_spec_rates.begin();
	double numerator = 0.0;
	double denominator = 0.0;
	for (; cit != counts.end(); ++cit, ++rit)
		{
		double count = (double)(*cit);
		double rate = *rit;
		numerator += (count*rate);
		denominator += count;
		}
	assert(numerator > 0.0);
	assert(denominator > 0.0);
	double norm_factor = numerator/denominator;
	//std::cerr << "norm_factor = " << norm_factor << std::endl;
	std::vector<double>::iterator it = rates.begin();
	rit = unnorm_pat_spec_rates.begin();
	for (; it != rates.end(); ++it, ++rit)
		{
		*it = (*rit)/norm_factor;
		}
	//std::transform(unnorm_pat_spec_rates.begin(), unnorm_pat_spec_rates.end(), rates.begin(), boost::lambda::_1/norm_factor);

	// former (wrong) code is below
	//double sum = std::accumulate(unnorm_pat_spec_rates.begin(), unnorm_pat_spec_rates.end(), 0.0);
	//assert(sum > 0.0);
	//std::transform(unnorm_pat_spec_rates.begin(), unnorm_pat_spec_rates.end(), rates.begin(), boost::lambda::_1/sum);
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
boost::python::numeric::array Model::getPMatrix(double edgeLength) const
	{
	//@POL this function, along with calcPMat and calcPMatrices, should be provided by a policy class
	// (which would be QMatrix for GTR model). This function is full of ad hoc nonsense at the moment!
	double * * pMat = NewTwoDArray<double>(num_states, num_states);
	calcPMat(pMat, edgeLength);

	// copy to a vector so that we can delete pMat
	unsigned flat_length = num_states*num_states;
	std::vector<double> p(flat_length, 0.0);
	double * pMat_begin = &pMat[0][0];
	double * pMat_end   = pMat_begin + flat_length;
	std::copy(pMat_begin, pMat_end, p.begin());

	DeleteTwoDArray<double>(pMat);

	// create vector of dimensions
	std::vector<int> dim_vect;
	dim_vect.push_back((int)num_states);
	dim_vect.push_back((int)num_states);

	return num_util::makeNum(&p[0], dim_vect);	//PELIGROSO
	}

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
	assert(state_freqs.size() == 4);
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
