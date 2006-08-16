#include <cmath>
#include <iostream>
#include "pyphy/src/likelihood_models.hpp"
#include "pyphy/src/cipres/AllocateMatrix.hpp"
#if defined(PYTHON_ONLY)
#	include <boost/python/numeric.hpp>
#	include "pyphy/src/thirdparty/num_util.h"
#endif
using std::cout;
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The base class version of this function should be called by all derived classes because it is where the edge 
|	length parameters, edge length hyperparameter and rate heterogeneity (gamma shape and pinvar) parameters are 
|	created.
*/
void Model::createParameters(
  TreeShPtr t,									/**< is the tree (the nodes of which are needed for creating edge length parameters) */
  MCMCUpdaterVect & edgelens_vect_ref,			/**< is the vector of edge length parameters to fill */
  MCMCUpdaterShPtr & edgelen_hyperparam_ref,	/**< is the edge length hyperparameter */
  MCMCUpdaterVect & parameters_vect_ref,		/**< is the vector of model-specific parameters to fill */
  bool separate_edgelens) const					/**< specifies (if true) that each edge should have its own parameter or (if false) that one edge length master parameter should be created */
	{
	assert(t);
	assert(edgelens_vect_ref.empty());
	assert(!edgelen_hyperparam_ref);
	assert(parameters_vect_ref.empty());

	// Add the edge length parameter(s)
	if (separate_edgelens)
		{
		// Add an edge length parameter for every node in the tree that has an edge
		TreeNode * nd = t->GetFirstPreorder();	// this is the root node, which has no edge
		nd = nd->GetNextPreorder();
		for (; nd != NULL; nd = nd->GetNextPreorder())
			{
			MCMCUpdaterShPtr p = MCMCUpdaterShPtr(new EdgeLenParam(nd));
			std::string nm = str(boost::format("edge length for node %d") % nd->GetNodeNumber());
			p->setName(nm);
			p->setTree(t);
			p->setPrior(edgeLenPrior);
			if (edge_lengths_fixed)
				p->fixParameter();
			edgelens_vect_ref.push_back(p);
			}
		}
	else
		{
		// Add one edge length parameter to manage the prior for all edge lengths in the tree;
		// This edge length master parameter does not actually update edge lengths; when 
		// separate_edgelen_params is false, this means that a Metropolis proposal such as
		// the LargetSimonMove is responsible for updating edge lengths.
		MCMCUpdaterShPtr p = MCMCUpdaterShPtr(new EdgeLenMasterParam());
		std::string nm = str(boost::format("edge length master parameter"));
		p->setName(nm);
		p->setTree(t);
		p->setPrior(edgeLenPrior);
		if (edge_lengths_fixed)
			p->fixParameter();
		edgelens_vect_ref.push_back(p);
		}

	// Save a vector of shared pointers to the edge length parameters so that we can modify their
	// fixed/free status if we need to
	edgelen_params.resize(edgelens_vect_ref.size());
	std::copy(edgelens_vect_ref.begin(), edgelens_vect_ref.end(), edgelen_params.begin());

	// Add the edge length hyperparameter if requested
	if (edgeLenHyperPrior)
		{
		MCMCUpdaterShPtr p = MCMCUpdaterShPtr(new HyperPriorParam());
		p->setName(std::string("edge length hyperprior"));
		p->setTree(t);
		p->setPrior(edgeLenHyperPrior);
		if (edgelen_hyperprior_fixed)
			p->fixParameter();
		edgelen_hyperparam_ref = p;

		// Retain a copoy of the shared pointer so that we can later modify the fixed/free status
		// of this parameter
		edgelen_hyper_param = p;
		}

	// Create any model-specific parameters and add to the parameters vector
	if (is_flex_model)
		{
		gamma_rates_unnorm.resize(num_gamma_rates, 0.0);
		assert(flex_rate_params.empty());
		assert(flex_prob_params.empty());
		for (unsigned i = 0; i < num_gamma_rates; ++i)
			{
			// start with rates drawn from Uniform(0.0, flex_upper_rate_bound)
			//@POL to do this right, need to draw from prior, but this will be close if number of spacers is small
			double u = flex_prob_param_prior->GetLot()->Uniform(FILE_AND_LINE);
			gamma_rates_unnorm[i] = flex_upper_rate_bound*u;
			//old way: gamma_rates_unnorm[i] = flex_upper_rate_bound*(double)(i + 1)/(double)(num_gamma_rates + 1);

			// start with probabilities all equal
			assert(flex_prob_param_prior);
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
		assert(num_gamma_rates > 1);
		assert(!gamma_shape_param);
		gamma_shape_param = MCMCUpdaterShPtr(new DiscreteGammaShapeParam(invert_shape));
		gamma_shape_param->setName("Discrete gamma shape"); //@POL shouldn't this be done in the constructor?
		gamma_shape_param->setTree(t);
		gamma_shape_param->setPrior(gamma_shape_prior);
		if (gamma_shape_fixed)
			gamma_shape_param->fixParameter();
		parameters_vect_ref.push_back(gamma_shape_param);
		}

	if (is_pinvar_model)
		{
		assert(!pinvar_param);
		pinvar_param = MCMCUpdaterShPtr(new PinvarParam());
		pinvar_param->setName("Proportion of invariable sites");
		pinvar_param->setTree(t);
		pinvar_param->setPrior(pinvar_prior);
		if (pinvar_fixed)
			pinvar_param->fixParameter();
		parameters_vect_ref.push_back(pinvar_param);
		}
	}
	
#if defined(PYTHON_ONLY)
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
