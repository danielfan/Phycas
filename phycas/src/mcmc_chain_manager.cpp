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
#include <cstdlib>

#include <fstream>//temp 

#include <iterator>		// for std::distance 
#include <algorithm>	// for std::lower_bound  

//#include "phycas/force_include.h"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/mapping_move.hpp"
#include "phycas/src/dirichlet_move.hpp"

extern "C"
{
//#include "phycas/src/thirdparty/praxis/machine.h"
//double praxis(double (*_fun)(double *, int), double * _x, int _n);
//double praxisWrapper(double * x, int n);
#include "phycas/src/thirdparty/praxis/dls_core.h"
#include "phycas/src/thirdparty/praxis/dls_brent.h"
}

//@POL should use ChainManagerShPtr here
phycas::MCMCChainManager * praxisMCMCChainManager = 0;

double praxisWrapper(double * x, void *data)
    {
    assert(praxisMCMCChainManager);
    return -1.0*praxisMCMCChainManager->praxisCalcLogPosterior(x);
    }

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of `samc_ref_tree', which should be a shared pointer to the best tree known.
*/
void MCMCChainManager::setRefTree(
  TreeShPtr t) 	/**< is the reference tree */
	{
	ref_tree = t;
	}
	
void MCMCChainManager::debugUpdaterReport(std::string s)
	{
	std::cerr << "\n>>>>> MCMCChainManager updater report (" << s << ") <<<<<" << std::endl;
	std::cerr << "all_updaters:" << std::endl;
	for (MCMCUpdaterIter iter = all_updaters.begin(); iter != all_updaters.end(); ++iter)
		{
		std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
		}
	std::cerr << "moves:" << std::endl;
	for (MCMCUpdaterIter iter = moves.begin(); iter != moves.end(); ++iter)
		{
		std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
		}
	std::cerr << "model_params:" << std::endl;
	for (MCMCUpdaterIter iter = model_params.begin(); iter != model_params.end(); ++iter)
		{
		std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
		}
	std::cerr << "edge_len_params:" << std::endl;
	for (MCMCUpdaterIter iter = edge_len_params.begin(); iter != edge_len_params.end(); ++iter)
		{
		std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
		}
	std::cerr << "edge_len_hyperparams:" << std::endl;
	for (MCMCUpdaterIter iter = edge_len_hyperparams.begin(); iter != edge_len_hyperparams.end(); ++iter)
		{
		std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
		}
	std::cerr << "\n" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor for 
*/
MCMCChainManager::~MCMCChainManager() 
	{
	//std::cerr << "\n>>>>> MCMCChainManager dying..." << std::endl;
	releaseUpdaters();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears vectors of shared pointers to MCMCUpdater objects. Called from destructor (or can be called at will). 
*/
void MCMCChainManager::releaseUpdaters()
    {
	//std::cerr << "\n>>>>> MCMCChainManager::releaseUpdaters()..." << std::endl;    
	//std::cerr << "\nBefore all_updaters cleared..." << std::endl;
	for (MCMCUpdaterIter iter = all_updaters.begin(); iter != all_updaters.end(); ++iter)
		{
		//std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
		(*iter)->releaseSharedPointers();
		}
	all_updaters.clear();

	//std::cerr << "\nBefore moves cleared..." << std::endl;
	//for (MCMCUpdaterIter iter = moves.begin(); iter != moves.end(); ++iter)
	//	{
	//	std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
	//	}
	moves.clear();

	//std::cerr << "\nBefore model_params cleared..." << std::endl;
	//for (MCMCUpdaterIter iter = model_params.begin(); iter != model_params.end(); ++iter)
	//	{
	//	std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
	//	}
	model_params.clear();

	//std::cerr << "\nBefore edge_len_params cleared..." << std::endl;
	//for (MCMCUpdaterIter iter = edge_len_params.begin(); iter != edge_len_params.end(); ++iter)
	//	{
	//	std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
	//	}
	edge_len_params.clear();

	//std::cerr << "\nBefore edge_len_hyperparams cleared..." << std::endl;
	//for (MCMCUpdaterIter iter = edge_len_hyperparams.begin(); iter != edge_len_hyperparams.end(); ++iter)
	//	{
	//	std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
	//	}
	edge_len_hyperparams.clear();
    }
    
/*----------------------------------------------------------------------------------------------------------------------
|	Locates the mode of the posterior distribution. The supplied `starting_point' is informational only; the actual 
|   starting point is represented by the current state of the model. It is thus the duty of the caller to ensure that the
|   model is setup to match the value of `starting_point'.
*/
void MCMCChainManager::praxisLocatePosteriorMode()
	{	
	// Might need to add more checks here
	if (dirty)
		{
		throw XLikelihood("cannot call refreshLastLnPrior() for chain manager before calling finalize()");
		}

    praxis_param_values.clear();
    praxis_stewards.clear();
    
    double_vect_t tmp_param_values;
    //double TL = 0.0;
    
    std::cerr << "\n@@@@@@@@@@ Building vector of parameter names in MCMCChainManager::praxisLocatePosteriorMode() @@@@@@@@@@" << std::endl;
    std::cerr << "@@@@@@@@@@ ignored steward names:" << std::endl;
    
	// Visit all updaters and let those who are prior stewards add their name to praxis_stewards.
	for (MCMCUpdaterConstIter it = all_updaters.begin(); it != all_updaters.end(); ++it)
		{
		const MCMCUpdaterShPtr s = *it;
		if (s->isPriorSteward() && !s->isFixed())
			{
            const std::string & nm = s->getName();
            if (nm.find("state_freqs") != std::string::npos) 
                {
                const StateFreqMove & state_freq_move = dynamic_cast<const StateFreqMove &>(*s);
                std::vector<double> state_freqs(4, 0.0);
                state_freq_move.getCurrValuesFromModel(state_freqs);
                tmp_param_values.push_back(state_freqs[0]);
                tmp_param_values.push_back(state_freqs[1]);
                tmp_param_values.push_back(state_freqs[2]);
                tmp_param_values.push_back(state_freqs[3]);
                double log_freqA = log(state_freqs[0]);
                state_freqs[0] = log(state_freqs[1]) - log_freqA;   // freqC
                state_freqs[1] = log(state_freqs[2]) - log_freqA;   // freqG
                state_freqs[2] = log(state_freqs[3]) - log_freqA;   // freqT
                state_freqs.resize(3);
                std::copy(state_freqs.begin(), state_freqs.end(), std::back_inserter(praxis_param_values));
                praxis_stewards.push_back(s);
                }
            else if (nm.find("relrates") != std::string::npos) 
                {
                const RelRatesMove & relrates_move = dynamic_cast<const RelRatesMove &>(*s);
                std::vector<double> relative_rates(6, 0.0);
                relrates_move.getCurrValuesFromModel(relative_rates);
                tmp_param_values.push_back(relative_rates[0]);
                tmp_param_values.push_back(relative_rates[1]);
                tmp_param_values.push_back(relative_rates[2]);
                tmp_param_values.push_back(relative_rates[3]);
                tmp_param_values.push_back(relative_rates[4]);
                tmp_param_values.push_back(relative_rates[5]);
                double log_rAC = log(relative_rates[0]);
                relative_rates[0] = log(relative_rates[1]) - log_rAC;   // rAG
                relative_rates[1] = log(relative_rates[2]) - log_rAC;   // rAT
                relative_rates[2] = log(relative_rates[3]) - log_rAC;   // rCG
                relative_rates[3] = log(relative_rates[4]) - log_rAC;   // rCT
                relative_rates[4] = log(relative_rates[5]) - log_rAC;   // rGT
                relative_rates.resize(5);
                std::copy(relative_rates.begin(), relative_rates.end(), std::back_inserter(praxis_param_values));
                praxis_stewards.push_back(s);
                }
            else if (nm.find("gamma_shape") != std::string::npos) 
                {
                const DiscreteGammaShapeParam & gamma_shape_parameter = dynamic_cast<const DiscreteGammaShapeParam &>(*s);
                double shape = gamma_shape_parameter.getCurrValueFromModel();
                double log_shape = log(shape);
                tmp_param_values.push_back(shape);
                praxis_param_values.push_back(log_shape);
                praxis_stewards.push_back(s);
                }
            else if (nm.find("master_edgelen") != std::string::npos) 
                {
                EdgeLenMasterParam & edgelen_master_parameter = dynamic_cast<EdgeLenMasterParam &>(*s);
                TreeShPtr tree = edgelen_master_parameter.getTree();
                TreeNode * nd = tree->GetFirstPreorder();
                for (nd = nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
                    {
                    double edge_len = (nd->GetEdgeLen());
                    //TL += edge_len;
                    tmp_param_values.push_back(edge_len);
                    double log_edge_len = log(edge_len);
                    praxis_param_values.push_back(log_edge_len);
                    }
                std::cerr << boost::str(boost::format("tree = %s") % tree->MakeNewick()) << std::endl;
                praxis_stewards.push_back(s);
                }
            else 
                {
                std::cerr << nm << " " << std::endl;
                }
			}
		}
    std::cerr << "\n@@@@@@@@@@ praxis_stewards names:" << std::endl;
    for (MCMCUpdaterVect::iterator s = praxis_stewards.begin(); s != praxis_stewards.end(); ++s)
        {
        std::cerr << (*s)->getName() << " ";
        }
    std::cerr << "\n@@@@@@@@@@ praxis_param_values:" << std::endl;
    std::copy(praxis_param_values.begin(), praxis_param_values.end(), std::ostream_iterator<double>(std::cerr, " "));
    
    std::cerr << "\n@@@@@@@@@@ tmp_param_values:" << std::endl;
    std::copy(tmp_param_values.begin(), tmp_param_values.end(), std::ostream_iterator<double>(std::cerr, " "));
    
    int n = (int)praxis_param_values.size();
    double fx = praxisCalcLogPosterior(&praxis_param_values[0]);
    std::cerr << boost::str(boost::format("\nPosterior at starting point = %g\n") % fx) << std::endl;    

    // Note: praxisMCMCChainManager is a global (see top of this file)
    praxisMCMCChainManager = this;
    //fx = praxis(praxisWrapper, &praxis_param_values[0], n);
    
    //|		============
    //|		praxis usage
    //|		============
    //|
    //|			 1. Allocate and initialize a praxis object using newPraxisData:
    //|
    //|					PraxisData * praxd = newPraxisData(x, n, single_precision);
    //|					
    //|				  		x = array for parameter vector (must be large enough to store n doubles)
    //|						n = number of parameters in function being optimized
    //|						single_precision = true if function evaluation is only done at single-precision accuracy,
    //|						                   false otherwise (i.e., function is evaluated in double precision).
    //|
    //|			 2.	Initialize x to contain the starting parameter values (in most cases, the better the guess, the faster
    //|			    the convergence will be)
    //|
    //|			 3. Invoke the praxis routine:
    //|		 		
    //|			 		fx_opt = praxis(tol, h, n, f, data, praxd);
    //|					
    //|						tol   = tolerance limit used to test for convergence (if passed as 0, a default of 1e-5
    //|						        will be used)
    //|						h     = maximum step size used to test for convergence (if passed as 0, a default of 1.0
    //|							    will be used)
    //|						n     = the number of parameters in the function being minimized
    //|						f	  = the function being minimized (see above)
    //|						data  = an option pointer to arbitrary data that the function may need to compute its value;
    //|						        if the function needs no additional data, pass NULL
    //|						praxd = a praxis object previously allocated by newPraxisData
    //|						
    //|					On return, fx_opt will be the value of the minimized function, and praxd->h will contain the
    //|					corresponding optimized parameter values.
    //|
    //|			 4. Destroy the praxis object (to avoid leaking memory);
    //|
    //|					deletePraxisData(praxd);
    //|
    //|	@note	'gAborted' is a variable that can be set to request cancellation without terminating the program (e.g., by
    //|			a signal handler invoked when ctrl-C is pressed).  By default it is just #define'd to 0 in the header; if
    //|			you want to provide this capability, you will need to provide a mechanism for setting gAborted externally.
    //PUBLIC double praxis(
    //  double		tol,		/* tolerance used for convergence criterion */
    //  double		h,			/* maximum step size (e.g., 1.0) */
    //  int			n,			/* number of parameters to function */
    //  MinimizeFxn	f,			/* the function to be minimized, declared as "double fxn(double *x, void *data)" */
    //  void *		data,		/* pointer to data to be passed to 'f' */
    //  PraxisData *	praxd)		/* an object already allocated by newPraxisData */

    if (true)
        {
        PraxisData * praxd = newPraxisData(&praxis_param_values[0], n, false);
        double tol = 1e-5;
        double h = 1.0;
        double fx_opt = praxis(tol, h, n, praxisWrapper, NULL, praxd);
        std::cerr << boost::str(boost::format("Minimum = %g\n") % fx_opt) << std::endl;    
        deletePraxisData(praxd);
        }
    }
    
/*----------------------------------------------------------------------------------------------------------------------
|	Refreshes `last_ln_like' and `last_ln_prior' by calling refreshLastLnLike() and refreshLastLnPrior(), respectively,
|   then returning `last_ln_like' + `last_ln_prior'.
*/
double MCMCChainManager::praxisCalcLogPosterior(double * x)
	{
    double * cursor = &x[0];
    //std::cerr << "\n--->\n";
    for (MCMCUpdaterVect::iterator s = praxis_stewards.begin(); s != praxis_stewards.end(); ++s)
        {
        const std::string & nm = (*s)->getName();
        MCMCUpdaterShPtr p = (*s);
        MCMCUpdater & u = (*p);
        if (nm.find("state_freqs") != std::string::npos) 
            {
            StateFreqMove & state_freq_move = dynamic_cast<StateFreqMove &>(u);
            std::vector<double> state_freqs(4, 1.0);
            std::vector<double>::iterator it = state_freqs.begin();
            std::copy(cursor, cursor + 3, ++it);
            double freqA = 1.0/(1.0 + exp(state_freqs[1]) + exp(state_freqs[2]) + exp(state_freqs[3]));
            state_freqs[0] = freqA;
            //std::cerr << boost::str(boost::format("freq[0] = %g") % state_freqs[0]) << std::endl;
            for (unsigned i = 1; i < 4; ++i)
                {
                state_freqs[i] = freqA*exp(state_freqs[i]);
                //std::cerr << boost::str(boost::format("freq[%d] = %g") % i % state_freqs[i]) << std::endl;
                }
            state_freq_move.sendCurrValuesToModel(state_freqs);
            cursor += 3;
            }
        else if (nm.find("relrates") != std::string::npos) 
            {
            RelRatesMove & relrates_move = dynamic_cast<RelRatesMove &>(u);
            std::vector<double> relative_rates(6, 1.0);
            std::vector<double>::iterator it = relative_rates.begin();
            std::copy(cursor, cursor + 5, ++it);
            double rAC = 1.0/(1.0 + exp(relative_rates[1]) + exp(relative_rates[2]) + exp(relative_rates[3]) + exp(relative_rates[4]) + exp(relative_rates[5]));
            relative_rates[0] = rAC;
            //std::cerr << boost::str(boost::format("relrate[0] = %g") % relative_rates[0]) << std::endl;
            for (unsigned i = 1; i < 6; ++i)
                {
                relative_rates[i] = rAC*exp(relative_rates[i]);
                //std::cerr << boost::str(boost::format("relrate[%d] = %g") % i % relative_rates[i]) << std::endl;
                }
            relrates_move.sendCurrValuesToModel(relative_rates);
            cursor += 5;
            }
        else if (nm.find("gamma_shape") != std::string::npos) 
            {
            DiscreteGammaShapeParam & gamma_shape_parameter = dynamic_cast<DiscreteGammaShapeParam &>(u);
            double shape = exp(*cursor);
            //std::cerr << boost::str(boost::format("shape = %g") % shape) << std::endl;
            gamma_shape_parameter.sendCurrValueToModel(shape);
            gamma_shape_parameter.recalcRelativeRates();
            cursor++;
            }
        else if (nm.find("master_edgelen") != std::string::npos) 
            {
            EdgeLenMasterParam & edgelen_master_parameter = dynamic_cast<EdgeLenMasterParam &>(u);
            TreeShPtr tree = edgelen_master_parameter.getTree();
            TreeNode * nd = tree->GetFirstPreorder();
            //double TL = 0.0;
            for (nd = nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
                {
                double edge_len = exp(*cursor);
                //TL += edge_len;
                nd->SetEdgeLen(edge_len);
                cursor++;
                }
            //std::cerr << boost::str(boost::format("\npraxisCalcLogPosterior tree = %s\n") % tree->MakeNewick()) << std::endl;
            //std::cerr << boost::str(boost::format("TL (after) = %g") % TL) << std::endl;
            }
        else
            {
            std::cerr << "\n@@@@@@@@@@ warning: unhandled prior steward (" << nm << ") in MCMCChainManager::praxisCalcLogPosterior" << std::endl;
            }
        }
    refreshLastLnLike();
    refreshLastLnPrior();
    //std::cerr << "lnP = " << (last_ln_like + last_ln_prior) << std::endl;
    return last_ln_like + last_ln_prior;
    }
    
/*----------------------------------------------------------------------------------------------------------------------
|	Refreshes `last_ln_like' by calling recalcLike() for the first updater in the `all_updaters' vector.
*/
void MCMCChainManager::refreshLastLnLike()
	{

	// Might need to add more checks here
	if (dirty)
		{
		throw XLikelihood("cannot call refreshLastLnLike() for chain manager before calling finalize()");
		}

	// Only need one updater to get the log-likelihood because all updaters must be able to do this calculation

    //@POL should only go to this trouble if univents are being used
    // Try to find a MappingMove updater to use for this. If one is present, its update method must
    // be called before the likelihood can be calculated in order to generate the latent variables used in
    // uniformized mapping MCMC
    bool found = false; 
	for (MCMCUpdaterIter it = moves_begin; it != moves_end; ++it)
        {
        MappingMove * m = dynamic_cast<MappingMove *>(it->get());
        if (m)
            {
            m->update();
            last_ln_like = m->recalcLike();
            found = true;
            break;
            }
        }
    if (!found)
        {
        MCMCUpdaterShPtr u = *params_begin;
        //std::cerr << "***>>> updater calculating likelihood is " << u->getName() << std::endl;
	    last_ln_like = u->recalcLike();
        }
#if 0 //@POL this is what should be executed if not doing unimap 
	MCMCUpdaterShPtr u = *params_begin;
	last_ln_like = u->recalcLike();
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates and returns the log of the joint working prior by calling recalcWorkingPrior() for all parameters.
*/
double MCMCChainManager::recalcLnWorkingPrior() const
	{
	// Might need to add more checks here
	if (dirty)
		{
		throw XLikelihood("cannot call recalcLnWorkingPrior() for chain manager before calling finalize()");
		}
	
	//std::cerr << "\nIn MCMCChainManager::recalcLnWorkingPrior()..." << std::endl;//temp
	
	// Visit all updaters and let those who are prior stewards update their component of the joint prior.
	double ln_ref_dist = 0.0;
	for (MCMCUpdaterConstIter it = all_updaters.begin(); it != all_updaters.end(); ++it)
		{
		const boost::shared_ptr<MCMCUpdater> s = *it;
		if (s->isPriorSteward())
			{
			double this_ln_prior = s->recalcWorkingPrior();
			ln_ref_dist += this_ln_prior;
			//std::cerr << boost::str(boost::format("@@@@@@@@@@ this working prior = %g, cumulative = %g, name = '%s'") % this_ln_prior % ln_ref_dist % s->getName()) << std::endl;
			}
		}
	return ln_ref_dist;
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Replaces current (exponential) prior distribution for all internal edge length parameters with a new one having
|	the supplied mean, and returns the current log joint prior over all internal edge lengths.
*/
double MCMCChainManager::reviseInternalEdgeLenPrior(
  double mu)	/**< is the new mean edge length */
	{
	double sum_log_prior = 0.0;
	for (MCMCUpdaterIter edgelen_param_iter = edgelens_begin; edgelen_param_iter != edgelens_end; ++edgelen_param_iter)
		{
		EdgeLenParam * p = dynamic_cast<EdgeLenParam *>((*edgelen_param_iter).get());
		if (p->isInternalEdge())
			{
			//std::cerr << ">>> reviseInternalEdgeLenPriors: " << std::endl;	//@@@~
			//std::cerr << "      node number       = " << p->my_node->GetNodeNumber() << std::endl;
			//std::cerr << "      brlen             = " << p->getCurrValueFromModel() << std::endl;
			//std::cerr << "      prior address     = " << p->getPriorDist().get() << std::endl;
			
			if (sum_log_prior == 0.0)
				{
				//std::cerr << "      old prior         = " << p->getPriorDist()->GetDistributionDescription() << std::endl;
				// set prior just for first internal edge (all other internal edges use the same distribution object)
				p->getPriorDist()->SetMeanAndVariance(mu, mu*mu);
				}
			sum_log_prior += p->recalcPrior();
			
			//std::cerr << "      current prior     = " << p->getPriorDist()->GetDistributionDescription() << std::endl;
			//std::cerr << "      log-prior-density = " << p->recalcPrior() << std::endl;
			}
		}
	return sum_log_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces current (exponential) prior distribution for all internal edge length parameters with a new one having
|	the supplied mean, and returns the current log joint prior over all external edge lengths.
*/
double MCMCChainManager::reviseExternalEdgeLenPrior(
  double mu)	/**< is the new mean edge length */
	{
	double sum_log_prior = 0.0;
	for (MCMCUpdaterIter edgelen_param_iter = edgelens_begin; edgelen_param_iter != edgelens_end; ++edgelen_param_iter)
		{
		EdgeLenParam * p = dynamic_cast<EdgeLenParam *>((*edgelen_param_iter).get());
		if (p->isExternalEdge())
			{
			//std::cerr << ">>> reviseExternalEdgeLenPriors: " << std::endl;	//@@@~
			//std::cerr << "      node number       = " << p->my_node->GetNodeNumber() << std::endl;
			//std::cerr << "      brlen             = " << p->getCurrValueFromModel() << std::endl;
			//std::cerr << "      prior address     = " << p->getPriorDist().get() << std::endl;
			
			if (sum_log_prior == 0.0)
				{
				//std::cerr << "      old prior         = " << p->getPriorDist()->GetDistributionDescription() << std::endl;
				// set prior just for first external edge (all other external edges use the same distribution object)
				p->getPriorDist()->SetMeanAndVariance(mu, mu*mu);
				}
			sum_log_prior += p->recalcPrior();
			
			//std::cerr << "      current prior     = " << p->getPriorDist()->GetDistributionDescription() << std::endl;
			//std::cerr << "      log-prior-density = " << p->recalcPrior() << std::endl;
			}
		}
	return sum_log_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces current (exponential) prior distribution for all internal edge length parameters with a new one having
|	the supplied mean, and returns the current log joint prior over all edge lengths.
*/
double MCMCChainManager::reviseAllEdgeLenPriors(
  double mu)	/**< is the new mean edge length */
	{
	double sum_log_prior = 0.0;
	for (MCMCUpdaterIter edgelen_param_iter = edgelens_begin; edgelen_param_iter != edgelens_end; ++edgelen_param_iter)
		{
		EdgeLenParam * p = dynamic_cast<EdgeLenParam *>((*edgelen_param_iter).get());

		//std::cerr << ">>> reviseAllEdgeLenPriors: " << std::endl;	//@@@~
		//std::cerr << "      node number       = " << p->my_node->GetNodeNumber() << std::endl;
		//std::cerr << "      brlen             = " << p->getCurrValueFromModel() << std::endl;
		//std::cerr << "      prior address     = " << p->getPriorDist().get() << std::endl;
		//std::cerr << "      old prior         = " << p->getPriorDist()->GetDistributionDescription() << std::endl;
			
		p->getPriorDist()->SetMeanAndVariance(mu, mu*mu);
		sum_log_prior += p->recalcPrior();
			
		//std::cerr << "      current prior     = " << p->getPriorDist()->GetDistributionDescription() << std::endl;
		//std::cerr << "      log-prior-density = " << p->recalcPrior() << std::endl;
		}
	return sum_log_prior;
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Refreshes `last_ln_prior' by calling recalcPrior() for all parameters.
*/
void MCMCChainManager::refreshLastLnPrior()
	{
	last_ln_prior = 0.0;
	
	// Might need to add more checks here
	if (dirty)
		{
		throw XLikelihood("cannot call refreshLastLnPrior() for chain manager before calling finalize()");
		}
		
	// Visit all updaters and let those who are prior stewards update their component of the joint prior.
	for (MCMCUpdaterConstIter it = all_updaters.begin(); it != all_updaters.end(); ++it)
		{
		const boost::shared_ptr<MCMCUpdater> s = *it;
		if (s->isPriorSteward() && !s->isFixed())
			{
			double this_ln_prior = s->recalcPrior();
			last_ln_prior += this_ln_prior;
			//std::cerr << boost::str(boost::format("@@@@@@@@@@ Current prior = %g, Cumulative = %g, name = '%s'") % this_ln_prior % last_ln_prior % s->getName()) << std::endl;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds edge length parameters (if `separate_edgelen_params' is true) and model-specific parameters to the 
|	`edge_len_params' and `model_params' data members, respectively. Also adds edge length hyperparameters to the 
|   `edge_len_hyperparams' data member if any hyperpriors have been specified. The model-specific parameters are added 
|   by calling the model's createParameters() member function. This function does not call the finalize() function, 
|   however, in case other parameters need to be added after this function is called.
*/
void MCMCChainManager::addMCMCUpdaters(
  ModelShPtr m,						/**< is the substitution model */
  TreeShPtr t,						/**< is the tree */
  TreeLikeShPtr like,				/**< is the likelihood calculator */
  LotShPtr r,						/**< is the pseudo-random number generator */
  unsigned max_units,				/**< is the maximum number of slice sampler units to use in each update */
  unsigned weight,					/**< is the weight to be used for all parameters added by this function */
  int subset_pos)					/**< if 0 (first subset) or -1 (only subset), edge length parameters and hyperparams will be added */
	{
	//std::cerr << "\n***** model.getExternalEdgeLenPrior() upon entering MCMCChainManager::addMCMCUpdaters: use count = " << m->getExternalEdgeLenPrior().use_count() << std::endl;
	
	if (!like)
		{
		throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no TreeLikelihood object defined");
		}
	if (!t)
		{
		throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no tree defined");
		}
	if (!r)
		{
		throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no pseudorandom number generator defined");
		}
	if (!m)
		{
		throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no model defined");
		}
	if (subset_pos <= 0)
		{
		if (!m->getInternalEdgeLenPrior())
			{
			throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no internal edge length prior defined");
			}
		if (!m->getExternalEdgeLenPrior())
			{
			throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no external edge length prior defined");
			}
		}

	MCMCUpdaterIter iter;
	MCMCUpdaterVect edgelens;
	MCMCUpdaterVect edgelen_hyperparams;
	MCMCUpdaterVect parameters;

	// Ask the model to create the edge length parameters, edge length hyperparameters, and 
	// its own model-specific parameters.
	m->createParameters(t, edgelens, edgelen_hyperparams, parameters, subset_pos);

	if (subset_pos <= 0)
		{
		// Add the edge length parameters. If length 2, then first will be the master parameter for external edges
		// and the second will be the master parameter for the external edges. If length 1, then it will be the
		// edge length master parameter that governs all edges in the tree.
		for (iter = edgelens.begin(); iter != edgelens.end(); ++iter)
			{
			MCMCUpdaterShPtr p = (*iter);
			p->setWeight(weight);
			p->setMaxUnits(max_units);
			p->setModel(m);
			p->setTreeLikelihood(like);
			p->setLot(r);
			addEdgeLenParam(p);
			}

		// Add the edge length hyperparameters (if any were created). If length 2, then first will be the hyperparameter 
		// governing the mean of the external edge length prior and the second will be the hyperparameter governing
		// the mean of the internal edge length prior. If length 1, then it will be the hyperparameter governing the
		// mean of the prior for all edge lengths.
		for (iter = edgelen_hyperparams.begin(); iter != edgelen_hyperparams.end(); ++iter)
			{
			MCMCUpdaterShPtr p = (*iter);
			p->setWeight(weight);
			p->setMaxUnits(max_units);
			p->setModel(m);
			p->setTreeLikelihood(like);
			p->setLot(r);
			addEdgeLenHyperparam(p);
			}
		}

	//std::cerr << boost::str(boost::format("~~~ length of parameters vector is %d") % parameters.size()) << std::endl;

	for (iter = parameters.begin(); iter != parameters.end(); ++iter)
		{
		(*iter)->setWeight(weight);
		(*iter)->setMaxUnits(max_units);
		(*iter)->setModel(m);
		(*iter)->setTreeLikelihood(like);
		(*iter)->setLot(r);
		//std::cerr << boost::str(boost::format("~~~ calling MCMCChainManager::addModelParam for parameter %s") % (*iter)->getName()) << std::endl;
		addModelParam(*iter);
		}
	}

////////////////////////////////////////// inserted contents of mcmc_chain_manager.inl here /////////////////////////////

/*----------------------------------------------------------------------------------------------------------------------
|	The MCMCChainManager constructor does nothing.
*/
MCMCChainManager::MCMCChainManager() 
  : last_ln_like(0.0)
  , last_ln_prior(0.0)
  , dirty(true)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates the Robinson-Foulds distance between the supplied reference tree and the current tree (i.e. tree held by
|	first updater in `params_begin' vector).
*/
unsigned MCMCChainManager::calcRFDistance(TreeShPtr ref_tree) const
	{
	MCMCUpdaterShPtr u = *params_begin;
	return u->getTree()->robinsonFoulds(ref_tree);
	}

#define DEBUG_UPDATERS 1
#define DEBUG_UPDATERS_ENV_SENS 1

/*----------------------------------------------------------------------------------------------------------------------
|	For all updaters stored in `all_updaters', obtain the weight w and call the update fuction of the updater w times.
*/
void MCMCChainManager::verboseUpdateAllUpdaters()
	{
	for (MCMCUpdaterVect::iterator it = all_updaters.begin(); it != all_updaters.end(); ++it)
		{
		std::string nm = (*it)->getName();
		unsigned w = (*it)->getWeight();
		for (unsigned i = 0; i < w; ++i)
			{
			std::cerr << "########## updating " << nm << "..." << std::endl;
			(*it)->setSaveDebugInfo(true);
			(*it)->update();
			std::cerr << boost::str(boost::format("%s | %s") % nm % (*it)->getDebugInfo()) << std::endl;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	For all updaters stored in `all_updaters', obtain the weight w and call the update fuction of the updater w times.
*/
void MCMCChainManager::quietUpdateAllUpdaters()
	{
	for (MCMCUpdaterVect::iterator it = all_updaters.begin(); it != all_updaters.end(); ++it)
		{
		std::string nm = (*it)->getName();
		unsigned w = (*it)->getWeight();
		for (unsigned i = 0; i < w; ++i)
			{
			(*it)->update();
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	For all updaters stored in `all_updaters', obtain the weight w and call the update fuction of the updater w times.
*/
void MCMCChainManager::updateAllUpdaters()
	{
#	if DEBUG_UPDATERS
#		if DEBUG_UPDATERS_ENV_SENS
			const char * upEnv = getenv("DEBUG_UPDATE_ALL_UPDATERS");
			if (upEnv && upEnv[0] == '1')
				verboseUpdateAllUpdaters();
			else 
				quietUpdateAllUpdaters();
#		else
			verboseUpdateAllUpdaters();
#		endif
#	else
		quietUpdateAllUpdaters();
#	endif
	}



/*----------------------------------------------------------------------------------------------------------------------
|	Returns `all_updaters' data member as a const MCMCUpdaterVect reference. Used so that Python programs can iterate 
|	through all updaters using a construct such as this: for x in chainMgr.getAllUpdaters(): ...
*/
const MCMCUpdaterVect & MCMCChainManager::getAllUpdaters() const 
	{
	return all_updaters;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `moves' data member as a const MCMCUpdaterVect reference. Used so that Python programs can iterate 
|	through the stored moves using a construct such as this: for x in chainMgr.getMoves(): ...
*/
const MCMCUpdaterVect & MCMCChainManager::getMoves() const 
	{
	return moves;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `model_params' data member as a const MCMCUpdaterVect reference. Used so that Python programs can iterate 
|	through the model-specific parameters using a construct such as this: for x in chainMgr.getModelParams(): ...
*/
const MCMCUpdaterVect & MCMCChainManager::getModelParams() const 
	{
	return model_params;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `edge_len_params' data member as a const MCMCUpdaterVect reference. Used so that Python programs can iterate 
|	through the edge length parameters using a construct such as this: for x in chainMgr.getEdgeLenParams(): ...
*/
const MCMCUpdaterVect & MCMCChainManager::getEdgeLenParams() const 
	{
	return edge_len_params;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `edge_len_hyperparams' data member, which is a shared pointer to a HyperPriorParam.
|	Returns the `edge_len_hyperparams' data member as a const MCMCUpdaterVect reference. Used so that Python programs 
|   can iterate through the edge length hyperparameters using a construct such as this: 
|   for x in chainMgr.getEdgeLenHyperparams(): ...
*/
MCMCUpdaterVect MCMCChainManager::getEdgeLenHyperparams() const
	{
	return edge_len_hyperparams;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds a pointer to an MCMCUpdater-derived object to the `moves' vector. The `moves' vector is a vector of shared 
|	pointers to ensure that each move object is eventually destroyed.
*/
void	MCMCChainManager::addMove(MCMCUpdaterShPtr p) 
	{
	moves.push_back(p);
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds a pointer to an MCMCUpdater-derived object to the `model_params' vector. The `model_params' vector is a 
|	vector of shared pointers to ensure that each parameter object is eventually destroyed.
*/
void	MCMCChainManager::addModelParam(MCMCUpdaterShPtr p) 
	{
	model_params.push_back(p);
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds an EdgeLenParam shared pointer to the `edge_len_params' vector.
*/
void	MCMCChainManager::addEdgeLenParam(MCMCUpdaterShPtr p) 
	{
	edge_len_params.push_back(p);
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds a HyperPriorParam shared pointer to the `edge_len_hyperparams' vector.
*/
void MCMCChainManager::addEdgeLenHyperparam(MCMCUpdaterShPtr p) 
	{
	edge_len_hyperparams.push_back(p);
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the specified value in the `edge_len_hyperparams' vector. If there are two edge length hyperparameters, the 
|   first applies to external edge lengths and the second to internal edge lengths.
*/
void MCMCChainManager::setEdgeLenHyperparam(
  unsigned i,   /**< is the index of the value to set */
  double v)     /**< is the value to set */
	{
    PHYCAS_ASSERT(i < edge_len_hyperparams.size());
	edge_len_hyperparams[i]->sendCurrValueToModel(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current log likelihood. Each call to an update() member function of one of the updaters refreshes the
|	`last_ln_like' data member.
*/
double MCMCChainManager::getLastLnLike() const 
	{
	return last_ln_like;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Used by an MCMCUpdater-derived object to set the current log prior just before it returns from its update()
|	function. Ordinarily, the calcJointLnPrior() would be used instead to set `last_ln_prior', but for moves that need
|	to revert this function is used because the previous log prior was saved and can thus easily be restored.
*/
void MCMCChainManager::setLastLnPrior(double ln_prior) 
	{
	last_ln_prior = ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Used by an MCMCUpdater-derived object to set the current log likelihood just before it returns from its update()
|	function.
*/
void MCMCChainManager::setLastLnLike(double ln_like) 
	{
	last_ln_like = ln_like;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current log joint prior. Each call to the calcJointLnPrior() member function refreshes the 
|	`last_ln_prior' data member.
*/
double MCMCChainManager::getLastLnPrior() const 
	{
	return last_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears all vectors (all_updaters, moves, model_params, edge_len_params and edge_len_hyperparams).
*/
void MCMCChainManager::clear()
	{
	all_updaters.clear();
	moves.clear();
	model_params.clear();
	edge_len_params.clear();
	edge_len_hyperparams.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `chain_mgr' data member of all parameters so that they can call MCMCChainManager::getLnPrior when any one 
|	of them needs to compute the joint prior over all parameters (or wishes to refresh the MCMCChainManager's 
|	`last_ln_like' data member. The MCMCUpdater::`chain_mgr' data member is a weak_ptr (if a shared_ptr were used, we 
|	would have a cycle: MCMCChainManager holds a shared_ptr to the MCMCUpdater object, and MCMCUpdater would hold a 
|	shared_ptr to MCMCChainManager). The weak_ptr that is passed to the setChainManager member function of each 
|	MCMCUpdater object is created from a shared_ptr to this, which we can get because MCMCChainManager was derived from
|	boost::enable_shared_from_this<MCMCChainManager>. Using the shared_from_this() function presumes that there exists 
|	a shared_ptr already managing this. Thus, a MCMCChainManager object should not be created within C++ without using 
|	shared_ptr to manage it.
*/
void MCMCChainManager::finalize()
	{
	// Determine how big the all_updaters vector needs to be
	unsigned sz = (unsigned)(moves.size() + edge_len_params.size() + model_params.size() + edge_len_hyperparams.size());
	if (sz == 0)
		{
		throw XLikelihood("must add at least one updater object to the chain manager before calling finalize()");
		}

	// Clear all_updaters vector
	all_updaters.clear();
	all_updaters.resize(sz);

	// Add the moves first
	moves_end = std::copy(moves.begin(), moves.end(), all_updaters.begin());
	moves_begin = all_updaters.begin();

	// Next add the normal edge length parameters
	edgelens_end = std::copy(edge_len_params.begin(), edge_len_params.end(), moves_end);
	edgelens_begin = moves_end;     

	// Next add the edge length hyperparameters (if there are any)
	hyperparams_end = std::copy(edge_len_hyperparams.begin(), edge_len_hyperparams.end(), edgelens_end);
	hyperparams_begin = edgelens_end;

	// Next add the remaining model parameters
	model_params_end = std::copy(model_params.begin(), model_params.end(), hyperparams_end);
	model_params_begin = hyperparams_end;

	params_begin = edgelens_begin;
	params_end = model_params_end;
	
	// Call each updater's setChainManager member function, supplying a shared pointer to this object
	// This allows each updater to ask the MCMCChainManager to calculate the joint prior over all
	// parameters when needed for computing the posterior density during slice sampling
	//std::for_each(all_updaters.begin(), all_updaters.end(), 
	//	boost::lambda::bind(&MCMCUpdater::setChainManager, *boost::lambda::_1, ChainManagerWkPtr(shared_from_this())));
	ChainManagerWkPtr wptr(shared_from_this());
	for (MCMCUpdaterIter uit = all_updaters.begin(); uit != all_updaters.end(); ++uit)
		{
		//std::string nm = (*uit)->getName();
		//unsigned uc1 = (unsigned)(*uit).use_count();
		(*uit)->setChainManager(wptr);
		//unsigned uc2 = (unsigned)(*uit).use_count();
		//std::cerr << ">>>> updater " << nm << " (" << uc1 << " --> " << uc2 << ")" << std::endl;
		}
	
	for (MCMCUpdaterIter piter = params_begin; piter != params_end; ++piter)
		{
		// Call each parameter's recalcPrior member function to make sure their curr_ln_prior data members are
		// up to date
		(*piter)->recalcPrior();
		}

	//@POL Wondering why we need all these vectors? It is convenient to have one vector of updaters, but not all of 
	// these updaters are necessarily created in MCMCChainManager::addMCMCUpdaters. When finalize() is called, we 
	// can be sure that all the updaters that are going to be defined have been defined, and at that point we can 
	// quickly build the vector all_updaters. In the future, if it turns out we never need moves, model_params, or
	// edge_len_params after this point, we can just clear them at this point because we have iterators now set to
	// the beginning and end of each section of the all_updaters vector anyway.

	dirty = false;
	
#if 0
	std::cerr << std::endl;
	std::cerr << "##############################################" << std::endl;
	std::cerr << "Here are the updaters that can compute priors:" << std::endl;
	for (MCMCUpdaterIter it = all_updaters.begin(); it != all_updaters.end(); ++it)
		{
		MCMCUpdaterShPtr s = *it;
		if (s->computesUnivariatePrior())
			std::cerr << boost::str(boost::format("  '%s' (univariate prior)") % s->getName()) << std::endl;
		else if (s->computesMultivariatePrior())
			std::cerr << boost::str(boost::format("  '%s' (multivariate prior)") % s->getName()) << std::endl;
		else
			std::cerr << boost::str(boost::format("  '%s' (does not compute prior)") % s->getName()) << std::endl;
		}
		
	std::cerr << std::endl;
	std::cerr << "##############################################" << std::endl;
	std::cerr << "Here is the list of just the parameters:" << std::endl;
	for (MCMCUpdaterIter p = params_begin; p != params_end; ++p)
		{
		MCMCUpdaterShPtr s = *p;
		if (s->computesUnivariatePrior())
			std::cerr << boost::str(boost::format("  '%s' (univariate prior)") % s->getName()) << std::endl;
		else if (s->computesMultivariatePrior())
			std::cerr << boost::str(boost::format("  '%s' (multivariate prior)") % s->getName()) << std::endl;
		else
			std::cerr << boost::str(boost::format("  '%s' (does not compute prior)") % s->getName()) << std::endl;
		}
#endif
	}

// Using USE_PARTIAL_SUM_FUNCTOR is slightly faster, but makes the assumption that the edge length prior
// distribution is Exponential, which is not good. This assumption is not necessary, but without it
// speed suffers and the alternative is preferable.
// 
// Note: now that separate priors have been installed for internal and external edges, the PartialPriorSum 
// approach is even less useful and has been abandoned. I've left the code here in case in case something
// like this needs to be done elsewhere in the future.
//#define USE_PARTIAL_SUM_FUNCTOR

//#if defined(USE_PARTIAL_SUM_FUNCTOR)
//struct PartialPriorSum : public std::unary_function<double, void>
//	{
//	double sum;
//	ExponentialDistribution & distr;
//
//	PartialPriorSum(ExponentialDistribution & p) : sum(0.0), distr(p) 
//		{}
//
//	void operator()(const double & nxt)
//		{
//	try 
//		{
//		sum += distr.GetLnPDF(nxt);
//		}
//	catch(XProbDist & x)
//		{
//		sum += distr.GetRelativeLnPDF(nxt);
//		}
//		}
//
//	double result()
//		{
//		return sum;
//		}
//	};
//#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the unnormalized prior for the supplied edge length using the prior for external edge lengths.
*/
double MCMCChainManager::calcExternalEdgeLenPriorUnnorm(
  double v) const
	{
    // The first edge length parameter always governs external edge lengths
    PHYCAS_ASSERT((edgelens_end - edgelens_begin) >= 1);
	MCMCUpdaterShPtr u = *edgelens_begin;
	ProbDistShPtr p = u->getPriorDist();
	double tmp = 0.0;
	try 
		{
		tmp = p->GetLnPDF(v);
		}
	catch(XProbDist &)
		{
		PHYCAS_ASSERT(0);
		}
    return tmp;
    }
	
/*----------------------------------------------------------------------------------------------------------------------
|	Informs the EdgeLenMasterParam objects of a new minimum sample size `n' needed for constructing split-specific edge
|	length working priors. If `n' is zero, this effectively turns off the use of split-specific edge length working
|	priors, in which case a single working prior will be constructed and used for all edge lengths. This is only relevant
|	if a steppingstone analysis is being performed.
*/
void MCMCChainManager::setMinSSWPSampleSize(
  unsigned n)	/**< is the minimum sample size needed to create a split-specific edge length working prior */
	{
	// If there are two edge_length_parameters then the first is for external and second is for internal edges
	// (see MCMCChainManager::addMCMCUpdaters). If there is only one, it handles all edge lengths, so in this
	// case (external edge length) we know we should go with the first one.
    const MCMCUpdaterVect & edge_length_params = getEdgeLenParams();
	EdgeLenMasterParam * p = dynamic_cast<EdgeLenMasterParam *>(edge_length_params[0].get());
	
	// If p is NULL, then there are no EdgeLenMasterParam objects, which means that the topology is fixed and
	// an EdgeLenParam has been assigned to each edge. In this case we can just exit.
	if (p)
		{
		if (n > 0)
			{
			p->useEdgeSpecificWorkingPriors(true);
			p->setMinWorkingPriorSampleSize(n);
			}
		else 
			p->useEdgeSpecificWorkingPriors(false);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the log of the working prior for the supplied external edge length.
*/
double MCMCChainManager::calcExternalEdgeLenWorkingPrior(
  const TreeNode & nd,	/**< is the node for which the working prior log density should be computed */
  double v) const		/**< is the edge length for which the working prior log density should be computed */
	{
	// If there are two edge_length_parameters, first is for external and second is for internal edges
	// (see MCMCChainManager::addMCMCUpdaters). If only one, it handles all edge lengths, so in this
	// case (external edge length) we know we should go with the first one.
    const MCMCUpdaterVect & edge_length_params = getEdgeLenParams();
	EdgeLenMasterParam * p = dynamic_cast<EdgeLenMasterParam *>(edge_length_params[0].get());
	PHYCAS_ASSERT(p);
	return p->lnWorkingPriorOneEdge(nd, v);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Computes the log of the working prior for the supplied internal edge length.
*/
double MCMCChainManager::calcInternalEdgeLenWorkingPrior(
  const TreeNode & nd,		/**< is the node for which the working prior log density should be computed */
  double v) const			/**< is the edge length for which the working prior log density should be computed */
	{	
	// If there are two edge_length_parameters, first is for external and second is for internal edges
	// (see MCMCChainManager::addMCMCUpdaters). If only one, it handles all edge lengths.
	const MCMCUpdaterVect & edge_length_params = getEdgeLenParams();
	if (edge_length_params.size() == 2)
		{
		EdgeLenMasterParam * p = dynamic_cast<EdgeLenMasterParam *>(edge_length_params[1].get());
		PHYCAS_ASSERT(p);
		return p->lnWorkingPriorOneEdge(nd, v);
		//return edge_length_params[1]->lnWorkingPriorOneEdge(nd, false);
		}
	else 
		{
		EdgeLenMasterParam * p = dynamic_cast<EdgeLenMasterParam *>(edge_length_params[0].get());
		PHYCAS_ASSERT(p);
		return p->lnWorkingPriorOneEdge(nd, v);
		//return edge_length_params[0]->lnWorkingPriorOneEdge(nd, false);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the unnormalized prior for the supplied edge length using the prior for internal edge lengths.
*/
double MCMCChainManager::calcInternalEdgeLenPriorUnnorm(double v) const
	{
    // The first edge length parameter governs internal edge lengths if there is only one edge length parameter.
    // If there are two edge length parameters, then the second one always governs internal edge lengths
	bool two_edgelen_priors = (edgelens_end - edgelens_begin == 2);
    MCMCUpdaterShPtr u = (two_edgelen_priors ? (*(edgelens_begin + 1)) : (*edgelens_begin));
	ProbDistShPtr p = u->getPriorDist();
	double tmp = 0.0;
	try 
		{
		tmp = p->GetLnPDF(v);
		}
	catch(XProbDist &)
		{
		PHYCAS_ASSERT(0);
		}
    return tmp;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Visits each parameter and calls its getLnPrior function to compute the joint log prior over all parameters. The
|	accumulate algorithm requires its fourth argument to be a binary functor, for which an unnamed boost lambda function
|	object is used. This function refreshes the value of the data member `last_ln_prior'.
*/
double MCMCChainManager::calcJointLnPrior()
	{
	//	if (all_updaters.empty())
	//		{
	//		throw XLikelihood("should not call MCMCChainManager::calcJointLnPrior before calling MCMCChainManager::finalize");
	//		}
	//
	//	last_ln_prior = std::accumulate(params_begin, params_end, 0.0, 
	//		boost::lambda::_1 += boost::lambda::bind(&MCMCUpdater::getLnPrior, *boost::lambda::_2));
	refreshLastLnPrior();
	return last_ln_prior;
	}

#if 0
/*----------------------------------------------------------------------------------------------------------------------
|	Uses the two edge length parameters in the `all_updaters' vector to compute the log-prior for each value in the 
|	supplied vector of edge lengths, returning the sum. Useful for efficiency purposes for moves that modify only a few 
|	edge lengths.
*/
double MCMCChainManager::partialEdgeLenPrior(const std::vector<double> & edge_len_vect) const
	{
    // To do next: replace all calls to this function (partialEdgeLenPrior) with calls to new functions
    // named calcInternalEdgeLenPriorUnnorm and calcExternalEdgeLenPriorUnnorm. There are only a couple of places
    // where the edge_len_vect vector has more than one element, so this function is not saving much
    // and we don't know which of the edge lengths supplied belong to internal or external edges.

	if (edgelens_end - edgelens_begin == 0)
		{
		throw XLikelihood("partialEdgeLenPrior requires a master edge length prior to have been added to the chain manager");
		}

	MCMCUpdaterShPtr u = *edgelens_begin;     
	ProbDistShPtr p = u->getPriorDist();

//#if defined(USE_PARTIAL_SUM_FUNCTOR)
//	ExponentialDistribution * d = dynamic_cast<ExponentialDistribution *>(p.get());
//	PHYCAS_ASSERT(d);
//	//@POL assuming edge length distribution is Exponential for now (allows GetLnPDF to be inlined)
//	return std::for_each(edge_len_vect.begin(), edge_len_vect.end(), PartialPriorSum(*d)).result();
//#else
	double partial_prior_sum = 0.0;
	for (std::vector<double>::const_iterator it = edge_len_vect.begin(); it != edge_len_vect.end(); ++it)
		{
		double tmp = 0.0;
		try 
			{
			tmp = p->GetLnPDF(*it);
			}
		catch(XProbDist &)
			{
		    PHYCAS_ASSERT(0);
			}
		partial_prior_sum += tmp;
		}
	return partial_prior_sum;
//#endif
	}
#endif

}	// namespace phycas
