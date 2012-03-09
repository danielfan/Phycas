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

#if ! defined(GTR_MODEL_HPP)
#define GTR_MODEL_HPP

#include "phycas/src/states_patterns.hpp"
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <vector>
#include <numeric>
#include <algorithm>
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/mcmc_param.hpp"
#include "phycas/src/basic_cdf.hpp"
#include "phycas/src/q_matrix.hpp"
#include "phycas/src/likelihood_models.hpp"

namespace phycas{
	
/*----------------------------------------------------------------------------------------------------------------------
|	Specialization of the base class Model that represents the General Time Reversible (GTR) model.
*/
class GTR: public Model
	{
	public:
									GTR();
									~GTR()
										{
										//std::cerr << "GTR dying..." << std::endl;
										}

		virtual std::string			getModelName() const;
		virtual void				createParameters(TreeShPtr t, MCMCUpdaterVect & edgelens, MCMCUpdaterVect & edgelen_hyperparams, MCMCUpdaterVect & parameters, int subset_pos);
		double						calcUniformizationLambda() const;
        double					    calcLMat(double * * lMat) const;
        double					    calcUMat(double * * uMat) const;
		void						calcPMat(double * * pMat, double edgeLength) const;

        void						fixRelRates();
		void						freeRelRates();

        std::vector<double>			getRelRates();
		void						setRelRates(const std::vector<double> & rates);	
		void						setRelRateUnnorm(unsigned param_index, double value);
		double						getRelRateUnnorm(unsigned param_index);
		
		void						setRelRateParamPrior(ProbDistShPtr d);
		ProbDistShPtr				getRelRateParamPrior();
		
		void						setRelRatePrior(MultivarProbDistShPtr d);
		MultivarProbDistShPtr		getRelRatePrior();

        void						setNucleotideFreqs(double freqA, double freqC, double freqG, double freqT);
        void						setAllFreqsEqual();
		virtual void				setStateFreqUnnorm(unsigned param_index, double value);
        virtual void                setStateFreqsUnnorm(const std::vector<double> & values);

        void						setStateFreqPrior(MultivarProbDistShPtr d);
		MultivarProbDistShPtr		getStateFreqPrior();

        void						setStateFreqParamPrior(ProbDistShPtr d);
		ProbDistShPtr				getStateFreqParamPrior();

        virtual std::string			paramHeader() const;
		virtual std::string			paramReport(unsigned ndecimals, bool include_edgelen_hyperparams) const;
		double						calcTRatio();

		void						beagleGetStateFreqs(std::vector<double> & freqs);
		void						beagleGetEigenValues(std::vector<double> & eigenValues);
		void						beagleGetEigenVectors(std::vector<double> & eigenVectors);
		void						beagleGetInverseEigenVectors(std::vector<double> & inverseEigenVectors);
		double						beagleGetEdgelenScaler();

	protected:

		MultivarProbDistShPtr		rel_rate_prior;		    /**< The prior distribution governing each relative rate (usually a gamma distribution with scale 1 and shape equal to the desired Dirichlet parameter) */
		ProbDistShPtr				rel_rate_param_prior;	/**< The joint prior distribution governing all six relative rates */
		ProbDistShPtr				freq_param_prior;	    /**< The prior distribution governing each frequency parameter (usually a gamma distribution with scale 1 and shape equal to the desired Dirichlet parameter; used if frequencies are updated separately by slice sampling) */
    	MultivarProbDistShPtr		freq_prior;	            /**< The prior distribution governing the vector of frequencies (used if frequencies are updated jointly by StateFreqMove) */
		bool						rel_rates_fixed;	    /**< If true, the relative rate values will not change during MCMC updates */
		mutable MCMCUpdaterVect		rel_rate_params;	    /**< A vector containing copies of all six relative rate parameters (saved so that fixed/free status can be changed) */
		mutable QMatrix				q_matrix;			    /**< A QMatrix object used to compute transition probabilities */

		string_vect_t				relrate_name;			/**< Holds names of the 6 relative rates (e.g. rAC, rAG, rAT, rCG, rCT, rGT) used as headers in the param file */
		string_vect_t				freq_name;				/**< Holds names of the 4 state frequencies (e.g. freqA, freqC, freqG, freqT) used as headers in the param file */
		
		// Below here are quantities that directly affect likelihood calculations and which should increment time_stamp when modified
		std::vector<double>			rel_rates;			    /**< A vector containing the six relative rates */
	};

typedef boost::shared_ptr<GTR> GTRShPtr;

} // namespace phycas

//#include "phycas/src/gtr_model.inl"

#endif

