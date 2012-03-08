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

#if ! defined(LIKELIHOOD_MODELS_HPP)
#define LIKELIHOOD_MODELS_HPP

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

namespace phycas{
	
//@POL I am not entirely happy about the way the Model hierarchy is structured. Might be better to have
// a Model class that derives from several policy classes supplied as template arguments. One policy
// might govern state frequencies, another relative rates, a third the calculation of transition matrices,
// and a fourth rate heterogeneity. This would avoid, for example, storing base frequencies in two places
// (Model has one copy and the QMatrix data member in GTR must store a second copy). It would also avoid,
// I think, having to have pointers to GTR and HKY objects just so functions not in the Model base class
// can be called.

/*----------------------------------------------------------------------------------------------------------------------
|	Base class for all models used by Phycas. This class defines pure virtual functions, and thus only classes derived
|	from Model can be instantiated.
*/
class Model	{
#if defined(FLEXCAT_MOVE)		
	friend class NCatMove;
#endif

	public:

	// Destructor
		virtual							~Model();

		// Utilities
		virtual void					releaseUpdaters();
		virtual std::string				getModelName() const = 0;
		virtual void					calcPMat(double * * pMat, double edgeLength) const = 0;
		void							calcPMatrices(double * * * pMat, const double * edgeLength, unsigned numRates) const;
		virtual std::string				lookupStateRepr(int state) const;
        virtual void					createParameters(TreeShPtr t, MCMCUpdaterVect & edgelens, MCMCUpdaterVect & edgelen_hyperparams, MCMCUpdaterVect & parameters, int subset_pos);
        virtual void					buildStateList(state_list_t &, state_list_pos_t &) const;

		virtual std::string				paramHeader() const;
		virtual std::string				paramReport(unsigned ndecimals, bool include_edgelen_hyperparams) const;

		virtual double					calcUniformizationLambda() const = 0;
		virtual double					calcLMat(double * * lMat) const = 0;
		virtual double					calcUMat(double * * uMat) const = 0;
		
		virtual void					Clear();

		// Query functions
		bool							isCodonModel() const;
		unsigned						getTimeStamp() const
											{return time_stamp;}

		// Member functions related to relative state frequencies
		unsigned						getNumStates() const;
		const std::vector<double>	&	getStateFreqs() const;
		virtual void					setNucleotideFreqs(double freqA, double freqC, double freqG, double freqT);
		virtual void					setAllFreqsEqual();
        virtual void                    setStateFreqsUnnorm(const std::vector<double> & values);
		virtual void					setStateFreqUnnorm(unsigned param_index, double value);
        double						    getStateFreqUnnorm(unsigned param_index);
		void							normalizeFreqs();

		// Member functions related to relative rates
		unsigned						getNRatesTotal() const;
		unsigned						getNGammaRates() const;
		void							setNGammaRates(unsigned nRates);
		const std::vector<double>	&	getGammaRateProbs() const;
		void							setAllGammaRateProbsEqual();
		void							recalcRatesAndProbs(std::vector<double> & rates, std::vector<double> & probs) const;
		void							recalcGammaRatesAndBoundaries(std::vector<double> & rates, std::vector<double> & boundaries) const;
		
		// Member functions related to state frequencies
        bool							stateFreqsFixed() const;
		void							fixStateFreqs();
		void							freeStateFreqs();
		
		// Member functions related to proportion of invariable sites
        bool							isPinvarModel();
        void							setPinvarModel();
		void							setNotPinvarModel();
		bool							pinvarFixed() const;
		void							fixPinvar();
		void							freePinvar();
		double							getPinvar();
		void							setPinvar(double pinv);
		void							setPinvarPrior(ProbDistShPtr d);
		ProbDistShPtr					getPinvarPrior();

		// Member functions related to discrete gamma shape
		bool							shapeFixed() const;
		void							fixShape();
		void							freeShape();
		double							getShape();
		void							setShape(double alpha);
		void							setPriorOnShapeInverse(bool invert);
		void							setDiscreteGammaShapePrior(ProbDistShPtr d);
		ProbDistShPtr					getDiscreteGammaShapePrior();

		// Member functions related to edge lengths
        bool							edgeLengthsFixed() const;
		void							fixEdgeLengths();
		void							freeEdgeLengths();
		bool							edgeLenHyperParamFixed() const;
		void							fixEdgeLenHyperprior();
		void							freeEdgeLenHyperprior();
		void							setEdgeLenHyperPrior(ProbDistShPtr d);
		//void							setEdgeLenPrior(ProbDistShPtr d);
		void							setExternalEdgeLenPrior(ProbDistShPtr d);
		void							setInternalEdgeLenPrior(ProbDistShPtr d);
		ProbDistShPtr					getEdgeLenHyperPrior();
		//ProbDistShPtr					getEdgeLenPrior();
		bool							hasEdgeLenHyperPrior();
		void							setInternalEdgelenHyperparam(double mu) {internal_edgelen_hyperparam = mu;}
		double							getInternalEdgelenHyperparam() {return internal_edgelen_hyperparam;}

		void							setExternalEdgelenHyperparam(double mu) {external_edgelen_hyperparam = mu;}
		double							getExternalEdgelenHyperparam() {return external_edgelen_hyperparam;}
		ProbDistShPtr					getExternalEdgeLenPrior();
		ProbDistShPtr					getInternalEdgeLenPrior();
        void                            separateInternalExternalEdgeLenPriors(bool separate);
        bool                            isSeparateInternalExternalEdgeLenPriors() const;
        void                            setEdgeSpecificParams(bool param_for_each_edgelen);

		// Utility functions
		void flattenTwoDMatrix(VecDbl & p, double * * twoDarr, unsigned dim) const;


		virtual void					beagleGetStateFreqs(std::vector<double> & freqs) {}
		virtual void					beagleGetEigenValues(std::vector<double> & eigenValues){}
		virtual void					beagleGetEigenVectors(std::vector<double> & eigenVectors){}
		virtual void					beagleGetInverseEigenVectors(std::vector<double> & inverseEigenVectors) {}
		virtual double					beagleGetEdgelenScaler() {return 0;}
	
	
	// Python-specific utilities
#if defined(PYTHON_ONLY)
#	if defined(USING_NUMARRAY)
		boost::python::numeric::array	getPMatrix(double edgeLength) const;
#	else
		VecDbl getPMatrix(double edgeLength) const;
#	endif
#endif

protected:

	// Constructor
									Model(unsigned numStates);

protected:

	int								subset_index;				/**< The index of the partition subset to which this parameter belongs (set in createParameters function and mostly used in forming parameter names) */
	double							internal_edgelen_hyperparam;	/**< The internal edge length hyperparameter */
	double							external_edgelen_hyperparam;	/**< The external edge length hyperparameter */
	ProbDistShPtr					edgeLenHyperPrior;			/**< The prior distribution governing the mean of the edge length prior if a hierarchical model is used */
	ProbDistShPtr					internalEdgeLenPrior;		/**< The prior distribution governing internal edge lengths */
	ProbDistShPtr					externalEdgeLenPrior;		/**< The prior distribution governing external edge lengths */
    bool                            separate_edgelen_params;    /**> If true, each edge length in the (fixed) tree topology will have its own updater */
    bool                            separate_int_ext_edgelen_priors;    /**> If true, internal edge lengths have a different prior than external edge lengths */
	std::vector<double>				state_freqs;				/**< A vector of relative state frequencies (length is `num_states') */
	std::vector<std::string>		state_repr;					/**< A vector of strings representing the states allowed by this model */
	bool							state_freq_fixed;			/**< If true, the values in `state_freq_params' will not change during MCMC updates */
	mutable MCMCUpdaterVect			freq_params;				/**< Vector of shared pointers to the state frequency parameters (need to retain pointers to these so that the fixed/free status can be changed) */	
    mutable MCMCUpdaterVect		    edgelen_hyper_params;		/**< Vector of shared pointers to the edge length hyperparameters (need to retain a pointer so that the fixed/free status can be changed) */
    mutable MCMCUpdaterVect			edgelen_params;				/**< Vector of shared pointers to the edge length parameters (need to retain pointers to these so that their fixed/free status can be changed) */
	bool							edge_lengths_fixed;			/**< If true, the value of the edge lengths will not change during MCMC updates */
	bool							edgelen_hyperprior_fixed;	/**< If true, the value of the edge length hyperprior will not change during MCMC updates */
	bool							pinvar_fixed;				/**< If true, the value of pinvar will not change during MCMC updates */
	mutable MCMCUpdaterShPtr		pinvar_param;				/**< Shared pointer to the proportion of invariable sites parameter (need to retain a pointer so that the fixed/free status can be changed) */
	ProbDistShPtr					pinvar_prior;				/**< The prior distribution governing the proportion of invariable sites parameter */
	bool							gamma_shape_fixed;			/**< If true, the value of gamma_shape will not change during MCMC updates */
	mutable MCMCUpdaterShPtr		gamma_shape_param;			/**< Shared pointer to the gamma shape parameter (need to retain a pointer so that the fixed/free status can be changed) */
	ProbDistShPtr					gamma_shape_prior;			/**< The prior distribution governing the discrete gamma shape parameter */
	bool							invert_shape;				/**< If true, gamma_shape_param will hold inverse of shape rather than shape itself */
	CDF								cdf;						/**< Provides cumulative gamma distribution function */
	
	unsigned						time_stamp;					/**< Provides a way for methods that depend on the model to know if the model has changed */
	
	// Below here are quantities that directly affect likelihood calculations and which should increment time_stamp when modified
	unsigned						num_states;					/**< The number of states (e.g. 4 for DNA) */
	unsigned						num_gamma_rates;			/**< The number of discrete gamma rate categories. If greater than 1, the model becomes a discrete gamma rate heterogeneity ("G") model. */
	std::vector<double>				state_freq_unnorm;			/**< A vector of quantities that yield the values in `state_freqs' when normalized using the normalizeFreqs member function (length is `num_states') */
	mutable std::vector<double>		gamma_rates_unnorm;			/**< A vector of quantities that yield the relative rates when normalized in recalcRatesAndProbs (length is `num_gamma_rates') */
	mutable std::vector<double>		gamma_rate_probs;			/**< A vector of probabilities that a site falls in any given rate category (length is `num_gamma_rates') */
	double							gamma_shape;				/**< Used for discrete gamma rate heterogeneity */
	double							pinvar;						/**< The proportion of invariable sites. If non-zero, the model becomes an invariable-sites ("I") model. */
	bool							is_codon_model;				/**< If true, nucleotide states will be interpreted as triplets when creating TipData structures for tree */
	bool							is_pinvar_model;			/**< If true, a parameter for pinvar will be added to MCMC analysis (pinvar_fixed determines whether it is updated or not) */
	
};

typedef boost::shared_ptr<Model> ModelShPtr;
typedef std::vector<ModelShPtr> ModelVect;

/*----------------------------------------------------------------------------------------------------------------------
|	Specialization of the base class Model that represents the simple irreversible model.
*/
class Irreversible: public Model
{
public:
    Irreversible();
    ~Irreversible()
    {
        //std::cerr << "~Irreversible dying..." << std::endl;
    }
    
    virtual std::string		getModelName() const;
    double					calcUniformizationLambda() const;
    void					calcPMat(double * * pMat, double edgeLength) const;
    double					calcLMat(double * * lMat) const;
    double					calcUMat(double * * uMat) const;
    virtual std::string		paramHeader() const;
    virtual std::string		paramReport(unsigned ndecimals, bool include_edgelen_hyperparams) const;
    
    void					fixScalingFactor();
    void					freeScalingFactor();
    void                    setScalingFactor(double sf);
    double                  getScalingFactor();
    
    void                    setScalingFactorPrior(ProbDistShPtr d);
    ProbDistShPtr           getScalingFactorPrior();
    
    void                    setGainOnly();
    void                    setLossOnly();
    
protected:
	ProbDistShPtr           scaling_factor_prior;		/**< The prior distribution governing scaling_factor */
    double                  scaling_factor;            /**< Scaling factor parameter. Allows rate of evolution to depart from that implied by the edge lengths (i.e. scaling_factor = 1.0) */
    bool                    scaling_factor_fixed;      /**< If true (default), scaling_factor parameter will not be modified during MCMC analyses */
    bool                    root_present;   /**< If true (default), root state is assumed to be 1 and only losses are allowed; if false, root state is assumed to be 0 and only gains are allowed */
};

typedef boost::shared_ptr<Irreversible> IrreversibleShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Specialization of the base class Model that represents the simple 2-state Markov model.
*/
class Binary: public Model
{
public:
    Binary();
    ~Binary()
    {
        //std::cerr << "Binary dying..." << std::endl;
    }
    
    virtual std::string		getModelName() const;
    double					calcUniformizationLambda() const;
    void					calcPMat(double * * pMat, double edgeLength) const;
    double					calcLMat(double * * lMat) const;
    double					calcUMat(double * * uMat) const;
    virtual std::string		paramHeader() const;
    virtual std::string		paramReport(unsigned ndecimals, bool include_edgelen_hyperparams) const;
    
    void					fixScalingFactor();
    void					freeScalingFactor();
    void                    setScalingFactor(double sf);
    double                  getScalingFactor();
        
    void                    setScalingFactorPrior(ProbDistShPtr d);
    ProbDistShPtr           getScalingFactorPrior();

    void					fixKappa();
    void					freeKappa();
    void                    setKappa(double sf);
    double                  getKappa();

    void                    setKappaPrior(ProbDistShPtr d);
    ProbDistShPtr           getKappaPrior();
    
protected:
	ProbDistShPtr   scaling_factor_prior;   /**< The prior distribution governing scaling_factor */
    double          scaling_factor;         /**< Allows rate of evolution to depart from that implied by the edge lengths */
    bool            scaling_factor_fixed;   /**< If true (default), scaling_factor parameter will not be modified during MCMC analyses */

	ProbDistShPtr   kappa_prior;            /**< The prior distribution governing kappa, the forward/reverse rate ratio */
    double          kappa;                  /**< Forward/reverse rate ratio parameter */
    bool            kappa_fixed;            /**< If true (default), kappa will not be modified during MCMC analyses */
};

typedef boost::shared_ptr<Binary> BinaryShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Specialization of the base class Model that represents the Jukes and Cantor (1969) model.
*/
class JC: public Model
{
public:
    JC();
    ~JC()
    {
        //std::cerr << "JC dying..." << std::endl;
    }
    
    virtual std::string		getModelName() const;
    double					calcUniformizationLambda() const;
    void					calcPMat(double * * pMat, double edgeLength) const;
    double					calcLMat(double * * lMat) const;
    double					calcUMat(double * * uMat) const;
    virtual std::string		paramHeader() const;
    virtual std::string		paramReport(unsigned ndecimals, bool include_edgelen_hyperparams) const;
};

typedef boost::shared_ptr<JC> JCShPtr;
    
/*----------------------------------------------------------------------------------------------------------------------
|	Specialization of the base class Model that represents the Hasegawa-Kishino-Yano (1985) model.
*/
class HKY: public Model
	{
	public:
									HKY();
									~HKY();
									
		virtual std::string			getModelName() const;
        virtual void				createParameters(TreeShPtr t, MCMCUpdaterVect & edgelens, MCMCUpdaterVect & edgelen_hyperparams, MCMCUpdaterVect & parameters, int subset_pos);
		double						calcUniformizationLambda() const;
        double					    calcLMat(double * * lMat) const;
        double					    calcUMat(double * * uMat) const;
        void						calcPMat(double * * pMat, double edgeLength) const;

        void						fixKappa();
		void						freeKappa();

        double						getKappa();
		void						setKappa(double k);

        void						setKappaPrior(ProbDistShPtr d);
		ProbDistShPtr				getKappaPrior();

        void						setKappaFromTRatio(double tratio);
		double						calcTRatio();

		void						setStateFreqPrior(MultivarProbDistShPtr d);
		MultivarProbDistShPtr		getStateFreqPrior();

		void						setStateFreqParamPrior(ProbDistShPtr d);
		ProbDistShPtr				getStateFreqParamPrior();

		virtual std::string			paramHeader() const;
		virtual std::string			paramReport(unsigned ndecimals, bool include_edgelen_hyperparams) const;

protected:
	
	ProbDistShPtr				kappa_prior;		/**< The prior distribution governing kappa */
	ProbDistShPtr				freq_param_prior;	/**< The prior distribution governing each frequency parameter (used if frequencies are updated separately by slice sampling) */
	MultivarProbDistShPtr		freq_prior;	        /**< The prior distribution governing the vector of frequencies (used if frequencies are updated jointly by StateFreqMove) */
	bool						kappa_fixed;		/**< If true, the value of kappa will not change during MCMC updates */
	mutable MCMCUpdaterShPtr	kappa_param;		/**< Copy of the kappa parameter (saved so that fixed/free status can be changed) */

	string_vect_t				freq_name;				/**< Holds names of the 4 state frequencies (e.g. freqA, freqC, freqG, freqT) used as headers in the param file */

	// Below here are quantities that directly affect likelihood calculations and which should increment time_stamp when modified
	double						kappa;				/**< The transition/transversion rate ratio */
	};

typedef boost::shared_ptr<HKY> HKYShPtr;

} // namespace phycas

#include "phycas/src/likelihood_models.inl"

#endif

