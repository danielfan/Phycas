#if ! defined(LIKELIHOOD_MODELS_HPP)
#define LIKELIHOOD_MODELS_HPP

#include "pyphy/src/states_patterns.hpp"
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <vector>
#include <numeric>
#include <algorithm>
#include "pyphy/src/xlikelihood.hpp"
#include "pyphy/src/mcmc_param.hpp"
#include "pyphy/src/basic_cdf.hpp"
#include "pyphy/src/q_matrix.hpp"

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
	friend class NCatMove;

	public:

	// Destructor
		virtual							~Model();

		// Utilities
		virtual std::string				getModelName() const = 0;
		virtual void					calcPMat(double * * pMat, double edgeLength) const = 0;
		void							calcPMatrices(double * * * pMat, const double * edgeLength, unsigned numRates) const;
		virtual std::string				lookupStateRepr(int state) const;
		virtual void					createParameters(TreeShPtr t, MCMCUpdaterVect & edgelens, MCMCUpdaterShPtr & edgelen_hyperparam, MCMCUpdaterVect & parameters, bool separate_edgelens) const;
		virtual void					buildStateList(VecStateList &, VecStateListPos &) const;
		virtual std::string				paramHeader() const = 0;
		virtual std::string				paramReport() const = 0;

		// Member functions related to relative state frequencies
		unsigned						getNStates() const;
		const std::vector<double>	&	getStateFreqs() const;
		virtual void					setNucleotideFreqs(double freqA, double freqC, double freqG, double freqT);
		virtual void					setAllFreqsEqual();
		virtual void					setStateFreqUnnorm(unsigned param_index, double value);
		void							normalizeFreqs();

		// Member functions related to relative rates
		unsigned						getNRatesTotal() const;
		unsigned						getNGammaRates() const;
		void							setNGammaRates(unsigned nRates);
		const std::vector<double>	&	getGammaRateProbs() const;
		void							setAllGammaRateProbsEqual();
		void							recalcRatesAndProbs(std::vector<double> & rates, std::vector<double> & probs) const;
		void							recalcGammaRatesAndBoundaries(std::vector<double> & rates, std::vector<double> & boundaries) const;
		bool							edgeLengthsFixed() const;
		void							fixEdgeLengths();
		void							freeEdgeLengths();
		bool							edgeLenHyperParamFixed() const;
		void							fixEdgeLenHyperprior();
		void							freeEdgeLenHyperprior();
		bool							stateFreqsFixed() const;
		void							fixStateFreqs();
		void							freeStateFreqs();
		void							setFlexRateUpperBound(double new_upper_bound);
		void							setNumFlexSpacers(unsigned s);
		void							setFlexModel();
		void							setNotFlexModel();
		void							fixFlexProbs();
		void							freeFlexProbs();
		void							fixFlexRates();
		void							freeFlexRates();
		void							setFLEXProbParamPrior(ProbDistShPtr d);
		ProbDistShPtr					getFLEXProbParamPrior();
		virtual void					setFlexRateUnnorm(unsigned param_index, double value);
		virtual void					setFlexProbUnnorm(unsigned param_index, double value);
		void							normalizeRatesAndProbs(std::vector<double> & rates, std::vector<double> & probs) const;
		void							setPinvarModel();
		void							setNotPinvarModel();
		bool							pinvarFixed() const;
		void							fixPinvar();
		void							freePinvar();
		double							getPinvar();
		void							setPinvar(double pinv);
		void							setPinvarPrior(ProbDistShPtr d);
		ProbDistShPtr					getPinvarPrior();
		bool							shapeFixed() const;
		void							fixShape();
		void							freeShape();
		double							getShape();
		void							setShape(double alpha);
		void							setPriorOnShapeInverse(bool invert);
		void							setDiscreteGammaShapePrior(ProbDistShPtr d);
		ProbDistShPtr					getDiscreteGammaShapePrior();

		// Member functions related to edge lengths
		void							setEdgeLenHyperPrior(ProbDistShPtr d);
		void							setEdgeLenPrior(ProbDistShPtr d);
		ProbDistShPtr					getEdgeLenHyperPrior();
		ProbDistShPtr					getEdgeLenPrior();

		// Python-specific utilities
#if defined(PYTHON_ONLY)
		boost::python::numeric::array	getPMatrix(double edgeLength) const;
#endif

protected:

	// Constructor
									Model(unsigned numStates);

protected:

	ProbDistShPtr					edgeLenHyperPrior;			/**< The prior distribution governing the mean of the edge length prior if a hierarchical model is used */
	ProbDistShPtr					edgeLenPrior;				/**< The prior distribution governing all edge lengths */
	unsigned						num_states;					/**< The number of states (e.g. 4 for DNA) */
	unsigned						num_gamma_rates;			/**< The number of discrete gamma rate categories. If greater than 1, the model becomes a discrete gamma rate heterogeneity ("G") model. */
	mutable std::vector<double>		gamma_rates_unnorm;			/**< A vector of quantities that yield the relative rates when normalized in recalcRatesAndProbs (length is `num_gamma_rates') */
	mutable std::vector<double>		gamma_rate_probs;			/**< A vector of probabilities that a site falls in any given rate category (length is `num_gamma_rates') */
	std::vector<double>				state_freqs;				/**< A vector of relative state frequencies (length is `num_states') */
	std::vector<double>				state_freq_unnorm;			/**< A vector of quantities that yield the values in `state_freqs' when normalized using the normalizeFreqs member function (length is `num_states') */
	std::vector<std::string>		state_repr;					/**< A vector strings representing the states allowed by this model */
	bool							state_freq_fixed;			/**< If true, the values in `state_freq_params' will not change during MCMC updates */
	mutable MCMCUpdaterVect			freq_params;				/**< Vector of shared pointers to the state frequency parameters (need to retain pointers to these so that the fixed/free status can be changed) */	
	mutable MCMCUpdaterShPtr		edgelen_hyper_param;		/**< Shared pointer to the edge length hyperparameter (need to retain a pointer so that the fixed/free status can be changed) */
	mutable MCMCUpdaterVect			edgelen_params;				/**< Vector of shared pointers to the edge length parameters (need to retain pointers to these so that their fixed/free status can be changed) */
	bool							edge_lengths_fixed;			/**< If true, the value of the edge lengths will not change during MCMC updates */
	bool							edgelen_hyperprior_fixed;	/**< If true, the value of the edge length hyperprior will not change during MCMC updates */
	bool							is_pinvar_model;			/**< If true, a parameter for pinvar will be added to MCMC analysis (pinvar_fixed determines whether it is updated or not) */
	bool							is_flex_model;				/**< If true, the FLEX model of rate heterogeneity will be used instead of the discrete gamma model */
	mutable double					flex_upper_rate_bound;		/**< Largest possible unnormalized relative rate parameter value (lower bound is always 0.0) */
	mutable unsigned				num_flex_spacers;			/**< The number of spacers between rates in the FLEX model rate prior. Spacers act like repelling magnets, keeping adjacent rates from getting too close together. Adding more spacers between each pair of adjacent rates increases the repulsive force. Changing the number of spacers is the only modification allowed to the FLEX model rate prior. */
	bool							flex_probs_fixed;			/**< If true, the values in `flex_prob_params' will not change during MCMC updates */
	bool							flex_rates_fixed;			/**< If true, the values in `flex_rate_params' will not change during MCMC updates */
	mutable MCMCUpdaterVect			flex_rate_params;			/**< Vector of shared pointers to the relative rate parameters used in the FLEX model (need to retain pointers to these so that the fixed/free status can be changed) */	
	mutable MCMCUpdaterVect			flex_prob_params;			/**< Vector of shared pointers to the rate probability parameters used in the FLEX model (need to retain pointers to these so that the fixed/free status can be changed) */	
	ProbDistShPtr					flex_prob_param_prior;		/**< The prior distribution governing all `num_gamma_rate' FLEX probability parameters */
	ProbDistShPtr					flex_rate_param_prior;		/**< The prior distribution governing all `num_gamma_rate' FLEX rate parameters */
	bool							pinvar_fixed;				/**< If true, the value of pinvar will not change during MCMC updates */
	mutable MCMCUpdaterShPtr		pinvar_param;				/**< Shared pointer to the proportion of invariable sites parameter (need to retain a pointer so that the fixed/free status can be changed) */
	ProbDistShPtr					pinvar_prior;				/**< The prior distribution governing the proportion of invariable sites parameter */
	double							pinvar;						/**< The proportion of invariable sites. If non-zero, the model becomes an invariable-sites ("I") model. */
	bool							gamma_shape_fixed;			/**< If true, the value of gamma_shape will not change during MCMC updates */
	mutable MCMCUpdaterShPtr		gamma_shape_param;			/**< Shared pointer to the gamma shape parameter (need to retain a pointer so that the fixed/free status can be changed) */
	ProbDistShPtr					gamma_shape_prior;			/**< The prior distribution governing the discrete gamma shape parameter */
	double							gamma_shape;				/**< Used for discrete gamma rate heterogeneity */
	bool							invert_shape;				/**< If true, gamma_shape_param will hold inverse of shape rather than shape itself */
	CDF								cdf;						/**< Provides cumulative gamma distribution function */
};

typedef boost::shared_ptr<Model> ModelShPtr;

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
		void					calcPMat(double * * pMat, double edgeLength) const;
		virtual std::string		paramHeader() const;
		virtual std::string		paramReport() const;
	};

typedef boost::shared_ptr<JC> JCShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Specialization of the base class Model that represents the Hasegawa-Kishino-Yano (1985) model.
*/
class HKY: public Model
	{
	public:
									HKY();
									~HKY()
										{
										//std::cerr << "HKY dying..." << std::endl;
										}

		virtual std::string			getModelName() const;
		virtual void				createParameters(TreeShPtr t, MCMCUpdaterVect & edgelens, MCMCUpdaterShPtr & edgelen_hyperparam, MCMCUpdaterVect & parameters, bool separate_edgelens) const;
		void						calcPMat(double * * pMat, double edgeLength) const;
		void						fixKappa();
		void						freeKappa();
		double						getKappa();
		void						setKappa(double k);
		void						setKappaPrior(ProbDistShPtr d);
		ProbDistShPtr				getKappaPrior();
		void						setBaseFreqParamPrior(ProbDistShPtr d);
		ProbDistShPtr				getBaseFreqParamPrior();
		void						setKappaFromTRatio(double tratio);
		double						calcTRatio();
		virtual std::string			paramHeader() const;
		virtual std::string			paramReport() const;

protected:
	
	double						kappa;				/**< The transition/transversion rate ratio */
	ProbDistShPtr				kappa_prior;		/**< The prior distribution governing kappa */
	ProbDistShPtr				freq_param_prior;	/**< The prior distribution governing each frequency parameter */
	bool						kappa_fixed;		/**< If true, the value of kappa will not change during MCMC updates */
	mutable MCMCUpdaterShPtr	kappa_param;		/**< Copy of the kappa parameter (saved so that fixed/free status can be changed) */
	};

typedef boost::shared_ptr<HKY> HKYShPtr;

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
		virtual void				createParameters(TreeShPtr t, MCMCUpdaterVect & edgelens, MCMCUpdaterShPtr & edgelen_hyperparam, MCMCUpdaterVect & parameters, bool separate_edgelens) const;
		void						calcPMat(double * * pMat, double edgeLength) const;
		void						fixRelRates();
		void						freeRelRates();
		std::vector<double>			getRelRates();
		void						setRelRates(const std::vector<double> & rates);	
		void						setRelRateUnnorm(unsigned param_index, double value);
		void						setRelRatePrior(ProbDistShPtr d);
		ProbDistShPtr				getRelRatePrior();
		void						setNucleotideFreqs(double freqA, double freqC, double freqG, double freqT);
		void						setAllFreqsEqual();
		void						setStateFreqUnnorm(unsigned param_index, double value);
		void						setBaseFreqParamPrior(ProbDistShPtr d);
		ProbDistShPtr				getBaseFreqParamPrior();
		virtual std::string			paramHeader() const;
		virtual std::string			paramReport() const;
		double						calcTRatio();

	protected:

		std::vector<double>			rel_rates;			/**< A vector containing the six relative rates */
		ProbDistShPtr				rel_rate_prior;		/**< The prior distribution governing each relative rate (usually a gamma distribution with scale 1 and shape equal to the desired Dirichlet parameter) */
		ProbDistShPtr				freq_param_prior;	/**< The prior distribution governing each frequency parameter (usually a gamma distribution with scale 1 and shape equal to the desired Dirichlet parameter) */
		bool						rel_rates_fixed;	/**< If true, the relative rate values will not change during MCMC updates */
		mutable MCMCUpdaterVect		rel_rate_params;	/**< A vector containing copies of all six relative rate parameters (saved so that fixed/free status can be changed) */
		mutable QMatrix				q_matrix;			/**< A QMatrix object used to compute transition probabilities */
	};

typedef boost::shared_ptr<GTR> GTRShPtr;

} // namespace phycas

#include "pyphy/src/likelihood_models.inl"

#endif