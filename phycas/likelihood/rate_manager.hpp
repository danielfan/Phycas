#ifndef PHO_RATE_MANAGER_H
#define PHO_RATE_MANAGER_H

#include <boost/shared_ptr.hpp>

class GibbsParameter;
typedef GibbsParameter PhoParameter;
typedef boost::shared_ptr<GibbsParameter> ConstParamShPtr;
#include "phycas/modules/mcmc/gibbs_param.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|	This class is used as a field of model so that any implemented model can be used with any combination of three 
|	different types of rate heterogeneity among sites: invariant sites, discrete approximation to a gamma-
|	distributed rates across sites, and differing rates for different subsets of the data.
|
|	Essentially the job of the class is to maintain an array of doubles, `rates', which holds the rate multiplier
| 	for every variable rate category. Model's CalculatePrMat(double ***,double) creates an array of Prob. of change 
|	matrices for each of the i rate categories by multiplying the value in PhoRateManager->rate[i] times the branch 
|	length.
|
|	The mean rate across all sites will be 1 (and the Parameter * meanRate will be NULL) unless subset specific rates 
|	are being used (in which case the value of meanRate will hold the mean rate for this model (each subset will have 
|	its own model).
|
|	In the case of the invariant sites rate heterogeneity, the rate manager simply makes sure that the value in rates 
|	times the proportion of VARYING sites equals the mean rate for the model.
|
|	For gamma distributed rates the distribution of rates among the categories is drawn from a gamma with parameter
|	gammaParam (either the shape parameter or the variance depending on conditional compilation). The rates are then 
|	scaled so the mean rate = the value in meanRate.
*/
class PhoRateManager : public boost::noncopyable
	{
	public:
	
								PhoRateManager(ConstParamShPtr pinv, ConstParamShPtr gamParam, ConstParamShPtr meanRateParam, unsigned nGamCats);
		
		void					SetNumVariableRateCategories(unsigned nr);
		unsigned				GetNumVariableRateCategories() const;

		unsigned				GetNumParameters() const;
		const PhoParameter		*GetParameter(unsigned ind) const;
		bool					OwnsParameter(const PhoParameter *p) const;
		bool					UsingPInvar() const;
		double					GetPInvar() const;
		std::vector<double>		GetRepRates() const;

		class					XIterMax{}; /**< Possibly thrown by a function that causes gamma rates recalculation with extreme parameters */

	protected:
	
								~PhoRateManager();  // make this virtual if we delete other than via derived class (e.g. PhoRateHetModel)
		
		void					CalculateRates();
		
		ConstParamShPtr const	pInvar;			/**< Proportion of invariant sites parameter. */	
		ConstParamShPtr const	gammaParam;		/**< Gamma shape parameter. */
		ConstParamShPtr const	meanRate;		/**< Mean relative rate. */

		double					*rates;			/**< Array of relative rates. */
		unsigned				nVarRateCats;	/**< Number of variable rate categories. */
		
	private:

		STATELESS_FUNC double			point_chi2 (double pr, double v);
		STATELESS_FUNC double			gammq(double a, double x);
		STATELESS_FUNC void				gser(double& gamser, double a, double x, double& gln);
		STATELESS_FUNC void				gcf(double& gammcf, double a, double x, double& gln);
		STATELESS_FUNC double			gammln(double xx);
		STATELESS_FUNC double			point_normal(double pr);

		void 					FlushRates();
		void					MultiplyRatesBy(double);
		void					ReallocateRates();
			
		double 					pinvBrCorrection; /**< number that must be multiplied by the branch length so that the exp. # changes = brlen */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of rate means or medians.
*/
inline std::vector<double> PhoRateManager::GetRepRates() const
	{
	return std::vector<double>(rates, rates + nVarRateCats);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls FlushRates().
*/
inline PhoRateManager::~PhoRateManager()
	{
	FlushRates();	
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the numerical value of the `pInvar' parameter.
*/
inline double PhoRateManager::GetPInvar() const 
	{
	return (pInvar == NULL ? 0.0 : pInvar->GetValue());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if this object owns the parameter pointed to by `p'. \todo Need to make this clearer.
*/
inline bool PhoRateManager::OwnsParameter(const PhoParameter *p) const
	{
	return (p != NULL && (p == pInvar.get() || p == gammaParam.get()|| p == meanRate.get()));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the `nVarRateCats' data member.
*/
inline unsigned PhoRateManager::GetNumVariableRateCategories() const
	{
	return nVarRateCats;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of the data member `nVarRateCats' and then reallocates the `rates' array and calls CalculateRates()
|	to recompute the new representative rate values.
*/
inline void PhoRateManager::SetNumVariableRateCategories(unsigned nr)
	{
	assert(nr > 0);
	assert(nr < 100);	// the 100 is arbitrary, but some upper limit is needed to catch grievous errors
	nVarRateCats = nr;
	FlushRates();
	ReallocateRates();
	CalculateRates();
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Returns the parameter specified by the index `ind'. The ordering is shape, pinv, and then mean, but if one of the 
|	parameters is not being used it is removed from the list (so 0 could be mean rate if rates = equal and pinvar is 
|	set to 0).
*/
inline const PhoParameter *PhoRateManager::GetParameter(
  unsigned ind) const	/**< is the index */
	{
	if (gammaParam && nVarRateCats > 1)
		{
		if (ind == 0)
			return gammaParam.get();
		--ind;
		}
	if (UsingPInvar())
		{
		if (ind == 0)
			return pInvar.get();
		--ind;
		}
	if (ind == 0)
		return meanRate.get();
	return NULL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if invariant sites model is in effect.
*/
inline bool PhoRateManager::UsingPInvar() const
	{
	return (pInvar);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Multiplies each rate by the same amount. Called when `pInvar' changes (with the ratio of new to old pinvBrCorrection 
|	as the argument) and at the end of calculate rates (with pinvBrCorrection as the argument).
*/
inline void PhoRateManager::MultiplyRatesBy(
  double m)	/**< the number to be multiplied to every rate */
	{
	for (unsigned i = 0; i < nVarRateCats; ++i)
		rates[i] *= m;
	}

#endif
