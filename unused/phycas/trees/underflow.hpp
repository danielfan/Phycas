#ifndef UNDERFLOW_HPP
#define UNDERFLOW_HPP

#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>

/*----------------------------------------------------------------------------------------------------------------------
|	Manages underflow protection for site likelihoods, maintaining an array of unsigned values representing the number
|	of times each site has had to be rescaled during the last computation of the log-likelihood.
*/
class UnderflowManager
	{
	public:
								UnderflowManager();
								~UnderflowManager();

		// Accessors
		//
		double					GetLikeCutoff() const;
		double					GetLnLikeCutoff() const;
		unsigned				*GetBouncesArray();

		// Modifiers
		//
		void					SetLikeCutoff(double c);
		void					ProtectSite(unsigned pattern);
		void					BounceSite(unsigned pattern, unsigned num_times);
		void					InitBouncesVector(unsigned npatterns);

		// Queries
		//
		bool					HasBeenCorrected(unsigned i);
		bool					HasUnderflowed(unsigned i);
		bool					AnySiteHasUnderflowed();

		// Utilities
		//
		double					GetCorrectedLnLike(double x, unsigned i);
		unsigned				CalcNumCorrectionsNeeded(double x);
		double					CalcCorrectionFactor(unsigned n);
		void					Reset();

		// Debugging
		//
		double					DebugGetAvgNumBounces();
		std::string				DebugCreateBounceVectorRepresentation();

	protected:
		bool					any_site_has_underflowed;	/**< False initially, but becomes true whenever ProtectSite is called */
		double					kLikeCutoff;	/**< Site likelihood is rescaled when it becomes smaller than this value */
		double					kLnLikeCutoff;	/**< Used to rescale the final log-likelihood of the tree based on the number of times kLikeCutoff was breached */
		double					kSquaredCutoff;	/**< Used in CalcNumCorrectionsNeeded() */
		std::vector<unsigned>	vBounces;		/**< Element i is the number of times site i needed to be rescaled (bounced back up) */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Sets default value for `kLikeCutoff' and sets `any_site_has_underflowed` to its initial value of false.
*/
inline UnderflowManager::UnderflowManager()
	{
#if 0 //defined(POL_TEMP)
	SetLikeCutoff(0.00001);
#else
	SetLikeCutoff(1.0e-100);
#endif
	Reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the Reset() member function.
*/
inline UnderflowManager::~UnderflowManager()
	{
	Reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the data member `kLikeCutoff'.
*/
inline double UnderflowManager::GetLikeCutoff() const
	{
	return kLikeCutoff;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the data member `kLnLikeCutoff'.
*/
inline double UnderflowManager::GetLnLikeCutoff() const
	{
	return kLnLikeCutoff;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns pointer to first element of the `vBounces' vector, allowing efficient traversal of the underlying array.
*/
inline unsigned *UnderflowManager::GetBouncesArray()
	{
	return &vBounces[0];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces current value of `kLikeCutoff' with supplied value `c' and recomputes `kLnLikeCutoff'. Asserts `c' > 0.0.
*/
inline void UnderflowManager::SetLikeCutoff(
  double c)	/**< is the new value of `kLikeCutoff' */
	{
	assert(c > 0.0);
	kLikeCutoff = c;
	kLnLikeCutoff = std::log(c);
	kSquaredCutoff = c*c;
	assert(kSquaredCutoff > 0.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If the length of the `vBounces' vector is equal to `npatterns', this function resets all elements of `vBounces' that
|	are not UINT_MAX to 0, otherwise it resizes the `vBounces' vector to have length 'npatterns' and sets every element
|	to UINT_MAX. Thus, to force a reinitialization of `vBounces' (e.g. if the model or data set changes), call the 
|	Reset() function, which will clear the `vBounces' vector.
*/
inline void UnderflowManager::InitBouncesVector(
  unsigned npatterns)	/**< is the length of the new `vBounces' vector */
	{
	if (npatterns == (unsigned)vBounces.size())
		{
		// Initialize the sites that are being protected from underflow to 0
		//
		std::vector<unsigned>::iterator it = vBounces.begin();
		for (; it != vBounces.end(); it++)
			{
			if ((*it) < UINT_MAX)
				(*it) = 0;
			}
		}
	else
		{
		// Initialize all sites to UINT_MAX to indicate that they have not proven a 
		// need for underflow protection
		//
		vBounces.resize(npatterns);
		std::fill_n(vBounces.begin(), npatterns, UINT_MAX);
		}

#if 0  //defined(POL_TEMP)
		// for testing purposes, pretend like every third site needs protecting
		unsigned i = 0;
		unsigned count = 0;
		std::vector<unsigned>::iterator it = vBounces.begin();
		for (; it != vBounces.end(); it++)
			{
			if (i % 3 == 0)
				{
				(*it) = 0;
				count++;
				}
			i++;
			}
		//std::cout << "*** pretending that " << count << " sites had underflow problems" << std::endl;
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Bounces site `i', which simply means that the `i'th element of the `vBounces' vector is incremented by 1.
|	Asserts that `i' is a valid index into `vBounces'.
*/
inline void UnderflowManager::BounceSite(
  unsigned pattern,		/**< is the 0-based index of the pattern to bounce */
  unsigned num_times)	/**< is the number of times to bounce it */
	{
	assert(pattern >= 0);
	assert(pattern < vBounces.size());
	assert(vBounces[pattern] < UINT_MAX);	// shouldn't be trying to bounce a site that had no underflow problems
	vBounces[pattern] += num_times;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if pattern `i' is known to underflow and has been corrected at least once. Specifically, returns true 
|	if `vBounces'[i] is not UINT_MAX (indicating the site has underflowed in the past) and is greater than 0 (indicating
|	that at least one correction has been applied).
*/
inline bool UnderflowManager::HasBeenCorrected(
  unsigned i)	/**< is the 0-based index of site to check */
	{
	assert(i >= 0);
	assert(i < vBounces.size());
	return (vBounces[i] < UINT_MAX && vBounces[i] > 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if pattern `i' is known to underflow (i.e. it underflowed the last time CalcLogLikelihood was called).
|	Specifically, returns true if `vBounces'[i] is not UINT_MAX.
*/
inline bool UnderflowManager::HasUnderflowed(
  unsigned i)	/**< is the 0-based index of site to check */
	{
	assert(i >= 0);
	assert(i < vBounces.size());
	return (vBounces[i] < UINT_MAX);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `any_site_has_underflowed'. Use to avoid unnecessary computations for data sets that never have 
|	underflow problems.
*/
inline bool UnderflowManager::AnySiteHasUnderflowed()
	{
	return any_site_has_underflowed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Corrects the site log-likelihood by dividing out the factor (or factors) used to keep it from underflowing while 
|	computing the conditional likelihood arrays. The number of times the conditional likelihoods for site `i' were
|	bounced is stored in the vector `vBounces', and each time a bounce was necessary, the conditinal likelihoods were
|	divided by the amount kLikeCutoff. This means at the end, an amount `vBounces'[i]*kLnLikeCutoff must be added to the
|	log of the site likelihood to return the site log-likelihood to the correct value. Asserts `x' > 0.0 and that `i'
|	is valid index into `vBounces'.
*/
inline double UnderflowManager::GetCorrectedLnLike(
  double x,		/**< the uncorrected site likelihood */
  unsigned i)	/**< is the 0-based index of the site pattern */
	{
	assert(x > 0.0);
	double lnx = std::log(x);
	unsigned n_times_bounced = vBounces[i];
	assert(n_times_bounced < UINT_MAX);
	double lnf = kLnLikeCutoff*(double)n_times_bounced;

#	if 0 // defined(POL_TEMP)
	//std::ofstream f("doof.txt");
	std::ostream &f = std::cerr;
	f << "\nGetCorrectedLnLike(" << x << ", " << i << ")\n";
	f << "\n  lnx = " << lnx;
	f << "\n  lnf = " << lnf;
	f << "\n  n_times_bounced = " << n_times_bounced;
	f << std::endl;
	//f.close();
#	endif

	return lnx + lnf;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value z such that x/(kLikeCutoff^z) >= kLikeCutoff. During computation of conditional likelihood arrays,
|	if a value is less than kLikeCutoff then it needs to be multiplied by some factor (1/kLikeCutoff) enough times that
|	it is no longer less than kLikeCutoff. This function computed the number of multiplications necessary.
*/
inline unsigned UnderflowManager::CalcNumCorrectionsNeeded(
  double x)	/**< is the value needing to be corrected */
	{
	assert(x > 0.0);
	unsigned z = 1 + (unsigned)(kSquaredCutoff/x);
	return z;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns (1/kLikeCutoff) raised to the power `n', where `n' is normally computed using the CalcNumCorrectionsNeeded
|	function.
*/
inline double UnderflowManager::CalcCorrectionFactor(
  unsigned n)	/**< is the number of times value needs to be corrected */
	{
	assert(n > 0);
	double f = -(double)n*kLnLikeCutoff;
	return exp(f);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Changes value stored in `vBounces' vector for site whose index is given by `pattern' to 0. Asserts that the original
|	value was UINT_MAX, which is the value used to initialize every element of `vBounces'. Setting an element to 0 
|	indicates that this site is to be protected from underflow in the future. Also sets `any_site_has_underflowed' to 
|	true to indicate that at least one site has had underflow problems.
*/
inline void UnderflowManager::ProtectSite(
  unsigned pattern)	/**< is the number of times value needs to be corrected */
	{
	assert(pattern >= 0);
	//@mth: paul sz was unused 
	// unsigned sz = (unsigned)vBounces.size();
	assert(pattern < vBounces.size());
#if !defined(POL_TEMP)
	assert(vBounces[pattern] == UINT_MAX);
#endif
	vBounces[pattern] = 0;
	any_site_has_underflowed = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears the `vBounces' vector, which forces a complete reinitialization when InitBouncesVector() is next called. Also
|	sets `any_site_has_underflowed' to its initial value of false.
*/
inline void UnderflowManager::Reset()
	{
	any_site_has_underflowed = false;
	vBounces.clear();
	}

#endif

