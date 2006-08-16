#if ! defined(TOPO_PRIOR_CALCULATOR_INL)
#define TOPO_PRIOR_CALCULATOR_INL

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `is_rooted' to true, `is_resolution_class_prior' to true, `C' to 1.0 and `ntax' to 4. Then calls
|	Reset(), which calculates the `counts' and `topology_prior' vectors to reflect the resolution class topology priors
|	for unrooted trees with 4 taxa. 
*/
inline TopoPriorCalculator::TopoPriorCalculator() 
	{
	topo_priors_dirty			= true;
	is_rooted					= false;
	is_resolution_class_prior	= true;
	C							= 1.0;
	ntax						= 4;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor clears the vectors `counts' and `topology_prior'.
*/
inline TopoPriorCalculator::~TopoPriorCalculator() 
	{
	counts.clear();
	topology_prior.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `new_ntax' differs from `ntax', sets `ntax' to `new_ntax' and sets `topo_priors_dirty'. Returns immediately 
|	without taking any action if `ntax' equals `new_ntax'. Assumes that `new_ntax' is greater than 1 if `is_rooted' is 
|	true, or that `new_ntax' is greater than 2 if `is_rooted' is false.
*/
inline void TopoPriorCalculator::SetNTax(
  unsigned new_ntax)	/**< is the new number of taxa */
	{
	if (ntax != new_ntax)
		{
		// Set ntax to the new value
		assert(new_ntax > (unsigned)(is_rooted ? 1 : 2));
		ntax = new_ntax;

		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Forces recalculation of `counts' and `polytomy_prior' vectors using current values of `ntax', `is_rooted' and 
|	`is_resolution_class_prior'.
*/
inline void TopoPriorCalculator::Reset()
	{
	unsigned num_internal_nodes = (is_rooted ? (ntax - 1) : (ntax - 2));
	RecalcCountsAndPriorsImpl(num_internal_nodes);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `ntax'.
*/
inline unsigned TopoPriorCalculator::GetNTax() const
	{
	return ntax;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `is_rooted' data member to true. There are more rooted than unrooted trees for the same value of `ntax', so
|	this setting is important when asking questions that require knowledge of the numbers of possible trees.
*/
inline void TopoPriorCalculator::ChooseRooted()
	{
	if (!is_rooted)
		{
		is_rooted = true;
		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `is_rooted' data member to false. There are more rooted than unrooted trees for the same value of `ntax', so
|	this setting is important when asking questions that require knowledge of the numbers of possible trees.
*/
inline void TopoPriorCalculator::ChooseUnrooted()
	{
	if (is_rooted)
		{
		is_rooted = false;
		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns copy of the `counts' vector, which contains in its first element the total number of tree topologies having
|	any number of internal nodes, and which contains in its mth element (m > 0) the number of tree topologies having 
|	exactly m internal nodes.
*/
inline std::vector<double> TopoPriorCalculator::GetCountsVect()
	{
	//@POL this function could be const were it not for lazy evaluation
	if (topo_priors_dirty)
		Reset();
	return counts;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `is_resolution_class_prior' data member to true.
*/
inline void TopoPriorCalculator::ChooseResolutionClassPrior()
	{
	if (!is_resolution_class_prior)
		{
		is_resolution_class_prior = true;
		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `is_resolution_class_prior' data member to false.
*/
inline void TopoPriorCalculator::ChoosePolytomyPrior()
	{
	if (is_resolution_class_prior)
		{
		is_resolution_class_prior = false;
		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `C' data member to `c'. Assumes `c' is greater than 0.0.
*/
inline void TopoPriorCalculator::SetC(double c)
	{
	assert(c > 0.0);
	if (c != C)
		{
		C = c;
		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `C' data member.
*/
inline double TopoPriorCalculator::GetC() const
	{
	return C;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns copy of the `topology_prior' vector, which contains in its mth element the unnormalized prior for tree 
|	topologies having exactly m internal nodes. The 0th element of `topology_prior' holds the normalizing constant.
*/
inline std::vector<double> TopoPriorCalculator::GetTopoPriorVect()
	{
	//@POL this function could be const were it not for lazy evaluation
	if (topo_priors_dirty)
		Reset();
	return topology_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the unnormalized topology prior. This represents the resolution class prior if
|	`is_resolution_class_prior' is true, otherwise it represents the polytomy prior. Assumes `m' is less than the 
|	length of the `topology_prior' vector.
*/
inline double TopoPriorCalculator::GetLnTopologyPrior(
  unsigned m)	/**< is the number of internal nodes */
	{
	//@POL currently using lazy evaluation to avoid recalculating counts and topology_prior vectors more than necessary, 
	// but this may actually slow things down in practice because functions like this one that are called often
	// have to always check topo_priors_dirty to see if they need to recalculate the two vectors. Might try experimenting
	// with greedy evaluation to see if that speeds up MCMC noticeably

	if (topo_priors_dirty)
		Reset();
	assert(m < topology_prior.size());
	return topology_prior[m];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the normalized topology prior. This represents the resolution class prior if
|	`is_resolution_class_prior' is true, otherwise it represents the polytomy prior. Assumes `m' is less than the length
|	of the `topology_prior' vector. The log of the normalized topology prior is obtained as `topology_prior'[m] minus
|	`topology_prior'[0] (the 0th element of `topology_prior' holds the log of the normalization constant).
*/
inline double TopoPriorCalculator::GetLnNormalizedTopologyPrior(
  unsigned m)	/**< is the number of internal nodes */
	{
	if (topo_priors_dirty)
		Reset();
	assert(m < topology_prior.size());
	return (topology_prior[m] - topology_prior[0]);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the normalizing constant for the topology prior. This value is stored in
|	`topology_prior'[0].
*/
inline double TopoPriorCalculator::GetLnNormConstant()
	{
	if (topo_priors_dirty)
		Reset();
	return topology_prior[0];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `is_resolution_class_prior'.
*/
inline bool TopoPriorCalculator::IsResolutionClassPrior() const
	{
	return is_resolution_class_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `is_resolution_class_prior' is true, returns false; if `is_resolution_class_prior' is false, returns true.
*/
inline bool TopoPriorCalculator::IsPolytomyPrior() const
	{
	return !is_resolution_class_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `is_rooted'.
*/
inline bool TopoPriorCalculator::IsRooted() const
	{
	return is_rooted;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `is_rooted' is true, returns false; if `is_rooted' is false, returns true.
*/
inline bool TopoPriorCalculator::IsUnrooted() const
	{
	return !is_rooted;
	}

} // namespace phycas

#endif
