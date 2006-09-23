#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#define NO_IMPORT_ARRAY

#include "probability_distribution.hpp"
#include "ncl/nxs_exception.hpp"

namespace phycas
{

DirichletDistribution::DirichletDistribution(const std::vector<double> &params) 
	{
	unsigned params_size = (unsigned)params.size();
	//std::cerr << "params_size = " << params_size << std::endl;
	PHYCAS_ASSERT(params_size > 1);

	// Reserve necessary space to save copying
	//
	dirParams.reserve(params_size);
	scratchSpace.reserve(params_size);
	paramDistributions.reserve(params_size);

#	if defined NCL_NXS_THROW_UNDEFINED
		if (params.size() < 2)
			throw NxsX_UndefinedException("Illegal Dirichlet", __FILE__, __LINE__);
#	endif

	std::vector<double>::const_iterator iter;
	for (iter = params.begin(); iter != params.end(); ++iter)
		{
		double c_i = *iter;
		dirParams.push_back(c_i);
		scratchSpace.push_back(0.0);
		GammaDistribution gamma_i = GammaDistribution(c_i, 1.0);
		gamma_i.SetLot(&myLot);	// ensure that they all use same random number generator
		paramDistributions.push_back(gamma_i);
		}
	}

#if 0 && defined(PYTHON_ONLY)
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the variance-covariance matrix in the form of a numarray object.
*/
boost::python::numeric::array DirichletDistribution::GetVarCovarMatrix()
	{
	unsigned i, j;
	unsigned dim = (unsigned)dirParams.size();

	std::vector<int> dims;
	dims.push_back((int)dim);
	dims.push_back((int)dim);

	double c = 0.0;
	for (i = 0; i < dim; ++i)
		{
		c += dirParams[i];
		}
	double denom = c*c*(c + 1.0);

	VecDbl V;
	for (i = 0; i < dim; ++i)
		{
		for (j = 0; j < dim; ++j)
			{
			if (i == j)
				{
				double var_i = dirParams[i]*(c - dirParams[i])/denom;
				V.push_back(var_i);
				}
			else
				{
				double cov_ij = -dirParams[i]*dirParams[j]/denom;
				V.push_back(cov_ij);
				}
			}
		}

	return num_util::makeNum<double>(&V[0], dims);
	}
#endif

} // namespace phycas
