#include "pyphy/src/q_matrix.hpp"
#include "pyphy/src/xlikelihood.hpp"
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Uses `pi' and `rr' vectors to create the two-dimensional array `qmat'. The data member `qmat' is reallocated if 
|	changes in `pi' and `rr' imply a new dimension (or if `qmat' is NULL). Recomputes eigenvalues and eigenvectors 
|	corresponding to the new Q matrix. Assumes `q_dirty' is true, and sets `q_dirty' to false when finished. Throws an 
|	exception if lengths of the `pi' and `rr' vectors are incompatible.
*/
void QMatrix::recalcQMatrixImpl()
	{
	assert(q_dirty);

	// pi, sqrtPi and rr should all be non-empty
	assert(pi.size() > 0);
	assert(sqrtPi.size() == pi.size());
	assert(rr.size() > 0);

	// First check to make sure lengths of `rr' and `pi' are compatible with each other
	unsigned dim_pi = (unsigned)pi.size();
	unsigned dim_rr = (unsigned)((1.0 + std::sqrt(1.0 + 8.0*rr.size()))/2.0); // rr.size() = (n^2 - n)/2, so use quadratic formula to get n
	if (dim_pi != dim_rr)
		{
		throw XLikelihood("Number of relative rates and number of state frequencies specified are incompatible");
		}

	// If qmat is NULL or dimension differs from dim_pi and dim_rr, reallocate qmat, z, w and dv
	if (!qmat || (dimension != dim_pi))
		redimension(dim_pi);

	// This vector will hold the row sums
	std::vector<double> row_sum(dimension, 0.0);

	unsigned i, j, k = 0;
	double sum_for_scaling = 0.0;
	for (i = 0; i < dimension; ++i)
		{
		double pi_i = pi[i];
		double sqrtPi_i = sqrtPi[i];
		for (j = i + 1; j < dimension; ++j, ++k)
			{
			double rr_k = rr[k];
			double pi_j = pi[j];
			double sqrtPi_j = sqrtPi[j];

			// set value in upper triangle
			qmat[i][j] = rr_k*sqrtPi_i*sqrtPi_j;

			// set value in lower triangle
			qmat[j][i] = qmat[i][j];

			// add to relevant row sums
			row_sum[i] += rr_k*pi_j;
			row_sum[j] += rr_k*pi_i;

			// add to total expected number of substitutions 
			sum_for_scaling += 2.0*pi_i*rr_k*pi_j;
			}
		qmat[i][i] = -row_sum[i];
		}

	assert(sum_for_scaling > 0.0);
	edgelen_scaler = 1.0/sum_for_scaling;

	// Calculate eigenvalues and eigenvectors
	int err_code = EigenRealSymmetric(dimension, qmat, w, z, fv);
	if (err_code != 0)
		{
		clearAllExceptFreqsAndRates();
		throw XLikelihood("Error in the calculation of eigenvectors and eigenvalues of the Q matrix");
		}

	q_dirty = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
void QMatrix::recalcPMat(
  double * * pmat,		/**< */
  double edgelen) 		/**< */
	{
	assert(edgelen >= 0.0);
	recalcQMatrix();

	// Adjust the supplied edgelen to account for the fact that the expected number of substitutions
	// implied by the Q matrix is not unity
	double v = edgelen*edgelen_scaler;

	// Exponentiate eigenvalues and put everything back together again
	for (unsigned i = 0; i < dimension; ++i)
		{
		double sqrtPi_i = sqrtPi[i];
		for (unsigned j = 0; j < dimension; ++j)
			{
			double factor = sqrtPi[j]/sqrtPi_i;
			double Pij = 0.0;
			for (unsigned k = 0; k < dimension; ++k)
				{
				double tmp = z[i][k]*z[j][k]*std::exp(w[k]*v);
				Pij +=  tmp;
				}
			pmat[i][j] = Pij*factor;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
void QMatrix::recalcPMatrix(
  std::vector<double> & P,	/**< */
  double edgelen)			/**< */
	{
	assert(edgelen >= 0.0);
	recalcQMatrix();

	// Adjust the supplied edgelen to account for the fact that the expected number of substitutions
	// implied by the Q matrix is not unity
	double v = edgelen*edgelen_scaler;

	P.clear();
	P.reserve(flat_length);

	// Exponentiate eigenvalues and put everything back together again
	for (unsigned i = 0; i < dimension; ++i)
		{
		double sqrtPi_i = sqrtPi[i];
		for (unsigned j = 0; j < dimension; ++j)
			{
			double factor = sqrtPi[j]/sqrtPi_i;
			double Pij = 0.0;
			for (unsigned k = 0; k < dimension; ++k)
				{
				double tmp = z[i][k]*z[j][k]*std::exp(w[k]*v);
				Pij +=  tmp;
				}
			P.push_back(Pij*factor);
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
inline std::string QMatrix::showMatrixImpl(const double * q) const
	{
	unsigned i, j;

	// Output one column for row labels and a label for every column of qmat
	std::string s = str(boost::format("%12s ") % " ");
	for (i = 0; i < dimension; ++i)
		{
		s += str(boost::format("%12d ") % i);
		}
	s += "\n";

	// Output rows of q
	double * pq = (double *)q;	//PELIGROSO
	for (i = 0; i < dimension; ++i)
		{
		s += str(boost::format("%12d ") % i);
		for (j = 0; j < dimension; ++j)
			{
			s += str(boost::format("%12.5f ") % *pq++);
			}
		s += "\n";
		}

	return s;
	}
