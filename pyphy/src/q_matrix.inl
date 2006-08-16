#include <cmath>								// std::sqrt
#include <boost/lambda/lambda.hpp>				// for boost::lambda::_1
#include <boost/lambda/bind.hpp>				// for boost::lambda::bind
#include "pyphy/src/q_matrix.hpp"
#include "pyphy/src/cipres/AllocateMatrix.hpp"
#include <functional>
#include <algorithm>
#include <boost/format.hpp>

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|
*/
inline QMatrix::QMatrix()
  : dimension(0), flat_length(0) 
  , qmat(0), qmat_begin(0), qmat_end(0) 
  , w(0), w_begin(0), w_end(0)
  , z(0), z_begin(0), z_end(0)
  , fv(0) 
  , q_dirty(true)
	{
	// Set up Q matrix representing the JC69 model by default
	rr.assign(6, 1.0);
	pi.assign(4, 0.25);
	sqrtPi.assign(4, 0.5);
	recalcQMatrix();
	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
inline QMatrix::~QMatrix()
	{
	clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns object to its just-constructed state, with the exception that the vectors `rr', `pi' and `sqrtPi' are not 
|	cleared.
*/
inline void QMatrix::clearAllExceptFreqsAndRates()
	{
	dimension	= 0;
	flat_length	= 0; 
	qmat_begin	= 0; 
	qmat_end	= 0; 
	w_begin		= 0; 
	w_end		= 0; 
	z_begin		= 0; 
	z_end		= 0; 

	dim_vect.clear();

	if (qmat)
		{
		DeleteTwoDArray<double>(qmat);
		qmat = 0;
		}
	
	if (w)
		{
		delete [] w;
		w = 0;
		}

	if (fv)
		{
		delete [] fv;
		fv = 0;
		}

	if (z)
		{
		DeleteTwoDArray<double>(z);
		z = 0;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns object to its just-constructed state by calling clearAllExceptFreqsAndRates and then also clearing the 
|	`rr', `pi' and `sqrtPi' vectors.
*/
inline void QMatrix::clear()
	{
	clearAllExceptFreqsAndRates();
	pi.clear();
	sqrtPi.clear();
	rr.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the current value of the data member `dimension', which is the number of rows (and
|	columns) in the Q matrix.
*/
inline unsigned QMatrix::getDimension()
	{
	return dimension;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Copies elements from the supplied vector `rates' to the data member `rr' and sets `q_dirty' to true.
*/
inline void QMatrix::setRelativeRates(const std::vector<double> & rates)
	{
	rr.resize(rates.size());
	std::copy(rates.begin(), rates.end(), rr.begin());
	q_dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Copies elements from the supplied vector `rates' to the data member `rr' and sets `q_dirty' to true.
*/
inline void QMatrix::setStateFreqs(const std::vector<double> & freqs)
	{
	pi.resize(freqs.size());
	std::copy(freqs.begin(), freqs.end(), pi.begin());

	sqrtPi.resize(freqs.size());
	std::transform(pi.begin(), pi.end(), sqrtPi.begin(), boost::lambda::bind(static_cast<double(*)(double)>(&std::sqrt), boost::lambda::_1));

	q_dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reallocates all arrays and vectors (except `rr', `pi' and `sqrtPi') according to the value of `new_dim' and sets 
|	`dimension' equal to `new_dim' before returning. Assumes `new_dim' is greater than zero. Note that `qmat', `w', `z'
|	and `fv' are allocated but not filled with any values when this function returns.
*/
inline void QMatrix::redimension(unsigned new_dim)
	{
	assert(new_dim > 0);
	clearAllExceptFreqsAndRates();

	dimension = new_dim;
	flat_length = new_dim*new_dim;

	// Set the shape of a NumArray object that represents Q, P (transition probs), E (eigenvectors) or V (eigenvalues)
	dim_vect.push_back((int)dimension);
	dim_vect.push_back((int)dimension);

	// Create and fill qmat with default (Mk model) values 
	qmat = NewTwoDArray<double>(dimension, dimension);
	assert(qmat);

	qmat_begin		= &qmat[0][0];
	qmat_end		= qmat_begin + flat_length;

	// Create w (array of eigenvalues)
	w = new double[dimension];
	w_begin = &w[0];
	w_end = w_begin + dimension;

	// Create z (two-dimensional array of eigenvectors)
	z = NewTwoDArray<double>(dimension, dimension);
	z_begin		= &z[0][0];
	z_end		= z_begin + flat_length;

	// Create fv (workspace)
	fv = new double[dimension];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `q_dirty' is true, calls recalcQMatrixImpl to recompute the `qmat', `z' and `w' data members. If `q_dirty' is 
|	false, however, this function returns immediately and thus (because it is declared inline) is computationally 
|	inexpensive.
*/
inline void QMatrix::recalcQMatrix()
	{
	if (!q_dirty)
		return;
	recalcQMatrixImpl();
	}

#if defined(PYTHON_ONLY)
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the entries in the data member `qmat' as a NumArray. Calls recalcQMatrix first to ensure that `qmat' is up
|	to date.
*/
inline boost::python::numeric::array QMatrix::getQMatrix()
	{
	recalcQMatrix();
	return num_util::makeNum(qmat_begin, dim_vect);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the eigenvectors in the data member `z' as a NumArray. Calls recalcQMatrix first to ensure that `z' is up
|	to date.
*/
inline boost::python::numeric::array QMatrix::getEigenVectors()
	{
	recalcQMatrix();
	return num_util::makeNum(z_begin, dim_vect);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the eigenvectors in the data member `z' as a NumArray. Calls recalcQMatrix first to ensure that `z' is up
|	to date.
*/
inline boost::python::numeric::array QMatrix::getPMatrix(double edgelen)
	{
	std::vector<double> p;
	recalcPMatrix(p, edgelen);
	return num_util::makeNum(&p[0], dim_vect);	//PELIGROSO
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the eigenvalues in the data member `w' as a vector. Calls recalcQMatrix first to ensure that `w' is up to 
|	date.
*/
inline std::vector<double> QMatrix::getEigenValues()
	{
	recalcQMatrix();
	std::vector<double> v(dimension, 0.0);
	std::copy(w_begin, w_end, v.begin());
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
inline std::string QMatrix::showQMatrix()
	{
	recalcQMatrix();
	return showMatrixImpl(qmat_begin);
	}

} // namespace phycas