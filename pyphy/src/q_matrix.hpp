#if ! defined(Q_MATRIX_HPP)
#define Q_MATRIX_HPP

extern "C"
{
#include "linalg.h"
}

#if defined(_MSC_VER)
#	pragma warning(disable: 4267)	// warning about loss of data when converting size_t to int
#endif

#include <cmath>
#include <vector>
#if defined(PYTHON_ONLY)
#	include <boost/python/tuple.hpp>
#	include <boost/python/numeric.hpp>
#	include "pyphy/src/thirdparty/num_util.h"
#endif

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	The QMatrix class manages a matrix of relative instantaneous rates. The primary purpose of the class is to compute
|	transition probabilities for the GTR model (and other models in which there do not exist closed-form expressions for
|	transition probabilities). The Q matrix is a square stochastic (i.e. rows sum to 0.0) matrix that need not be 
|	symmetric itself, but has the property that the matrix product Pi*Q is symmetric (where Pi is a diagonal matrix of 
|	state frequencies). When setQMatrix is used to set the Q matrix (or setQMatrixElement is used to change one element
|	of the Q matrix), the QMatrix class computes the eigenvalue and eigenvectors of Q. The eigenvalues and eigenvectors
|	are used by calcPMatrix to compute transition matrices needed for likelihood calculations.
*/
class QMatrix
	{
	friend class GTR;

	public:
		
										QMatrix();
										~QMatrix();
		
		// Accessors
		//
		unsigned						getDimension();

		// Modifiers
		//
		void							setRelativeRates(const std::vector<double> & rates);
		void							setStateFreqs(const std::vector<double> & freqs);

		// Utilities
		//
#if defined(PYTHON_ONLY)
		boost::python::numeric::array	getQMatrix();
		boost::python::numeric::array	getEigenVectors();
		boost::python::numeric::array	getPMatrix(double edgelen);
#endif
		std::vector<double>				getEigenValues();

	private:

		void							recalcPMat(double * * pmat, double edgelen);	// used by GTR
		void							recalcPMatrix(std::vector<double> & P, double edgelen);
		std::string						showQMatrix();
		void							clear();
		void							recalcQMatrix();
		void							recalcQMatrixImpl();
		void							clearAllExceptFreqsAndRates();
		void							redimension(unsigned new_dim);
		std::string						showMatrixImpl(const double * q) const;

	private:

		unsigned						dimension;		/**< Number of elements in any given row or column of Q matrix (e.g. 4 for a 4x4 Q matrix) */
		std::vector<int>				dim_vect;		/**< Vector of length 2 describing the shape of the Q matrix for purposes of sending matrices across to Python as NumArray objects. For a 4x4 Q matrix, dim_vect[0] would be 4 and dim_vect[1] would be 4. */

		unsigned						flat_length;	/**< Number of elements in Q matrix (e.g. 16 for a 4x4 Q matrix) */

		std::vector<double>				pi;				/**< The state frequencies (length equals dimension) */
		std::vector<double>				sqrtPi;			/**< The square roots of the state frequencies (length equals dimension) */

		std::vector<double>				rr;				/**< The relative rates (elements in the upper diagonal of the R matrix). If the R matrix is 4x4, the order of the six elements in the relrates vector should be R[0][1], R[0][2], R[0][3], R[1][2], R[1][3] and R[2][3]. The R matrix is combined with the pi vector to create the Q matrix. */

		double							edgelen_scaler;	/**< */
		double * *						qmat;			/**< */

		double *						qmat_begin;		/**< */
		double *						qmat_end;		/**< */

		double *						w;				/**< */
		double *						w_begin;		/**< */
		double *						w_end;			/**< */

		double * *						z;				/**< */
		double *						z_begin;		/**< */
		double *						z_end;			/**< */

		double *						fv;				/**< */

		bool							q_dirty;		/**< If true, then the Q-matrix has been changed and eigenvalues and eigenvectors need to be recalculated */
	};
}	// namespace phycas

#include "pyphy/src/q_matrix.inl"

#endif
	