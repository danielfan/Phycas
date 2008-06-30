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

#if ! defined(SQUARE_MATRIX_HPP)
#define SQUARE_MATRIX_HPP

#include "phycas/src/cipres/AllocateMatrix.hpp"
	
namespace phycas
{

void fillTranspose(double ** p_mat_trans_scratch, const double * const *p_mat, const unsigned);

#if POLPY_NEWWAY
class SquareMatrix
    {
    public:
                                        SquareMatrix();
                                        SquareMatrix(unsigned sz, double value);
                                        SquareMatrix(const SquareMatrix & other);
                                        ~SquareMatrix();
        void                            CreateMatrix(unsigned sz, double value);
        void                            Identity();
        void                            ScalarMultiply(double scalar);
        void                            Subtract(const SquareMatrix & other);
        void                            SetToScalar(double scalar);
        unsigned                        GetDimension() const;
        double * *                      GetMatrix() const;
        double *                        operator[](unsigned i) const;
        void                            MatrixToString(std::string & s, std::string fmt) const;
                                        
    protected:
        static unsigned                 k;      //temporary
        unsigned                        id;     //temporary
        double * *                      m;      /**< the two-dimensional matrix of doubles */
        unsigned                        dim;    /**< dimension of both rows and columns */ 
    };

inline void fillTranspose(double ** p_mat_trans_scratch, const double * const *p_mat, const unsigned dim)
	{
	for (unsigned j = 0; j < dim; ++j)
		{
		for (unsigned k = 0 ; k < dim; ++k)
			p_mat_trans_scratch[j][k] = p_mat[k][j];
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that simply returns the value of `dim' (the row and column dimension).
*/
inline unsigned SquareMatrix::GetDimension() const 
    {
    return dim;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that simply returns `m'. This allows `m' to be passed to functions that expect a two-dimensional
|   array of doubles, but keep in mind that this is somewhat unsafe.
*/
inline double * * SquareMatrix::GetMatrix() const 
    {
    return m;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Operator that allows access to this SquareMatrix object as if it were a two-dimensional array of doubles.
*/
inline double * SquareMatrix::operator[](
  unsigned i) const 
    {
    PHYCAS_ASSERT(i < dim); 
    return m[i];
    }

#endif
}   // namespace phycas
#endif
