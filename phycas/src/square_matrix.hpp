/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
|																			  |
|  This program is free software; you can redistribute it and/or modify		  |
|  it under the terms of the GNU General Public License as published by		  |
|  the Free Software Foundation; either version 2 of the License, or		  |
|  (at your option) any later version.										  |
|																			  |
|  This program is distributed in the hope that it will be useful,			  |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of			  |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the			  |
|  GNU General Public License for more details.								  |
|																			  |
|  You should have received a copy of the GNU General Public License along	  |
|  with this program; if not, write to the Free Software Foundation, Inc.,	  |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.				  |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#if ! defined(SQUARE_MATRIX_HPP)
#define SQUARE_MATRIX_HPP

//#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include <boost/shared_ptr.hpp>

extern "C"
{
#include "linalg.h"
}

namespace phycas
{

void fillTranspose(double ** p_mat_trans_scratch, const double * const *p_mat, unsigned dim);

class SquareMatrix
	{
	public:
										SquareMatrix();
										SquareMatrix(unsigned sz, double value);
										SquareMatrix(const SquareMatrix & other);
										~SquareMatrix();
		void							CreateMatrix(unsigned sz, double value);
		void							Identity();
		void							ScalarMultiply(double scalar);
		void							Subtract(const SquareMatrix & other);
		void							SetToScalar(double scalar);
		unsigned						GetDimension() const;
		double                          Trace() const;
		double * *						GetMatrixAsRawPointer() const;
		double *						operator[](unsigned i) const;
		void							MatrixToString(std::string & s, std::string fmt) const;
		void							Fill(double value);
        std::string                     GetStringRepresentation() const;
        std::string                     GetFormattedStringRepresentation(std::string fmt) const;
        double                          GetElement(unsigned i, unsigned j) const;
        void                            SetElement(unsigned i, unsigned j, double v);
        std::vector<double>             GetMatrix() const;
        void                            SetMatrix(unsigned sz, std::vector<double>);
        SquareMatrix *                  Inverse() const;
        SquareMatrix *                  LeftMultiply(SquareMatrix & matrixOnLeft) const;
        SquareMatrix *                  RightMultiply(SquareMatrix & matrixOnRight) const;
                
	protected:
		static unsigned					k;		//temporary
		unsigned						id;		//temporary
		double * *						m;		/**< the two-dimensional matrix of doubles */
		unsigned						dim;	/**< the dimension of the matrix */
	};
    
typedef boost::shared_ptr<SquareMatrix> SquareMatrixShPtr;


}	// namespace phycas
#endif
