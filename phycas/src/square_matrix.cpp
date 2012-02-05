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

#include "phycas/src/square_matrix.hpp"

#include <boost/format.hpp>

#include "ncl/nxsallocatematrix.h"


namespace phycas
{

unsigned SquareMatrix::k = 0;	//temporary!

/*----------------------------------------------------------------------------------------------------------------------
|	SquareMatrix default constructor. Sets `dim' and `m' data members to 0.
*/
SquareMatrix::SquareMatrix()
  :	 id(++k),  m(0), dim(0)
	{
	//std::cerr << "----> constructing default SquareMatrix " << id << " <----" << std::endl;
	} 

/*----------------------------------------------------------------------------------------------------------------------
|	SquareMatrix constructor. Creates square matrix `m' with dimension `sz' and sets all elements of `m' to `value'.
*/
SquareMatrix::SquareMatrix(
  unsigned sz,		/**< the row and column dimension of the matrix */
  double value)		/**< the value to which to set all elements */
  : id(++k), m(0), dim(0)
	{
	//std::cerr << "----> constructing SquareMatrix " << id << " of size " << sz << " with all values set to " << value << " <----" << std::endl;
	CreateMatrix(sz, value);
	} 

/*----------------------------------------------------------------------------------------------------------------------
|	SquareMatrix copy constructor. This is necessary because the standard vector resize function in some implementations
|	of the standard template library creates only one element with the default constructor, then creates the remaining
|	elements using the copy constructor.
*/
SquareMatrix::SquareMatrix(
  const SquareMatrix & other)	/**< is the SquareMatrix to copy */
  : id(++k), m(0), dim(other.dim)
	{
	//std::cerr << "----> copy constructing SquareMatrix " << id << " from " << other.id << " <----" << std::endl;
	if (dim > 0)
		{
		CreateMatrix(dim, 0.0);
		double * pother = &other.m[0][0];
		double * p = &m[0][0];
		unsigned last = dim*dim;
		for (unsigned i = 0; i < last; ++i)
			*p++ = *pother++;
		}
	} 

/*----------------------------------------------------------------------------------------------------------------------
|	SquareMatrix destructor. If `m' exists, deletes it.
*/
SquareMatrix::~SquareMatrix() 
	{
	//std::cerr << "----> destroying SquareMatrix " << id << " <----" << std::endl;
	if (m) 
		DeleteTwoDArray<double>(m);
	m = 0;
	dim = 0;
	id = UINT_MAX;
	}
    
//static double ** CreateNewTwoDArray(unsigned f, unsigned s)
//	{
//	double **ptr;
//	ptr = new double *[f];
//	*ptr = new double [f * s];
//	for (unsigned fIt = 1 ; fIt < f ; fIt ++)
//		ptr[fIt] = ptr[fIt -1] +  s ;
//	return ptr;
//	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates memory for `m' data member and calls SquareMatrix::Identity.
*/
void SquareMatrix::CreateMatrix(
  unsigned sz,		/**< is the row and column dimension of the square matrix to be created */
  double value)		/**< is the value to which each element will be set */ 
	{
	//std::cerr << "----> SquareMatrix::CreateMatrix " << id << ", sz = " << sz << ", value = " << value << " <----" << std::endl;
    dim = sz;
	if (m)
		DeleteTwoDArray<double>(m);
	m = 0;
	if (sz > 0)
		{
		m = NewTwoDArray<double>(sz, sz);
		Fill(value);
		}
        
	//std::cerr << "----> SquareMatrix::CreateMatrix " << id << ", sz = " << sz << ", value = " << value << " <----\n";
    //for (unsigned i = 0; i < sz; ++i)
    //    {
    //    for (unsigned j = 0; j < sz; ++j) 
    //        {
    //        std::cerr << boost::str(boost::format("%12.1f") % m[i][j]);
    //        }
    //    std::cerr << "\n";
    //    }
    //std::cerr << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Fills `p_mat_trans_scratch' with transpose of `p_mat'. Assumes both `p_mat_trans_scratch' and `p_mat' are allocated
|   square matrices of doubles and are of the same dimension. Overwrites any values currently in `p_mat_trans_scratch'.
*/
void fillTranspose(double ** p_mat_trans_scratch, const double * const *p_mat, unsigned dim)
	{
	for (unsigned j = 0; j < dim; ++j)
		{
		for (unsigned k = 0 ; k < dim; ++k)
			p_mat_trans_scratch[j][k] = p_mat[k][j];
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the dimension of the data member `m' (the row and column dimension are the same).
*/
unsigned SquareMatrix::GetDimension() const 
	{
	return dim;
	}
    
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the sum of the elements on the main diagonal.
*/
double SquareMatrix::Trace() const 
	{
    double trace = 0.0;
    for (unsigned i = 0; i < dim; ++i)
        trace += m[i][i];
	return trace;
	}
    
/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that simply returns `m'. This allows `m' to be passed to functions that expect a two-dimensional
|	array of doubles, but keep in mind that this is somewhat unsafe.
*/
double * * SquareMatrix::GetMatrixAsRawPointer() const 
	{
	return m;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Operator that allows access to this SquareMatrix object as if it were a two-dimensional array of doubles.
*/
double * SquareMatrix::operator[](
  unsigned i) const 
	{
	PHYCAS_ASSERT(i < GetDimension()); 
	return m[i];
	}
    
/*----------------------------------------------------------------------------------------------------------------------
|	Returns string representation of the matrix.
*/
std::string SquareMatrix::GetStringRepresentation() const 
	{
    std::string s;
    std::string fmt = "\t%g";
    MatrixToString(s, fmt);
    return s;
	}
    
/*----------------------------------------------------------------------------------------------------------------------
|	Returns string representation of the matrix.
*/
std::string SquareMatrix::GetFormattedStringRepresentation(std::string fmt = "%g") const 
	{
    std::string s;
    MatrixToString(s, fmt);
    return s;
	}
    
/*----------------------------------------------------------------------------------------------------------------------
|	Returns element at row `i', column `'j of the matrix.
*/
double SquareMatrix::GetElement(unsigned i, unsigned j) const 
	{
	PHYCAS_ASSERT(i < dim); 
	PHYCAS_ASSERT(j < dim); 
    return m[i][j];
	}
    
/*----------------------------------------------------------------------------------------------------------------------
|	Sets element at row `i', column `j' of the matrix to the value `v'.
*/
void SquareMatrix::SetElement(unsigned i, unsigned j, double v) 
	{
	PHYCAS_ASSERT(i < dim); 
	PHYCAS_ASSERT(j < dim); 
    m[i][j] = v;
	}
    
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the entire matrix as a vector by row.
*/
std::vector<double> SquareMatrix::GetMatrix() const 
	{
    std::vector<double> v;
    double * ptr = (double *)&m[0];
    for (unsigned i = 0; i < dim*dim; ++i)
        v.push_back(ptr[i]);
    return v;
	}
    
/*----------------------------------------------------------------------------------------------------------------------
|	Copies supplied vector `v' into this matrix. Vector `v' is expected to have length `sz'*`sz', and matrix is read
|   from `v' by row. The data member `dim' is set equal to `sz'.
*/
void SquareMatrix::SetMatrix(unsigned sz, std::vector<double> v) 
	{
	PHYCAS_ASSERT(sz == (int)sqrt(v.size())); 
    CreateMatrix(sz, 0.0);
    double * ptr = (double *)&m[0][0];
    for (std::vector<double>::iterator it = v.begin(); it != v.end(); ++it)
        {
        double val = *it;
        *ptr++ = val;
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Fills all cells of matrix `m' with `value'.
*/
void SquareMatrix::Fill(
  double value) /**< is the value to which every element in matrix is to be set */
	{
    unsigned dim = GetDimension();
	unsigned last = dim*dim;
	double * p = &m[0][0];
	for (unsigned i = 0; i < last; ++i)
		*p++ = value;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assumes dimension > 0 and `m' exists. Converts `m' to the identity matrix.
*/
void SquareMatrix::Identity()
	{
	//std::cerr << "----> SquareMatrix::Identity " << id << " <----" << std::endl;
	PHYCAS_ASSERT(dim > 0);
	unsigned i;
	unsigned last = dim*dim;
	double * p = &m[0][0];
	for (i = 0; i < last; ++i)
		*p++ = 0.0;
	for (i = 0; i < dim; ++i)
		m[i][i] = 1.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assumes dimension > 0 and `m' exists. Multiplies each element of `m' by the supplied `scalar'.
*/
void SquareMatrix::ScalarMultiply(
  double scalar)	/**< is the scalar value multiplied by each element */
	{
	//std::cerr << "----> SquareMatrix::ScalarMultiply " << id << ", scalar = " << scalar << " <----" << std::endl;
	PHYCAS_ASSERT(dim > 0);
	unsigned last = dim*dim;
	double * p = &m[0][0];
	for (unsigned i = 0; i < last; ++i)
		*p++ *= scalar;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assumes dimension > 0 and `m' exists. Subtracts each element of `other' from the corresponding element of this.
*/
void SquareMatrix::Subtract(
  const SquareMatrix & other)	 /**< is the matrix to subtract from this */
	{
	PHYCAS_ASSERT(dim > 0);
	unsigned last = dim*dim;
	double * otherp = &other.m[0][0];
	double * p = &m[0][0];
	for (unsigned i = 0; i < last; ++i)
		*p++ -= *otherp++;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assumes dimension > 0 and `m' exists. Sets each element of `m' to the supplied `scalar'.
*/
void SquareMatrix::SetToScalar(
  double scalar)	/**< is the scalar value to which each element is set */
	{
	//std::cerr << "----> SquareMatrix::SetToScalar " << id << ", scalar = " << scalar << " <----" << std::endl;
	PHYCAS_ASSERT(dim > 0);
	unsigned last = dim*dim;
	double * p = &m[0][0];
	for (unsigned i = 0; i < last; ++i)
		*p++ = scalar;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Saves a representation of the matrix to the supplied string for use in debugging. The caller should supply a format
|	string suitable for representing a single value of the matrix: e.g. "%12.5f\t".
*/
void SquareMatrix::MatrixToString(
  std::string & s,		 /**< the string to which a representation will be saved */
  std::string fmt) const /**< the format string */
	{
    unsigned dim = GetDimension();
	//s += boost::str(boost::format("This is SquareMatrix %d\n") % id);
	for (unsigned i = 0; i < dim; ++i)
		{
		for (unsigned j = 0; j < dim; ++j)
			{
			s += str(boost::format(fmt) % m[i][j]);
			}
		s += '\n';
		}
	}
    
/*----------------------------------------------------------------------------------------------------------------------
|	Creates and returns a matrix that represents the inverse of this SquareMatrix.
*/
SquareMatrix * SquareMatrix::Inverse() const
    {
    double * col = new double[dim];
    int * permutation = new int[dim];
    SquareMatrix * tmp = new SquareMatrix(*this);
    double ** a = tmp->GetMatrixAsRawPointer();
    SquareMatrix * inv = new SquareMatrix(*this);
    double ** a_inv = inv->GetMatrixAsRawPointer();

    // double **  a           matrix represented as vector of row pointers 
    // int        n           order of matrix
    // double *   col         work vector of size n 
    // int *      permutation work vector of size n 
    // double **  a_inv       inverse of input matrix a (matrix a is destroyed) 
    int result = InvertMatrix(a, dim, col, permutation, a_inv);
    delete tmp;
    delete [] col;
    delete [] permutation;
    if (result != 0)
        {
        delete inv;
        return 0;
        }
    return inv;
    }
    
/*----------------------------------------------------------------------------------------------------------------------
|	Creates and returns a matrix that represents the product of `matrixOnLeft' with this matrix.
*/
SquareMatrix * SquareMatrix::LeftMultiply(SquareMatrix & matrixOnLeft) const
    {
    double ** l = matrixOnLeft.GetMatrixAsRawPointer();
    double ** r = GetMatrixAsRawPointer();
    SquareMatrix * tmp = new SquareMatrix(*this);
    double ** t = tmp->GetMatrixAsRawPointer();
    
    for (unsigned i = 0; i < dim; ++i)
        {
        for (unsigned j = 0; j < dim; ++j)
            {
            double vsum = 0.0;
            for (unsigned k = 0; k < dim; ++k)
                {
                vsum += l[i][k]*r[k][j];
                }
            t[i][j] = vsum;
            }
        }
    return tmp;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Creates and returns a matrix that represents the product of this matrix with `matrixOnRight'.
*/
SquareMatrix * SquareMatrix::RightMultiply(SquareMatrix & matrixOnRight) const
    {
    double ** l = GetMatrixAsRawPointer();
    double ** r = matrixOnRight.GetMatrixAsRawPointer();
    SquareMatrix * tmp = new SquareMatrix(*this);
    double ** t = tmp->GetMatrixAsRawPointer();
    
    for (unsigned i = 0; i < dim; ++i)
        {
        for (unsigned j = 0; j < dim; ++j)
            {
            double vsum = 0.0;
            for (unsigned k = 0; k < dim; ++k)
                {
                vsum += l[i][k]*r[k][j];
                }
            t[i][j] = vsum;
            }
        }
    return tmp;
    }
    
    
} // namespace phycas
