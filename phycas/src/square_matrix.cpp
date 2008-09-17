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
|	SquareMatrix constructor. Sets `dim' to `sz' and sets all elements of `m' to `value'.
*/
SquareMatrix::SquareMatrix(
  unsigned sz,		/**< the row and column dimension of the matrix */
  double value)		/**< the value to which to set all elements */
  : id(++k), m(0), dim(sz)
	{
	//std::cerr << "----> constructing SquareMatrix " << id << " of size " << dim << " with all values set to " << value << " <----" << std::endl;
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

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates memory for `m' data member and calls SquareMatrix::Identity.
*/
void SquareMatrix::CreateMatrix(
  unsigned sz,		/**< is the row and column dimension of the square matrix to be created */
  double value)		/**< is the value to which each element will be set */ 
	{
	//std::cerr << "----> SquareMatrix::CreateMatrix " << id << ", sz = " << sz << ", value = " << value << " <----" << std::endl;
	if (m)
		DeleteTwoDArray<double>(m);
	m = 0;
	dim = 0;
	if (sz > 0)
		{
		m = NewTwoDArray<double>(sz, sz);
		dim = sz;
		Fill(value);
		}
	}

void SquareMatrix::Fill(double value)
	{
	unsigned last = dim*dim;
	double * p = &m[0][0];
	for (unsigned i = 0; i < last; ++i)
		*p++ = value;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Assumes `dim' > 0 and `m' exists. Converts `m' to the identity matrix.
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
|	Assumes `dim' > 0 and `m' exists. Multiplies each element of `m' by the supplied `scalar'.
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
|	Assumes `dim' > 0 and `m' exists. Subtracts each element of `other' from the corresponding element of this.
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
|	Assumes `dim' > 0 and `m' exists. Sets each element of `m' to the supplied `scalar'.
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
	s += str(boost::format("This is SquareMatrix %d\n") % id);
	for (unsigned i = 0; i < dim; ++i)
		{
		for (unsigned j = 0; j < dim; ++j)
			{
			s += str(boost::format(fmt) % m[i][j]);
			}
		s += '\n';
		}
	}

} // namespace phycas
