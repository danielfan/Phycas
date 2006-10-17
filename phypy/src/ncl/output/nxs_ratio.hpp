/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

#ifndef NCL_NXSRATIO_H
#define NCL_NXSRATIO_H

#include <cfloat>

/*----------------------------------------------------------------------------------------------------------------------
|	NxsRatio provides a way to compute means of ratios that are numerically stable. A NxsRatio object stores the
|	numerator and denominator of the ratio separately, only dividing them when necessary.
*/
class NxsRatio
	{
	public:
					NxsRatio(double initTop = 0, double initBottom = 1.0);
					NxsRatio(const NxsRatio &r);

		NxsRatio	&operator +=(const NxsRatio &r);
		NxsRatio	&operator +=(double d);

		NxsRatio	&operator /=(const NxsRatio &r);
		NxsRatio	&operator /=(double d);

		bool		operator ==(const NxsRatio &r) const;
		bool		operator !=(const NxsRatio &r) const;
		
		double		GetRatioAsDouble() const;
		void		Invert();

		double		top;		/* the numerator */
		double		bottom;		/* the denominator */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `top' to 0.0 and `bottom' to 1.0.
*/
inline NxsRatio::NxsRatio(
  double initTop, 
  double initBottom)
  	:top(initTop),
  	bottom(initBottom)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `top' to `r.top' and `bottom' to `r.bottom'.
*/
inline NxsRatio::NxsRatio(
  const NxsRatio &r)	/* the NxsRatio object to be copied */
	{
	top		= r.top;
	bottom	= r.bottom;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds `r.top' to `top' and `r.bottom' to `bottom', returning reference to this object.
*/
inline NxsRatio &NxsRatio::operator +=(
  const NxsRatio &r)	/* the NxsRatio object to be added to this object */
	{
	top += r.top; 
	bottom += r.bottom; 
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds the supplied value `d' to both `top' and `bottom', returning a reference to this object.
*/
inline NxsRatio &NxsRatio::operator +=(
  double d)	/* the double value to be added to both numerator and denominator */
	{
	top += d; 
	bottom += d; 
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Divides `top' by `r.top' and `bottom' by `r.bottom'. returning reference to this object.
*/
inline NxsRatio &NxsRatio::operator /=(
  const NxsRatio &r)	/* the NxsRatio object used to perform the divisions */
	{
	top /= r.top; 
	bottom /= r.bottom; 
	return *this; 
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Divides both `top' and `bottom' by `d'. returning reference to this object.
*/
inline NxsRatio &NxsRatio::operator /=(
  double d)	/* the double value to be divided into both numerator and denominator */
	{ 
	top /= d; 
	bottom /= d; 
	return *this; 
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true iff `top' equals `r.top' and `bottom' equals `r.bottom'. Note that 1/2 is not equal to 2/4 if the 
|	NxsRatio == operator is used for the comparison! 
*/
inline bool NxsRatio::operator ==(
  const NxsRatio &r)	/* the NxsRatio object for comparison */
  const
	{
	return (top == r.top && bottom == r.bottom);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if either `top' is not equal to `r.top' or `bottom' is not equal to `r.bottom'.
*/
inline bool NxsRatio::operator !=(
  const NxsRatio &r)	/* the NxsRatio object for comparison */
  const
	{
	return (top != r.top || bottom != r.bottom);
	}
	

/*----------------------------------------------------------------------------------------------------------------------
|	Divides `top' by `bottom' and returns the floating point result. If `bottom' is exactly zero, then the largest 
|	possible positive (if `top' is positive) or negative (if `top' is negative) double value is returned.
*/
inline double NxsRatio::GetRatioAsDouble()
  const
	{
	double retval = 0.0;
	if (bottom != 0.0)
		retval = top / bottom;
	else if (top < 0.0)
		retval = -DBL_MAX;
	else if (top > 0.0)
		retval = DBL_MAX;

	return retval;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Swaps `top' and `bottom', converting the NxsRatio object to its inverse.
*/
inline void NxsRatio::Invert()
	{
	std:: swap<double>(top, bottom);
	}

#endif
