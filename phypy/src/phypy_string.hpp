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

#if !defined(PYPHY_STRING_HPP)
#define PYPHY_STRING_HPP

//POL 29-Nov-2005
//
// This file was originally intended to be a replacement for string_extensions.hpp and a workspace
// for trying new things (e.g. boost::format) without tampering with working code that depends on
// string_extensions.hpp. Unfortunately, there are too many places in the phypy code where I needed
// to include headers from ncl, which eventually leads to a conflict between string_extensions.hpp
// and phypy_string.hpp. So, now this file is largely useless and should be eliminated. I still like
// my DoubleFormatter function, however. It does the same job as string_extensions' DblFormatter,
// but uses boost::format (which is slower than sprintf but has the benefit of being typesafe).
// Thus, I have kept DoubleFormatter but now include string_extensions.hpp to get the other stuff.

#include <string>
#include <limits>
#include "boost/format.hpp"
#include "ncl/misc/string_extensions.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|	Functor that takes two std::string arguments (`a' and `b') and returns true if and only if the length of string `a'
|	is less than the length of string `b'. Can be used with algorithms std::min_element (or std::max_element) to find
|	the shortest (or longest) string in a vector, for example.
*/
class StringLengthLess
	{
	public:
		bool operator()(const std::string & a, const std::string & b) const
			{
			return (a.length() < b.length());
			}
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Class used with operator<<() to format doubles by providing the arguments needed by sprintf() to append a double 
|	value to the supplied std::string `s'. The first constructor parameter (`widthOfField') is the minimum width of the 
|	double; the second parameter (`nDigitsAfterDecimal') is the number of digits of precision to use.
*/
class DoubleFormatter
	{
	public:

						DoubleFormatter(unsigned widthOfField, unsigned nDigitsAfterDecimal);

		// Accessors
		//
		unsigned		GetWidth() const;
		unsigned		GetPrecision() const;

		// Modifiers
		//
		void			SetWidth(unsigned int widthOfField);
		void			SetPrecision(unsigned int nDigitsPrecision);

		// Utilites
		//
		std::string &	FormatDouble(std::string &s, double d);		

	protected:

		unsigned fieldWidth;					/**< is the number of characters to use when displaying the floating point value */
		unsigned digitsAfterDecimal;			/**< is the number of decimal places to show in the floating point value */
		boost::format formatter;				/**< is the formatting object */

	private:

		static const unsigned int	kMaxDoublePrintPrecision = std::numeric_limits<double>::digits10; 
		void						RecreateFormatCode();
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor that returns value of data member `fieldWidth'.
*/
inline unsigned DoubleFormatter::GetWidth() const
	{
	return fieldWidth;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor that returns value of data member `digitsAfterDecimal'.
*/
inline unsigned DoubleFormatter::GetPrecision() const
	{
	return digitsAfterDecimal;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier that sets the value of data member `fieldWidth' to the supplied value `widthOfField', then calls 
|	RecreateFormatCode().
*/
inline void DoubleFormatter::SetWidth(
  unsigned widthOfField)	/**< is the new field width */
	{
	fieldWidth = widthOfField;
	RecreateFormatCode();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier that sets the value of data member `digitsAfterDecimal' to the supplied value `nDigitsPrecision', then
|	calls RecreateFormatCode(). Note: this function ensures that `digitsAfterDecimal' will be no larger than
|	`kMaxDoublePrintPrecision'.
*/
inline void DoubleFormatter::SetPrecision(
  unsigned nDigitsPrecision)	/**< is the new precision */
	{
	PHYCAS_ASSERT(nDigitsPrecision < fieldWidth - 1);
	digitsAfterDecimal =  (unsigned) (nDigitsPrecision > kMaxDoublePrintPrecision ? kMaxDoublePrintPrecision : nDigitsPrecision);
	RecreateFormatCode();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The first constructor parameter (`widthOfField') is the minimum width of the double; the second parameter
|	(`nDigitsAfterDecimal') is the number of digits of precision to use.
*/
inline DoubleFormatter::DoubleFormatter(
  unsigned widthOfField,		/**< is the field width */
  unsigned nDigitsAfterDecimal)	/**< is the precision */
	: fieldWidth(widthOfField), digitsAfterDecimal(UINT_MAX)
	{
	SetPrecision(nDigitsAfterDecimal); //crops precision if needed and calls RecreateFormatCode()
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of the string `formatCode' sent to sprintf() function based on the values of `fieldWidth' and 
|	`digitsAfterDecimal'.
*/
inline void DoubleFormatter::RecreateFormatCode()
	{
	//@ perhaps we should be using %.g here?
	std::string formatStr;
	if (fieldWidth == UINT_MAX)
		{
		if (digitsAfterDecimal == UINT_MAX)
			formatStr = "%f"; //strcpy(formatCode, "%lf");
		else
			formatStr = str(boost::format("%%.%||f") % digitsAfterDecimal); //sprintf(formatCode, "%%.%dlf", digitsAfterDecimal);
		}
	else
		{
		if (digitsAfterDecimal == UINT_MAX)
			formatStr = str(boost::format("%%%||f") % fieldWidth); //sprintf(formatCode, "%%%dlf", fieldWidth);
		else
			formatStr = str(boost::format("%%%||.%||f") % fieldWidth % digitsAfterDecimal); //sprintf(formatCode, "%%%d.%dlf", fieldWidth, digitsAfterDecimal);
		}
	formatter = boost::format(formatStr);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline std::string &DoubleFormatter::FormatDouble(std::string &s, double d)
	{
	formatter % d;
	s.append(str(formatter));
	return s;
	}

std::string	&			append_unsigned(std::string &s, unsigned v);
std::string				abbreviate(const std::string &s, unsigned max_chars);

#if 0

template<typename T>
std::string &AppendNumber(std::string &s, const T i);

//class StringDblFormatterRef
//	{
//	public:
//		StringDblFormatterRef(std::string &s, DoubleFormatter &d) : strRef(s), dfRef(d) {}
//		std::string			&	strRef;
//		const DoubleFormatter		dfRef;
//	};

//StringDblFormatterRef	operator<<(std::string &, DoubleFormatter &);
//std::string				operator<<(StringDblFormatterRef &, const double d);

std::string &			operator<<(std::string &, int i);
std::string &			operator<<(std::string &, unsigned i);
std::string &			operator<<(std::string &, long l);
std::string &			operator<<(std::string &, unsigned long l);
std::string &			operator<<(std::string &, double d);
std::string &			operator<<(std::string &, const char *c);
std::string &			operator<<(std::string &, char c);
std::string &			operator<<(std::string &, const std::string &s);
bool					operator==(const std::string &s, const std::string &r);
bool					operator==(const std::string &s, const char *r);
bool					operator==(const std::string &s, const char &r);
bool					operator!=(const std::string &s, const char &r);
const std::string		operator+(const std::string &s, const char *r);
const std::string		operator+(const std::string &s, const std::string &r);

// moved up, out of #if 0
//std::string	&			append_unsigned(std::string &s, unsigned v);
//std::string				abbreviate(const std::string &s, unsigned max_chars);

template<>
inline std::string &AppendNumber<int>(std::string &s, const int i)
	{
	char tmp[81];
	std::sprintf(tmp, "%d", i);
	s.append(tmp);
	return s;
	}

template<>
inline std::string &AppendNumber<unsigned>(std::string & s, const unsigned i)
	{
	char tmp[81];
	std::sprintf(tmp, "%u", i);
	s.append(tmp);
	return s;
	}

template<>
inline std::string &AppendNumber<long>(std::string & s, const long l)
	{
	char tmp[81];
	std::sprintf(tmp, "%ld", l);
	s.append(tmp);
	return s;
	}

template<>
inline std::string &AppendNumber<unsigned long>(std::string & s, const unsigned long l)
	{
	char tmp[81];
	std::sprintf(tmp, "%lu", l);
	s.append(tmp);
	return s;
	}
	
template<>
inline std::string &AppendNumber<double>(std::string & s, const double d)
	{
	s = str(boost::format("%g") % d);
	return s;
	}

template<>
inline std::string & AppendNumber<float>(std::string & s, const float d)
	{
	s = str(boost::format("%g") % d);
	return s;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber(written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,	/**< is the string */
  int i)			/**< is the integer to append */
	{
	return AppendNumber<int>(s, i);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  unsigned i)   /* the unsigned integer to append */
	{
	return AppendNumber<unsigned>(s, i);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  long l)       /* the long integer to append */
	{
	return AppendNumber<long>(s, l);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  unsigned long l)      /* the unsigned long integer to append */
	{
	return AppendNumber<unsigned long>(s, l);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  double d)     /* the double floating point value to append */
	{
	return AppendNumber<double>(s, d);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the append operator (written to make it possible to use a std::string like an ostream)
*/
inline std::string & operator<<(
  std::string & s,
  const char * c)	/* the C-string to append */
	{
	s.append(c);
	return s;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the += operator (written to make it possible to use a std::string like an ostream)
*/
inline std::string & operator<<(
  std::string & s,
  char c)       /* the char to append */
	{
	s.push_back(c); //std namespace's operator+=
	return s;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the += operator (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & f,
  const std::string &s)   /* the std::string to append */
	{
	f.append(s); //std namespace's operator+=
	return f;
	}

#endif

#endif
