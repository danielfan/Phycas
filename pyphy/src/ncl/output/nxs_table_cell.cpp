//#include "phycas/force_include.h"
#include <climits>
#include <iostream>
#include <iomanip>
#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/misc/string_extensions.hpp"
#include "pyphy/src/ncl/output/nxs_ratio.hpp"
#include "pyphy/src/ncl/output/nxs_table_cell.hpp"
#include "pyphy/src/ncl/output/nxs_table.hpp"
#include "pyphy/src/ncl/output/nxs_output_stream.hpp"
using std::string;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::isspace;
	using std::sprintf;
	using std::ispunct;
#endif

unsigned 	NxsTableCell::def_width			=	0;
unsigned 	NxsTableCell::def_precision		=	1;
unsigned 	NxsTableCell::def_asterisks		=	3;
char		NxsTableCell::cell_workspace[ABSOLUTE_MAXIMUM_COLUMN_WIDTH+1];

/*----------------------------------------------------------------------------------------------------------------------
|	Sets 'justification' to NxsTableCell::nojustify, and `width', `precision', and `asterisks' to their default settings
|	(NxsTableCell::def_width, NxsTableCell::def_precision, and NxsTableCell::def_asterisks, respectively). Sets `s' to 
|	the empty string, `r' to the default ratio, `i' to 0, `l' to 0L, `d' to 0.0, and `resume_at' to 0.
*/
void NxsTableCell::FactorySettings()
	{
	resume_at = 0;
	s.clear();
	r = NxsRatio();
	l = 0L;
	i = 0;
	d = 0.0;
	justification	= NxsTableCell::nojustify;
	width			= NxsTableCell::def_width;
	precision		= NxsTableCell::def_precision;
	asterisks		= NxsTableCell::def_asterisks;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First calls FactorySettings(), then sets the cell type `tag' to NxsTableCell::I and `i' to 0.
*/
NxsTableCell::NxsTableCell()
	{
	FactorySettings();
	tag	= NxsTableCell::I;
	i	= 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First calls FactorySettings(), then sets the cell type `tag' to NxsTableCell::S and `s' to the supplied value `ss'.
*/
NxsTableCell::NxsTableCell(
  const string &ss)	/* the string to which this object will be set */
	{
	FactorySettings();
	tag	= NxsTableCell::S;
	s	= ss;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First calls FactorySettings(), then sets the cell type `tag' to D (double) and `d' to the supplied value `dd'.
*/
NxsTableCell::NxsTableCell(
  const double dd)	/* the floating point value to which this object will be set */
	{
	FactorySettings();
	tag	= NxsTableCell::D;
	d = dd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First calls FactorySettings(), then sets the cell type `tag' to NxsTableCell::I and `i' to the supplied value `ii'.
*/
NxsTableCell::NxsTableCell(
  const int ii)	/* the integer to which this object will be set */
	{
	FactorySettings();
	tag	= NxsTableCell::I;
	i	= ii;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First calls FactorySettings(), then sets the cell type `tag' to NxsTableCell::L and `l' to the supplied value `ll'.
*/
NxsTableCell::NxsTableCell(
  const long ll)	/* the long integer to which this object will be set */
	{
	FactorySettings();
	tag	= NxsTableCell::L;
	l	= ll;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First calls FactorySettings(), then sets the cell type `tag' to NxsTableCell::R and `r' to the supplied value `rr'.
*/
NxsTableCell::NxsTableCell(
  const NxsRatio rr)	/* the NxsRatio to which this object will be set */
	{
	FactorySettings();
	tag	= NxsTableCell::R;
	r	= rr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First calls FactorySettings(), then sets the cell type `tag' to NxsTableCell::S and `s' to the supplied value `ss'.
|	Returns *this.
*/
NxsTableCell &NxsTableCell::operator =(
  const string &ss)	/* the string to which this object will be set */
	{
	FactorySettings();
	tag	= NxsTableCell::S;
	s	= ss;
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First calls FactorySettings(), then sets the cell type `tag' to NxsTableCell::D and `d' to the supplied value `dd'.
|	Returns *this.
*/
NxsTableCell &NxsTableCell::operator =(
  const double dd)	/* the floating point value to which this object will be set */
	{
	FactorySettings();
	tag	= NxsTableCell::D;
	d	= dd;
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First calls FactorySettings(), then sets the cell type `tag' to NxsTableCell::I and `i' to the supplied value `ii'.
|	Returns *this.
*/
NxsTableCell &NxsTableCell::operator =(
  const int ii)	/* the integer value to which this object will be set */
	{
	FactorySettings();
	tag	= NxsTableCell::I;
	i	= ii;
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First calls FactorySettings(), then sets the cell type `tag' to NxsTableCell::L and `l' to the supplied value `ll'.
|	Returns *this.
*/
NxsTableCell &NxsTableCell::operator =(
  const long ll)	/* the long integer value to which this object will be set */
	{
	FactorySettings();
	tag	= NxsTableCell::L;
	l	= ll;
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First calls FactorySettings(), then sets the cell type `tag' to NxsTableCell::R and `r' to the supplied value `rr'.
|	Returns *this.
*/
NxsTableCell &NxsTableCell::operator =(
  const NxsRatio rr)	/* the ratio to which this object will be set */
	{
	FactorySettings();
	tag	= NxsTableCell::R;
	r	= rr;
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Compares this NxsTableCell with `cell' using a criterion based on the cell type (`tag'). If this NxsTableCell's 
|	`tag' is not NxsTableCell::D, NxsTableCell::I, NxsTableCell::L, or NxsTableCell::R, returns true. If the `tag' data
|	member of `cell' is not NxsTableCell::D, NxsTableCell::I, NxsTableCell::L, or NxsTableCell::R, returns false. 
|	Assuming the operator did not return because of either of the above tests, assigns the value of this NxsTableCell
|	to the double `thisvalue' and that of cell to the double `cellvalue'. The operator then returns true if `thisvalue'
|	is less than `cellvalue'.
*/
bool NxsTableCell::operator <(
  const NxsTableCell &cell)	/* the NxsTableCell for comparison */
	{
	//@pol need to find where this operator is used and see if a less confusing definition is possible
	//	(e.g. one that does a lexical comparison for strings, so that the operator works as expected no
	//	matter what the type of the cell
	//
	double this_value = 0.0;
	bool this_is_numeric = true;
	if (tag == NxsTableCell::D)
		this_value = d;
	else if (tag == NxsTableCell::I)
		this_value = static_cast<double>(i);
	else if (tag == NxsTableCell::L)
		this_value = static_cast<double>(l);
	else if (tag == NxsTableCell::R)
		this_value = r.GetRatioAsDouble();
	else
		this_is_numeric = false;

	if (!this_is_numeric) 
		return true;

	double cell_value = 0.0;
	bool cell_is_numeric = true;
	if (cell.tag == NxsTableCell::D)
		cell_value = cell.d;
	else if (cell.tag == NxsTableCell::I)
		cell_value = static_cast<double>(cell.i);
	else if (cell.tag == NxsTableCell::L)
		cell_value = static_cast<double>(cell.l);
	else if (cell.tag == NxsTableCell::R)
		cell_value = cell.r.GetRatioAsDouble();
	else
		cell_is_numeric = false;

	if (!cell_is_numeric)
		return false;

	return (this_value < cell_value ? true : false);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the data member `s', which is made to sontain a string of `n' characters, each one being 
|	the specified fill character `ch'. If `n' is greater than `width' minus `NxsTable::colspacer', the number of fill 
|	characters output is cut to `width' minus `NxsTable::colspacer'. 
*/
string &NxsTableCell::FillWith(
  unsigned n,	/* the number of fill characters */
  char ch)		/* the fill character */
	{
	unsigned effective_width = (width - NxsTable::colspacer);
	if (n > effective_width) 
		n = effective_width;

	s.clear();
	s.append(n, ch);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Outputs to `out' a string filled with exactly `width' blank space characters.
*/
void NxsTableCell::ShowEmpty(NxsOutputStream &out) const
	{
	out << std::string(width, ' ');
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Outputs to `out' the contents of this NxsTableCell, formatted appropriately. Returns true if the contents of this 
|	cell fit within the `width' allotted. If, for example, the cell type is a string and the string is longer than 
|	`width', the string is broken (at a punctuation character if possible), and the next character after the break 
|	position is stored in `resume_at'. The substring up to `resume_at' is output, and the function returns false to 
|	indicate that not all of the string could be printed. The remainder of the string will be output on a subsequent
|	line. The following list shows what is output for each possible cell type:
|~
|	o for NxsTableCell::D cells, the data member `d' is output
|	o for NxsTableCell::I cells, the data member `i' is output
|	o for NxsTableCell::L cells, the data member `l' is output
|	o for NxsTableCell::R cells, the value `r'.GetRatioAsDouble() is output
|	o for NxsTableCell::S cells, the data member `s' is output
|	o for NxsTableCell::E cells, the member function ShowEmpty(out) is called
|	o for NxsTableCell::B cells, the member function FillWith(`i', '*') is called
|	o for NxsTableCell::H cells, the member function FillWith(width - 1, '-') is called
|	o for NxsTableCell::N cells, the data member `l' is output
|	o for NxsTableCell::M cells, the data member `d' is output
|	o for NxsTableCell::SUM cells, the data member `d' is output
|	o for NxsTableCell::SS cells, the data member `d' is output
|	o for NxsTableCell::CV cells, the data member `d' is output
|	o for NxsTableCell::CIU cells, the data member `d' is output
|	o for NxsTableCell::CIL cells, the data member `d' is output
|	o for NxsTableCell::CIP cells, the data member `i' is output
|	o for NxsTableCell::SD cells, the data member `d' is output
|	o for NxsTableCell::V cells, the data member `d' is output
|~
*/
bool NxsTableCell::Show(NxsOutputStream &out)
{
	//@pol Previously, I had noted "While it is tempting to output simply s (rather than s.c_str()) for the S type,
	// the g++ compiler doesn't honor the width() call for std::string." Is this still true?
	//
	unsigned k, kk, m;
	unsigned punct_tol, slen, wslen, remaininglen, watch_after, last_char, ch, p;
	char cStringBuffer[256], f_str[256];
	//@pol Need to be able to set precision
	// out.precision(precision);
	// out.setf(ios::fixed, ios::floatfield);
	// out.setf(ios::showpoint);

	unsigned colwidth = width - NxsTable::colspacer;
	NXS_ASSERT(colwidth <= ABSOLUTE_MAXIMUM_COLUMN_WIDTH);

	// Return value is true unless we haven't finished with this cell
	// which is only the case if the cell contains a string that is
	// longer than the cell's width. A string will only be longer than
	// the cell's width if it is longer than ABSOLUTE_MAXIMUM_COLUMN_WIDTH.
	//
	bool done = true;

	switch(tag)
		{
		case NxsTableCell::S :
			// If a string is too long to print, we chop it into pieces, at punctuation
			// if possible. But we want to work back only so far from the chopping point in
			// search of punctuation. In the worst case scenario, there are no blank spaces,
			// colons, semicolons, commas, periods, etc. in the string at all, and we can't
			// afford to not print some of the string! punct_tol thus holds the number of
			// characters we are allowed to chop back from the break point in order to break
			// at a punctuation character.
			//
			punct_tol = (colwidth / 3);
			ch = 0;
			slen = (unsigned)s.length();
			last_char = resume_at + colwidth;
			watch_after = last_char - punct_tol;
			remaininglen = (unsigned)(s.length() - resume_at);
			if (remaininglen + resume_at < last_char)
				{
				last_char = remaininglen + resume_at;
				watch_after = last_char;
				}

			if (resume_at > (unsigned)0 || remaininglen > colwidth)
				{
				// First, adjust resume_at so that any leading blank space characters will be skipped
				//
				bool leading_spaces = true;
				for (m = resume_at; leading_spaces && (m < last_char); ++m)
					{
					if (s[m] == ' ')
						++resume_at;
					else
						leading_spaces = false;
					}
						
				// Fill the workspace with characters up to entire width of column
				//
				wslen = 0;
				cell_workspace[0] = '\0';
				for (m = resume_at; m < last_char; ++m)
					{
					cell_workspace[wslen++] = s[m];
					cell_workspace[wslen] = '\0';
					}

				// Now work backwards from the end looking for punctuation or space
				// characters, but don't look at more than punct_tol characters.
				// If no punctuation found in last punct_tol characters, leave
				// cell_workspace as it is.
				//
				p = m - resume_at - 1;
				for (k = 0; k < punct_tol; ++k, --p)
					{
					ch = cell_workspace[p];
					if (isspace((int) ch) || ispunct((int) ch))
						{
						last_char -= k;
						wslen -= k;
						cell_workspace[p+1] = '\0';
						break;
						}
					}
				
				// At this point, we've either found a punctuation or whitespace character
				// (in which case last_char is now set to point to the position in s just
				// beyond that punctuation or space character), or we've come up empty (in
				// which case last_char is still set to resume_at + colwidth.
				//
				if (justification == NxsTableCell::leftjustify)
					out << cell_workspace << std::string(width - last_char + resume_at, ' ');
				else if (justification == NxsTableCell::centerjustify)
					{
					k = wslen + (colwidth - wslen) / 2;
					out << MakeRightJustifiedString(cell_workspace, k);
					out << std::string((colwidth - k) + NxsTable::colspacer, ' ');
					}
				else
					{
					NXS_ASSERT(justification != NxsTableCell::nojustify);
					out << MakeRightJustifiedString(cell_workspace, colwidth);
					out << std::string(NxsTable::colspacer, ' ');
					}

				// If we've just output the last of the string, then reset
				// resume_at to 0, otherwise return false to indicate there
				// is more of the string left to print
				//
				if (m < (unsigned)s.length())
					{
					resume_at = last_char; 
					done = false;
					}
				else
					resume_at = 0;
				}

			// String fits into width, so go ahead and output the whole thing
			//
			else if (justification == NxsTableCell::leftjustify)
				{
				out << s;

				k = width - slen;
				out << string(k, ' ');
				}
			else if (justification == NxsTableCell::centerjustify)
				{
				k = slen + (colwidth - slen) / 2;
				out << MakeRightJustifiedString(s, k);
				out <<  string(width - k, ' ');
				}
			else
				{
				NXS_ASSERT(justification != NxsTableCell::nojustify);
				out << MakeRightJustifiedString(s, colwidth);
				out << string(NxsTable::colspacer, ' ');
				}
			break;
		
		case NxsTableCell::D :
		case NxsTableCell::SUM :
		case NxsTableCell::SS :
		case NxsTableCell::M :
		case NxsTableCell::V :
		case NxsTableCell::SD :
		case NxsTableCell::CV :
		case NxsTableCell::CIU :
		case NxsTableCell::CIL :
			sprintf(f_str, "%%.%df", precision);
			sprintf(cStringBuffer, f_str, d);
			kk = (unsigned)strlen(cStringBuffer);
			if (justification == NxsTableCell::leftjustify)
				out << cStringBuffer << string(width - kk, ' ');
			else if (justification == NxsTableCell::centerjustify)
				{
				k = kk + (colwidth - kk) / 2;
				out << MakeRightJustifiedString(cStringBuffer, k) << string(width - k, ' ');
				}
			else
				{
				NXS_ASSERT(justification != NxsTableCell::nojustify);
				out << MakeRightJustifiedString(cStringBuffer, colwidth) <<  string(NxsTable::colspacer, ' ');
				}
			break;
		
		case NxsTableCell::R :
			sprintf(f_str, "%%.%df", precision);
			sprintf(cStringBuffer, f_str, r.GetRatioAsDouble());
			kk = (unsigned)strlen(cStringBuffer);
			if (justification == NxsTableCell::leftjustify)
				out << cStringBuffer << string(width - kk, ' ');
			else if (justification == NxsTableCell::centerjustify)
				{
				k = kk + (colwidth - kk) / 2;
				out << MakeRightJustifiedString(cStringBuffer, k) << string(width - k, ' ');
				}
			else
				{
				NXS_ASSERT(justification != NxsTableCell::nojustify);
				out << MakeRightJustifiedString(cStringBuffer, colwidth) << string(NxsTable::colspacer, ' ');
				}
			break;
		
		case NxsTableCell::I :
		case NxsTableCell::CIP :
			sprintf(cStringBuffer, "%d", i);	//i not l
			kk = (unsigned)strlen(cStringBuffer);
			if (justification == NxsTableCell::leftjustify)
				out << cStringBuffer << string(width - kk, ' ');
			else if (justification == NxsTableCell::centerjustify)
				{
				k = kk + (colwidth - kk) / 2;
				out << MakeRightJustifiedString(cStringBuffer, k) << string(width - k, ' ');
				}
			else
				{
				NXS_ASSERT(justification != NxsTableCell::nojustify);
				out << MakeRightJustifiedString(cStringBuffer, colwidth) <<  string(NxsTable::colspacer, ' ');
				}
			break;
		
		case NxsTableCell::N :
		case NxsTableCell::L :
			sprintf(cStringBuffer, "%ld", l);
			kk = (unsigned)strlen(cStringBuffer);
			if (justification == NxsTableCell::leftjustify)
				out << cStringBuffer << string(width - kk, ' ');
			else if (justification == NxsTableCell::centerjustify)
				{
				k = kk + (colwidth - kk) / 2;
				out << MakeRightJustifiedString(cStringBuffer, k) << string(width - k, ' ');
				}
			else
				{
				NXS_ASSERT(justification != NxsTableCell::nojustify);
				out << MakeRightJustifiedString(cStringBuffer, colwidth) << string(NxsTable::colspacer, ' ');
				}
			break;
		
		case NxsTableCell::B :
			if (justification == NxsTableCell::leftjustify)
				{
				out << FillWith((unsigned) i, '*');
				out << FillWith(NxsTable::colspacer, ' ');
				}
			else if (justification == NxsTableCell::centerjustify)
				{
				k = i + (colwidth - i) / 2;
				out << MakeRightJustifiedString(FillWith((unsigned) i, '*'), k);
				out << FillWith(width - k, ' ');
				}
			else
				{
				NXS_ASSERT(justification != NxsTableCell::nojustify);
				out << MakeRightJustifiedString(FillWith((unsigned) i, '*'), colwidth);
				out << FillWith(NxsTable::colspacer, ' ');
				}
			break;
		
		case NxsTableCell::H :
				out << FillWith(colwidth, '-');
				out << FillWith(NxsTable::colspacer, ' ');
			break;
		
		case NxsTableCell::E :
			out << FillWith(width, ' ');
			break;

		default:
			out << "<error>";
		}
	return done;
	}
