#include "phycas/force_include.h"
#include <cmath>
#include "ncl/nxs_defs.hpp"
#include "ncl/misc/string_extensions.hpp"
#include "ncl/nxs_defs.hpp"
#include "ncl/output/nxs_ratio.hpp"
#include "ncl/misc/nxs_simple_stats.hpp"
#include "ncl/output/nxs_output_stream.hpp"
#include "ncl/output/nxs_table.hpp"
using std::vector;
using std::string;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::sprintf;
	using std::fabs;
#endif

unsigned NxsTable::rowincr		=	2;
unsigned NxsTable::colincr		=	2;
unsigned NxsTable::colspacer	=	2;

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `cells' and `activeCol' to NULL, `ncols' and `nrows' to 0, and calls Reset().
*/
NxsTable::NxsTable()
	{
	ncols		= 0;
	nrows		= 0;
	cells		= NULL;
	activeCol	= NULL;
	prnwidth		= 70; //@discuss with POL, had been in FactorySettings
	Reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Reset(`r', `c') after initializing `ncols' and `nrows' to 0 and `cells' and `activeCol' to NULL.
*/
NxsTable::NxsTable(
  unsigned r,	/* the initial number of rows in the table */
  unsigned c)	/* the initial number of columns in the table */
	{
	ncols		= 0;
	nrows		= 0;
	cells		= NULL;
	activeCol	= NULL;
	prnwidth		= 70; //@discuss with POL, had been in FactorySettings
	Reset(r, c);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Expands NxsTable by increasing the number of columns by `n'.
*/
void NxsTable::AddColumns(
  unsigned n)	/* the number of columns to add */
	{
	assert(n > 0);

	// Create a new activeCol array of the larger size and copy the first
	// ncols elements of the old activeCol array, initializing the remaining
	// elements to true.
	//
	unsigned k;
	bool *newActiveCol = new bool[ncols + n];
	for (k = 0; k < ncols; ++k)
		newActiveCol[k] = activeCol[k];
	for (k = 0; k < n; ++k)
		newActiveCol[ncols + k] = true;
	delete [] activeCol;
	activeCol = newActiveCol;
	newActiveCol = NULL;

	// Create nrows new rows (each ncols + n cells long) for the cells array, 
	// deleing the older, shorter row arrays.
	//
	for (unsigned i = 0; i < nrows; ++i)
		{
		NxsTableCell *newrow = new NxsTableCell[ncols + n];

		for (unsigned j = 0; j < ncols; ++j)
			newrow[j] = cells[i][j];

		delete [] cells[i];

		cells[i] = newrow;
		}

	ncols += n;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Expands NxsTable by increasing the number of rows.
*/
void NxsTable::AddRows(
  unsigned n)	/* the number of rows to add */
	{
	assert(n > 0);

	// Create a new cells array with n more rows than the current version
	//
	NxsTableCell **newcells = new NxsTableCell*[nrows + n];

	// Reassign the old rows to the new array, and create n new rows
	//
	for (unsigned i = 0; i < nrows; ++i)
		newcells[i] = cells[i];
	for (unsigned j = nrows; j < nrows + n; ++j)
		{
		newcells[j] = new NxsTableCell[ncols];
		}

	delete [] cells;
	cells = newcells;
	nrows += n;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First calls Clear to eliminate any existing data. Then allocates memory for the `cells' matrix and the `activeCol' 
|	array. Finally, called FactorySettings to initialize all the other data members to their default values. Called 
|	from the NxsTable constructor and can be called from outside to perform a Clear-and-fill. 
*/
void NxsTable::Reset(
  unsigned nr,	/* number of rows after (re)allocation */
  unsigned nc)	/* number of columns after (re)allocation */
	{
	Clear();
	
	nrows = nr;
	ncols = nc;
	
	activeCol = new bool[ncols];
	
	cells = new NxsTableCell*[nrows];
	for (unsigned i = 0; i < nrows; ++i)
		{
		cells[i] = new NxsTableCell[ncols];
		}
	
	FactorySettings();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sorts column `j' and fills in NxsTable::CIL, NxsTable::CIU, and NxsTable::CIP cells (if any).
*/
void NxsTable::CalcCI(
  unsigned j)	/* the column to be sorted */
	{
	// Sort cells from lowest value to highest (non-numeric cells first)
	//
	unsigned n = SortColumn(j);

	// Calculate k - the number to save from top and bottom
	//
	double d = (double)n * (double)nominal_CIpct / 100.0;
	unsigned k = (n - (unsigned)(d + 0.5)) / 2;
	d =  100.0 - 200.0 * (double)k / (double)n;
	unsigned actual = (unsigned)(d + 0.5);

	// Search for cell of type Upper and set its value if found
	//
	SetUpper(j, n, k);

	// Search for cell of type Lower and set its value if found
	//
	SetLower(j, k);

	// Search for cell of type Percent and set its value if found
	//
	SetPercent(j, actual);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds values of all cells in column `which_col' having type NxsTableCell::D, NxsTableCell::I, NxsTableCell::L, or 
|	NxsTableCell::R to the NxsSimpleStats object `stat'. If no value is specified for `which_col', it will be set to 
|	the current column (`curr_col').
*/
void NxsTable::CalcStats(
  NxsSimpleStats &stat,	/* the NxsSimpleStats object to receive values */
  unsigned which_col)	/* the column for which to obtain statistics */
	{
	if (which_col == UINT_MAX)
		which_col = curr_col;

	for (unsigned i = top; i < curr_row; ++i)
		{
		switch(cells[i][which_col].tag)
			{
			case NxsTableCell::D :
				stat += cells[i][which_col].d;
				break;
			case NxsTableCell::I :
				stat += static_cast<double>(cells[i][which_col].i);
				break;
			case NxsTableCell::L :
				stat += static_cast<double>(cells[i][which_col].l);
				break;
			case NxsTableCell::R :
				stat += cells[i][which_col].r;
				break;
			default:
				break; //default to avoid warning
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Makes sure all cells in column j have the same width. That common width is NxsTable::colspacer plus 
|~
|	o the user-specified width (if a width was specified for any cell in the column)
|	o the length of the	longest string (if user did not specified a width and the column contains at least one cell of
|	type NxsTableCell::S)
|	o the default column width (if no width was specified and the column contains no string cells)
|~
|	An exception is made for long strings. Strings longer than ABSOLUTE_MAXIMUM_COLUMN_WIDTH are truncated to 
|	ABSOLUTE_MAXIMUM_COLUMN_WIDTH. Also checks to see if any columns have had their justification set to something 
|	other than NxsTableCell::nojustify. If so, all columns are set to that justification; otherwise, all columns are 
|	set to the default, which is NxsTableCell::leftjustify.
*/
void NxsTable::CheckColumnWidth(
  unsigned j)	/* the column to check */
	{
	char tmp[34];
	unsigned lastRow = (bottom ? bottom : nrows);
	unsigned i, maxw = 0;
	unsigned iwidth, lwidth;
	unsigned slen, swidth = 0;
	unsigned fpwidth;
	NxsTableCell::TableColumnJustification which_just = NxsTableCell::leftjustify;

	for (i = 0; i < lastRow; ++i)
		{
		// See if justification has been specified for any cell in this column.
		// If so, remember the last one specified.
		//
		if (cells[i][j].justification != NxsTableCell::nojustify)
			which_just = cells[i][j].justification;

		// Accommodate largest specified cell width in the column
		//
		if (cells[i][j].width > maxw)
			maxw = cells[i][j].width;

		switch(cells[i][j].tag) 
			{
			case NxsTableCell::S :
				// Keep track of longest string in case user has specified
				// 0 for the width of all cells, in which case we want to
				// set maxw to min(swidth, ABSOLUTE_MAXIMUM_COLUMN_WIDTH) 
				//
				slen = (unsigned)cells[i][j].s.length();
				if (slen > swidth)
					swidth = slen;
				if (swidth > ABSOLUTE_MAXIMUM_COLUMN_WIDTH)
					swidth = ABSOLUTE_MAXIMUM_COLUMN_WIDTH;
				break;

			case NxsTableCell::D   :
			case NxsTableCell::SUM :
			case NxsTableCell::SS  :
			case NxsTableCell::M   :
			case NxsTableCell::V   :
			case NxsTableCell::SD  :
			case NxsTableCell::CV  :
			case NxsTableCell::CIU  :
			case NxsTableCell::CIL  :
				// Be sure to accommodate largest double in column
				//
				fpwidth = FloatStrLength(cells[i][j].d, cells[i][j].precision);
				if (fpwidth > maxw)
					maxw = fpwidth;
				break;
			case NxsTableCell::R   :
				// Be sure to accommodate largest ratio in column
				//
				fpwidth = FloatStrLength(cells[i][j].r.GetRatioAsDouble(), cells[i][j].precision);
				if (fpwidth > maxw)
					maxw = fpwidth;
				break;
			case NxsTableCell::I :
			case NxsTableCell::CIP :
				// Be sure to accommodate largest unsigned in column.
				//
				sprintf(tmp, "%d", cells[i][j].i);
				iwidth = (unsigned)strlen(tmp);
				if (iwidth > maxw)
					maxw = iwidth;
				break;
			case NxsTableCell::N :
			case NxsTableCell::L :
				// Be sure to accommodate largest long in column
				//
				sprintf(tmp, "%ld", cells[i][j].l);  
				lwidth = (unsigned)strlen(tmp);
				if (lwidth > maxw)
					maxw = lwidth;
				break;
			case NxsTableCell::E :
			case NxsTableCell::B :
			case NxsTableCell::H :
				break;
		}
	}

	// Set all cell widths to that of the widest cell, and set 
	// justification of each cell to the last one specified, or to
	// NxsTableCell::leftjustify if none were specified.
	//
	for (i = 0; i < lastRow; ++i)
		{
		cells[i][j].justification = which_just;

		// A column is at least as wide as the column spacer value
		//
		cells[i][j].width = NxsTable::colspacer;

		if (maxw == 0 && swidth == 0)
			{
			// Use the default column width because the user did not specify a width
			// for any cell in this column (maxw == 0) and there are no cells of type
			// S so we cannot use the width of the longest string as a guide (swidth == 0)
			//
			cells[i][j].width += NxsTableCell::def_width;
			}
		else if (maxw == 0)
			{
			// The user did not specify a width for any cell in this column, but the column
			// does contain string (type S) columns, so set all cells to the width of the 
			// longest string.
			//
			cells[i][j].width += swidth;
			}
		else
			{
			// The user specified a width for at least one cell in this column, so use 
			// that as the width of all cells in the column
			//
			cells[i][j].width += maxw;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the coefficient of variation for the observations in the current column. A NxsSimpleStats object (`stat')
|	is created and passed to the CalcStats method, which operates on the current column of the NxsTable. The value 
|	`stat'.CV() is returned.
*/
double NxsTable::CoeffVar()
	{
	NxsSimpleStats stat;
	CalcStats(stat);
	return stat.CV();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates an entire line of the NxsTable composed of only cells of type NxsTableCell::CV.
*/
void NxsTable::CoeffVarLine()
	{
	for (unsigned i = left; i <= right; ++i)
		this->AddCV();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recursive routine for sorting column `k' of the table. Note: should be called only through 
|	NxsTable::SortColumn. Originally written by Dmitri Zaykin for GDA.
*/
void NxsTable::QSort(
  unsigned k,		/* the column of the table to sort */
  unsigned dleft,	/* defines first row of sorted cells */
  unsigned dright)	/* defines last row of sorted cells */
	{
	unsigned i = dleft;
	unsigned j = dright;
	NxsTableCell x = cells[(dleft + dright) / 2][k];
	do
		{
		while (cells[i][k] < x  &&  i < right)
			++i;
		while (x < cells[j][k]  &&  j > left )
			--j;

		if (i <= j)
			{
			NxsTableCell y = cells[i][k];
			cells[i][k] = cells[j][k];
			cells[j][k] = y;
			++i;
			if (j)
				--j;
			}
		}
	while (i <= j);

	if (dleft < j)
		QSort(k, dleft, j);
	if (i < dright)
		QSort(k, i, dright);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Specifies initial settings for table variables, but does not affect the contents of cells. To "Clear-and-fill" the
|	table, clearing out all data, call Reset instead (Reset calls FactorySettings). 
*/
void NxsTable::FactorySettings()
	{
	nominal_CIpct	= 95;
	across			= true;
	//@temp discuss with pol prnwidth		= 70;
	curr_col		= 0;
	curr_row		= 0;
	left			= 0;
	right			= 0;
	top				= 0;
	bottom			= 0;
	hideFrom		= 0;
	hideTo			= 0;
	CItop			= 0;

	for (unsigned i = 0; i < ncols; ++i)
		{
		activeCol[i] = true;
		}

	for (unsigned i = 0; i < nrows; ++i) 
		{
		for (unsigned j = 0; j < ncols; ++j) 
			{
			cells[i][j].width = NxsTableCell::def_width;
			cells[i][j].precision = NxsTableCell::def_precision;
			}
		}
	rowincr = NxsTable::rowincr;
	colincr = NxsTable::colincr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Emulates the creation of a string representation of `d' using `p' to specify the numerical precision. Returns the 
|	length (in number of characters) that would be required by the string representation.
*/
unsigned NxsTable::FloatStrLength(
  double d,		/* the float */
  unsigned p)	/* the precision */
	{
	unsigned len = p;
	++p;	// add 1 for the decimal point
	if (d < 0.0)
		++len;		// add 1 for the minus sign
	long longd = (long)d;	// truncate decimals
	do {
		++len;
		longd /= 10;
	} while (longd > 0);
	return len;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes memory associated with the `cells' matrix and the `activeCol' array. Called from Reset, which is in turn 
|	called from the constructor.
*/
void NxsTable::Clear()
	{
	delete [] activeCol;
	activeCol = NULL;
	if (cells != NULL)
		{
		for (unsigned i = 0; i < nrows; ++i)
			delete [] cells[i];
		delete [] cells;
		cells = NULL;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Marks column `i' for hiding by setting `activeCol'[`i'] to false.
*/
void NxsTable::HideColumn(
  unsigned i)	/* the column to hide */
	{
	activeCol[i] = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates an entire line of the NxsTable composed of only cells of type NxsTableCell::H.
*/
void NxsTable::HyphenLine()
	{
	for (unsigned i = left; i <= right; ++i)
		this->AddHyphens();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates an entire line of the NxsTable composed of only cells of type NxsTableCell::CIL.
*/
void NxsTable::LowerLine()
	{
	for (unsigned i = left; i <= right; ++i)
		this->AddCIL();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the mean for the current column. A NxsSimpleStats object (`stat') is created and passed to the CalcStats 
|	method, which operates on the current column of the NxsTable. The value `stat'.Mean() is returned.
*/
double NxsTable::Mean()
	{
	NxsSimpleStats stat;
	CalcStats(stat);
	return stat.Mean();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates an entire line of the NxsTable composed of only cells of type NxsTableCell::M.
*/
void NxsTable::MeanLine()
	{
	for (unsigned i = left; i <= right; ++i)
		this->AddMean();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Moves cursor (`curr_col') to the right one cell. If that places it beyond the rightmost cell in the NxsTable body, 
|	the cursor is set to the first cell in the next row down.
*/
void NxsTable::NextAcross()
	{
	++curr_col;
	if (right && curr_col >= right)
		{
		curr_col = 0;
		NextDown();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Moves cursor to the next cell, which is to the right if `across' is true (down if `across' is false). Wrapping 
|	occurs if right (`bottom') is exceeded.
*/
void NxsTable::NextCell()
	{
	if (across)
		NextAcross();
	else
		NextDown();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `curr_row' to the first cell in the next row down.
*/
void NxsTable::NextDown()
	{
	++curr_row;
	if (bottom && curr_row >= bottom)
		{
		//@pol isn't this a situation where we should use an assert?
		curr_row = 0;
		NextAcross();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates an entire line of the NxsTable composed of only cells of type NxsTableCell::CIP.
*/
void NxsTable::PercentLine()
	{
	for (unsigned i = left; i <= right; ++i)
		this->AddCIP();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes sample size for the current column. A NxsSimpleStats object (`stat') is created and passed to the CalcStats
|	method, which operates on the current column of the NxsTable. The value `stat'.N() is returned.
*/
long NxsTable::SampleSize()
	{
	NxsSimpleStats stat;
	CalcStats(stat);
	return stat.N();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the lower bound for a cell of type NxsTableCell::CIL in the NxsTable object. Assumes column has already been 
|	sorted from low to high. Starting with the first cell of type double encountered, retrieves the kth cell's value - 
|	this is the lower bound `d'. Continuing on down the column until a cell having the NxsTableCell::CIL tag is found,
|	this cell's value is then set to `d'.
*/
void NxsTable::SetLower(
  unsigned col,	/* the column for which a lower CI bound is needed */
  unsigned k)	/* the kth valid double in col is the lower bound */
	{
	unsigned i = top;
	unsigned m = 0;
	while (cells[i][col].tag == NxsTableCell::B)
		++i;
	while (m < k)
		{
		++i; 
		++m;
		}
	double lowerbnd = cells[i][col].d;

	// Find the first cell of type Lower and set its value
	//
	while (i < bottom && cells[i][col].tag != NxsTableCell::CIL)
		++i;
	if (i < bottom)
		cells[i][col].d = lowerbnd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Finds the first cell of type NxsTableCell::CIP and sets its value.
*/
void NxsTable::SetPercent(
  unsigned col,		/* the column containing the Percent cell */
  unsigned actual)	/* the value to be stored in the Percent cell */
	{
	unsigned i = top;

	while (i < bottom && cells[i][col].tag != NxsTableCell::CIP) 
		++i;
	if (i < bottom)
		cells[i][col].i = (int) actual;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the floating-point precision for the entire table. Every NxsTableCell::D type cell will have `p' digits 
|	following the decimal point.
*/
void NxsTable::SetPrecision(
  unsigned p)	/* the precision to be applied to every NxsTableCell::D cell in the table */
	{
	for (unsigned i = 0; i < nrows; ++i) 
		{
		for (unsigned j = 0; j < ncols; ++j)
			{
			if (cells[i][j].tag == NxsTableCell::D)
				cells[i][j].precision = p;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the width (in characters) of the output. The table object will try to fit within this width by splitting the
|	table into several smaller tables, each containing a subset of columns.
*/
void NxsTable::SetPrnWidth(
  unsigned w)	/* the new print width (in characters) */
	{
	prnwidth = w;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the upper bound for a cell of type NxsTableCell::CIU in the NxsTable object. Assumes column has already been 
|	sorted from low to high. Starting with the first cell of type NxsTableCell::D encountered, retrieves the 
|	(`n' - `k')th cell's value - this is the upper bound `d'. Continuing on down the column until a cell having the 
|	NxsTableCell::CIU tag is found, this cell's value is then set to `d'.
*/
void NxsTable::SetUpper(
  unsigned col,	/* the column for which a lower CI bound is needed */
  unsigned n,	/* the number of entries in the column */
  unsigned k)	/* the k'th valid double from the end in col is the upper bound */
	{
	unsigned i = top;
	unsigned m = 0;
	while (cells[i][col].tag == NxsTableCell::B)
		++i;
	while (m < n-k-1)
		{
		++i;
		++m;
		}
	double upperbnd = cells[i][col].d;

	// Find the first cell of type NxsTableCell::CIU and set its value
	//
	while (i < bottom && cells[i][col].tag != NxsTableCell::CIU)
		++i;
	if (i < bottom)
		cells[i][col].d = upperbnd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the cell width for the entire table. Every cell will be `w' characters in width.
*/
void NxsTable::SetWidth(
  unsigned w)	/* the width to be applied to every cell in the table */
	{
	for (unsigned i = 0; i < nrows; ++i) 
		{
		for (unsigned j = 0; j < ncols; ++j)
			{
			cells[i][j].width = w;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets both the column width and numerical precision for the entire table. Every cell will be `w' characters wide 
|	and NxsTableCell::D type cells will have `p' digits following the decimal point.
*/
void NxsTable::SetWidthAndPrecision(
  unsigned w,	/* the width to be applied to every cell in the table */
  unsigned p)	/* the precision to be applied to every NxsTableCell::D cell in the table */
	{
	for (unsigned i = 0; i < nrows; ++i) 
		{
		for (unsigned j = 0; j < ncols; ++j)
			{
			cells[i][j].width = w;
			if (cells[i][j].tag == NxsTableCell::D)
				cells[i][j].precision = p;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sorts column `k', returning the number of cells having numeric types (e.g. NxsTableCell::D, NxsTableCell::I,
|	NxsTableCell::L, or NxsTableCell::R). Returns without doing anything if `CItop' has not been previously set.
|	Before returning, converts all numeric cells to type NxsTableCell::D and all non-numeric cells to type 
|	NxsTableCell::B so they will be displayed with the default number of asterisks.
*/
unsigned NxsTable::SortColumn(
  unsigned k)	/* the column of the NxsTable to sort */
	{
	if (CItop == UINT_MAX) 
		return 0;

	unsigned i;
	unsigned n = 0;

	// Quick sort
	//
	QSort(k, top, CItop - 1);

	// Convert all numeric cells to type NxsTableCell::D, all others to NxsTableCell::B
	//
	double tmp;
	for (i = top; i < CItop; ++i)
		{
		if (cells[i][k].tag == NxsTableCell::I)
			{
			tmp = (double)cells[i][k].i;
			cells[i][k] = tmp;
			}
		else if (cells[i][k].tag == NxsTableCell::L)
			{
			tmp = (double)cells[i][k].l;
			cells[i][k] = tmp;
			}
		else if (cells[i][k].tag == NxsTableCell::R)
			{
			tmp = cells[i][k].r.GetRatioAsDouble();
			cells[i][k] = tmp;
			}

		if (cells[i][k].tag == NxsTableCell::D)
			++n;
		else
			{
			NxsTableCell &c = cells[i][k];
			c.i = (int) NxsTableCell::def_asterisks;
			c.tag = NxsTableCell::B;
			}
		}

	return n;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the standard deviation of the observations in the current column. A NxsSimpleStats object (`stat') is 
|	created and passed to the CalcStats method, which operates on the current column of the NxsTable. The value 
|	`stat'.StdDev() is returned.
*/
double NxsTable::StdDev()
	{
	NxsSimpleStats stat;
	CalcStats(stat);
	return stat.StdDev();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates an entire line of the NxsTable composed of only cells of type NxsTableCell::SD.
*/
void NxsTable::StdDevLine()
	{
	for (unsigned i = left; i <= right; ++i)
		this->AddStdDev();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the sum of all the cells for the current column. A NxsSimpleStats object (`stat)' is created and passed to
|	the CalcStats method, which operates on the current column of the NxsTable. The value `stat'.Sum() is returned.
*/
double NxsTable::Sum()
	{
	NxsSimpleStats stat;
	CalcStats(stat);
	return stat.Sum();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates an entire line of the NxsTable composed of only cells of type NxsTableCell::SUM.
*/
void NxsTable::SumLine()
	{
	for (unsigned i = left; i <= right; ++i)
		this->AddSum();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the sum of squares for the current column. A NxsSimpleStats object (`stat') is created and passed to the 
|	CalcStats method, which operates on the current column of the NxsTable. The value `stat'.SumSq() is returned.
*/
double NxsTable::SumSq()
	{
	NxsSimpleStats stat;
	CalcStats(stat);
	return stat.SumSq();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates an entire line of the NxsTable composed of only cells of type NxsTableCell::SS.
*/
void NxsTable::SumSqLine()
	{
	for (unsigned i = left; i <= right; ++i)
		this->AddSumOfSquares();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates an entire line of the NxsTable composed of only cells of type NxsTableCell::CIU.
*/
void NxsTable::UpperLine()
	{
	for (unsigned i = left; i <= right; ++i)
		this->AddCIU();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the variance of the observations in the current column. A NxsSimpleStats object (`stat') is created and 
|	passed to the CalcStats method, which operates on the current column of the NxsTable. The value `stat'.Variance() 
|	is returned.
*/
double NxsTable::Variance()
	{
	NxsSimpleStats stat;
	CalcStats(stat);
	return stat.Variance();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates an entire line of the NxsTable composed of only cells of type NxsTableCell::V.
*/
void NxsTable::VarianceLine()
	{
	for (unsigned i = left; i <= right; ++i)
		this->AddVariance();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets justification data member of current column to NxsTableCell::leftjustify. When CheckColumnWidth is called, 
|	all cells in this column will be set to NxsTableCell::leftjustify if any of them are thus set. In case of 
|	conflicts, all cells are set to the the last justification specified.
*/
void NxsTable::LeftJustifyColumn()
	{
	NxsTableCell &c = cells[curr_row][curr_col];
	c.justification = NxsTableCell::leftjustify;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets justification data member of current column to NxsTableCell::centerjustify. When CheckColumnWidth is called, 
|	all cells in this column will be set to NxsTableCell::centerjustify if any of them are thus set. In case of 
|	conflicts, all cells are set to the the last justification specified.
*/
void NxsTable::CenterColumn()
	{
	NxsTableCell &c = cells[curr_row][curr_col];
	c.justification = NxsTableCell::centerjustify;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets justification data member of current column to NxsTableCell::rightjustify. When CheckColumnWidth is called, 
|	all cells in this column will be set to NxsTableCell::rightjustify if any of them are thus set. In case of 
\	conflicts, all cells are set to the the last justification specified.
*/
void NxsTable::RightJustifyColumn()
	{
	NxsTableCell &c = cells[curr_row][curr_col];
	c.justification = NxsTableCell::rightjustify;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Specifies the left boundary of the NxsTable body. All columns added before this point are considered row headers 
|	and will be repeated for each output page.
*/
void NxsTable::SetLeftMargin()
	{
	left = curr_col;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Specifies the right boundary of the NxsTable body. This fixes the number of columns in the NxsTable body. Cells 
|	added to the right of this point will wrap around and end up instead on the next row down.
*/
void NxsTable::SetRightMargin()
	{
	right = curr_col;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Specifies the upper boundary of the NxsTable body. All rows added before this point are considered column headers 
|	and will be repeated for each output page.
*/
void NxsTable::SetTopMargin()
	{
	top = curr_row;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Specifies the lower boundary of the NxsTable body. This fixes the number of rows in the NxsTable body. This function
|	should only be used when all cells have been added. Use of this function is necessary to prevent rows of all zeros 
|	from appearing at the bottom of the NxsTable when printed. Rows of zeros exist because the NxsTable may have 
|	automatically resized itself recently and added an extra row or two to the bottom.
*/
void NxsTable::SetBottomMargin()
	{
	bottom = curr_row;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Instructs NxsTable to begin hiding columns after this point and continue hiding until HideTo method is called.
*/
void NxsTable::HideFrom()
	{
	hideFrom = curr_row;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Instructs NxsTable to stop hiding columns at this point. The HideFrom method should have already been called.
*/
void NxsTable::HideTo()
	{
	hideTo = curr_row;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Call before accessing an element of the cells matrix. Adds row or column if necessary to cells matrix using 
|	AddRows or AddColumns method, respectively.
*/
void NxsTable::CheckExpand()
	{
	// Check whether we are in bounds; tried putting this check in the functions
	// AddColumns and AddRows, but that way you always end up with one more column
	// or row than you need because of the call to NextCell
	//
	if (curr_col >= ncols)
		AddColumns();
	if (curr_row >= nrows)
		AddRows();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::S into the NxsTable.
*/
void NxsTable::AddString(
  string s,	/* the string to be inserted */
  unsigned w)	/* the column width */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c = s;
	if (w > 0)
		c.width = w;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::D into the NxsTable.
*/
void NxsTable::AddDouble(
  double d,		/* the double value to be inserted */
  unsigned w,	/* the column width */
  unsigned p)	/* the numerical precision */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c = d;
	if (w > 0)
		c.width = w;
	c.precision = p;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::D into the NxsTable if the supplied array `d' has only a single element. 
|	Otherwise, inserts a cell of type NxsTableCell::S containing a parenthetical list of the double values composing
|	`d', with each value separated by a space character (e.g. "(0.4 0.3 0.2 0.1)").
*/
void NxsTable::AddDoubleArr(
  double *d,	/* pointer to the array of doubles to be inserted */
  unsigned len,	/* number of doubles to be inserted */
  unsigned w,	/* the column width */
  unsigned p)	/* the numerical precision */
	{
	if (len == 1)
		AddDouble(*d, w, p);
	else
		{
		CheckExpand();
		NxsTableCell &c = cells[curr_row][curr_col];
		char temp[100];
		temp[99] ='\0';
		string s;
		s = '(';
		for (unsigned i = 0; i < len; ++i)
			{
			if (fabs(d[i]) > 1.0e8 || fabs(d[i])<1.0e-6)
				sprintf(temp,"%e", d[i]);
			else
				sprintf(temp,"%f", d[i]);
			assert(!temp[99]);
			s << temp  << ' ';
			}
		s << ')';
		c = s;
		if (w > 0)
			c.width = w;
		c.precision = p;
		NextCell();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Serves as a front end for the function AddDoubleArr, for use when the array to be displayed is in the form of an
|	STL vector of doubles. Creates an array `d' containing the values in the vector `v', then passes `d' on to 
|	AddDoubleArr for processing.
*/
void NxsTable::AddDoubleVec(
  const vector<double> &v,	/* reference to the vector of doubles to be inserted */
  unsigned w,				/* the column width */
  unsigned p)				/* the numerical precision */
	{
	if (v.size() > 0)
		{
		unsigned s = (unsigned)v.size();
		double *d = new double [s];
		for (unsigned i = 0; i < s; ++i)
			d[i] = v[i];
		AddDoubleArr(d,s,w,p);
		delete [] d;
		}
	else
		AddString("()",w);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::I into the NxsTable if the supplied array `d' has only a single element. 
|	Otherwise, inserts a cell of type NxsTableCell::S containing a parenthetical list of the int values composing
|	`d', with each value separated by a space character (e.g. "(45 33 256 112)").
*/
void NxsTable::AddIntArr(
  int *d,		/* the int value to be inserted */
  unsigned len,	/* number of int values to be inserted */
  unsigned w)	/* the column width */
	{
	if (len == 1)
		AddInt(*d, w);
	else
		{
		CheckExpand();
		NxsTableCell &c = cells[curr_row][curr_col];
		char temp[100];
		temp[99] ='\0';
		string s;
		s = '(';
		for (unsigned i = 0; i < len; ++i)
			{
			sprintf(temp,"%d",d[i]);
			s << temp << ' ';
			}
		s << ')';
		c = s;
		if (w > 0)
			c.width = w;
		NextCell();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Serves as a front end for the function AddIntArr, for use when the array to be displayed is in the form of an
|	STL vector of ints. Creates an array `d' containing the values in the vector `v', then passes `d' on to 
|	AddIntArr for processing.
*/
void NxsTable::AddIntVec(
  const vector<int> &v,	/* reference to the vector of int values to be inserted */
  unsigned w)			/* the column width */
	{
	if (v.size() > 0)
		{
		unsigned s = (unsigned)v.size();
		int *d = new int[s];
		for (unsigned i = 0; i < s; ++i)
			d[i] = v[i];
		AddIntArr(d,s,w);
		delete [] d;
		}
	else
		AddString("()",w);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::I into the NxsTable
*/
void NxsTable::AddInt(
  int i,		/* the unsigned value to be inserted */
  unsigned w)	/* the column width */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c = i;
	if (w > 0)
		c.width = w;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::L into the NxsTable.
*/
void NxsTable::AddLong(
  long l,		/* the long value to be inserted */
  unsigned w)	/* the column width */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c = l;
	if (w > 0)
		c.width = w;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::R into the NxsTable.
*/
void NxsTable::AddRatio(
  NxsRatio r,	/* the ratio value to be inserted */
  unsigned w,	/* the column width */
  unsigned p)	/* the numerical precision */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c = r;
	if (w > 0)
		c.width = w;
	c.precision = p;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::B into the NxsTable.
*/
void NxsTable::AddAsterisks(
  unsigned num,	/* the number of asterisks */
  unsigned w)	/* the column width */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	c.i = (int) num;
	c.tag = NxsTableCell::B;
	if (w > 0)
		c.width = w;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::H into the NxsTable.
*/
void NxsTable::AddHyphens(
  unsigned num,	/* the number of hyphens */
  unsigned w)	/* the column width */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	c.i = (int) num;
	c.tag = NxsTableCell::H;
	if (w > 0)
		c.width = w;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::E into the NxsTable.
*/
void NxsTable::AddEmpty(
  unsigned w)	/* the column width */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	c.tag = NxsTableCell::E;
	if (w > 0)
		c.width = w;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::N into the NxsTable.
*/
void NxsTable::AddSampleSize(
  unsigned w)	/* the column width */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	c.l = SampleSize();
	c.tag = NxsTableCell::N;
	if (w > 0)
		c.width = w;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::SUM into the NxsTable.
*/
void NxsTable::AddSum(
  unsigned w,	/* the column width */
  unsigned p)	/* the numerical precision */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	c.d = Sum();
	c.tag = NxsTableCell::SUM;
	if (w > 0)
		c.width = w;
	c.precision = p;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::SS into the NxsTable.
*/
void NxsTable::AddSumOfSquares(
  unsigned w,	/* the column width */
  unsigned p)	/* the numerical precision */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	c.d = SumSq();
	c.tag = NxsTableCell::SS;
	if (w > 0)
		c.width = w;
	c.precision = p;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::M into the NxsTable.
*/
void NxsTable::AddMean(
  unsigned w,	/* the column width */
  unsigned p)	/* the numerical precision */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	if (SampleSize() > 0L)
		{
		c.tag = NxsTableCell::M;
		c.d = Mean();
		}
	else
		{
		c.tag = NxsTableCell::B;
		c.i = (int) NxsTableCell::def_asterisks;
		}
	if (w > 0)
		c.width = w;
	c.precision = p;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::V into the NxsTable.
*/
void NxsTable::AddVariance(
  unsigned w,	/* the column width */
  unsigned p)	/* the numerical precision */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	if (SampleSize() > 1L)
		{
		c.tag = NxsTableCell::V;
		c.d = Variance();
		}
	else
		{
		c.tag = NxsTableCell::B;
		c.i = (int) NxsTableCell::def_asterisks;
		}
	if (w > 0)
		c.width = w;
	c.precision = p;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::SD into the NxsTable.
*/
void NxsTable::AddStdDev(
  unsigned w,	/* the column width */
  unsigned p)	/* the numerical precision */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	if (SampleSize() > 1L)
		{
		c.tag = NxsTableCell::SD;
		c.d = StdDev();
		}
	else
		{
		c.tag = NxsTableCell::B;
		c.i = (int) NxsTableCell::def_asterisks;
		}
	if (w > 0)
		c.width = w;
	c.precision = p;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::CV into the NxsTable.
*/
void NxsTable::AddCV(
  unsigned w,	/* the column width */
  unsigned p)	/* the numerical precision */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	if (SampleSize() > 1L)
		{
		c.tag = NxsTableCell::CV;
		c.d = CoeffVar();
		}
	else
		{
		c.tag = NxsTableCell::B;
		c.i = (int) NxsTableCell::def_asterisks;
		}
	if (w > 0)
		c.width = w;
	c.precision = p;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::CIU into the NxsTable.
*/
void NxsTable::AddCIU(
  unsigned w,	/* the column width */
  unsigned p)	/* the numerical precision */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	c.tag = NxsTableCell::CIU;

	// CItop used when outputting to determine when it is safe to
	// do the sorting necessary to find confidence intervals
	//
	if (CItop == 0 || (curr_row < CItop))
		CItop = curr_row;

	if (w > 0)
		c.width = w;
	c.precision = p;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::CIL into the NxsTable.
*/
void NxsTable::AddCIL(
  unsigned w,	/* the column width */
  unsigned p)	/* the numerical precision */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	c.tag = NxsTableCell::CIL;

	// CItop used when outputting to determine when it is safe to
	// do the sorting necessary to find confidence intervals
	//
	if (CItop == 0 || (curr_row < CItop))
		CItop = curr_row;

	if (w > 0)
		c.width = w;
	c.precision = p;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts a cell of type NxsTableCell::CIP into the NxsTable.
*/
void NxsTable::AddCIP(
  unsigned w,	/* the column width */
  unsigned p)	/* the numerical precision */
	{
	CheckExpand();
	NxsTableCell &c = cells[curr_row][curr_col];
	c.FactorySettings();
	c.tag = NxsTableCell::CIP;

	// CItop used when outputting to determine when it is safe to
	// do the sorting necessary to find confidence intervals
	//
	if (CItop == 0 || (curr_row < CItop))
		CItop = curr_row;

	if (w > 0)
		c.width = w;
	c.precision = p;
	NextCell();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The only means of outputting the formatted table object.
*/
void NxsTable::Show(
  NxsOutputStream & out)	/* the output stream object to which the formatted table will be sent */
	{
	assert(hideFrom <= nrows);
	assert(hideTo <= nrows);
	assert(hideTo >= hideFrom);

	unsigned	lastRow	= (bottom > 0 ? bottom : nrows);
	unsigned	lastCol	= (right > 0 ? right : ncols);
	bool		hiding	= (hideTo - hideFrom > 0 ? true : false);
	bool		done	= false;
	unsigned	fst;	// fst is first column to be printed on current page
	unsigned	lst;	// lst is last column to be printed on current page
	unsigned	cumwid, cumlabelwid, i, j;

	// Ensure all cells in any one column have the same width. Hereafter, we can determine 
	// the column width by looking at the width of any cell in the column.
	//
	for (j = 0; j < lastCol; ++j)
		{
		CheckColumnWidth(j);
		}

	// Obtain cumulative width of label columns
	//
	cumlabelwid = 0;
	for (j = 0; j < left; ++j)
		{
		if (!activeCol[j]) 
			continue;
		cumlabelwid += cells[0][j].width /* NxsTableCell::width */;
		}

	// If a page with only one non-label column is too wide, throw an exception.
	//
	if (cumlabelwid + cells[0][left].width >= prnwidth)
		{
		throw NxsTable::NxsX_InsufficientWidth();
		}

	fst = left;
	while (!done)
		{
		// Let cumwid accumulate (active) column widths until cumwid is greater than prnwidth.
		//   fst is first column to be printed this page
		//   lst is last column to be printed this page
		//
		cumwid = cumlabelwid;
		for (j = fst; j < lastCol && cumwid <= prnwidth; ++j)
			{
			if (!activeCol[j]) 
				continue;
			cumwid += cells[0][j].width;
			}

		// At this point, either cumwid has exceeded prnwidth or all active columns have been included
		//
		lst = (cumwid > prnwidth ? j - 1 : j);

		// Print the column header rows.
		// Note that it is assumed that none of these are too long, so no wrapping is done
		//
		for (i = 0; i < top; ++i) 
			{
			// Skip hidden rows
			//
			if (hiding && i >= hideFrom && i < hideTo)
				continue;

			out << ncl::endl;

			// Print out the label columns
			//
			for (j = 0; j < lst; ++j) 
				{
				if (!activeCol[j])  
					continue;

				// If this is not the first page, need to account for the gap between left and fst
				//
				if (fst != left && j >= left && j < fst)
					continue;

				cells[i][j].Show(out);
				}
			}

		// Print non-header rows for all columns appearing on this page
		//
		for (i = top; i < lastRow; ++i) 
			{
			// Skip hidden rows
			//
			if (hiding && i >= hideFrom && i < hideTo)
				continue;

			// Calculate confidence intervals if requested
			//
			if (CItop > 0 && i == CItop) 
				{
				// Sort columns and set values for CI-related cells
				//
				for (j = fst; j < lst; ++j)
					CalcCI(j);
				}

			out << ncl::endl;

			// Print row label cells for row i
			//
			for (j = 0; j < left; ++j) 
				{
				if (!activeCol[j]) 
					continue;

				// If Show returns false, it means not all of the contents of the cell were able to be 
				// printed in the allotted space. In this case, keep cycling until all of the cell's contents
				// are able to be printed.
				//
				while (!cells[i][j].Show(out))
					{
					out << ncl::endl;

					// Need to spit out empty cells for the label cells to the left of this one
					// before spitting out the next part of this label cell's text.
					//
					for (unsigned jj = 0; jj < j; ++jj)
						{
						if (!activeCol[jj]) 
							continue;

						// If this is not the first page, need to account for the gap between left and fst
						//
						if (fst != left && jj >= left && jj < fst)
							continue;

						cells[i][jj].ShowEmpty(out);
						}
					}
				}

			// Print non-label cells for row i
			//
			for (j = fst; j < lst; ++j) 
				{
				if (!activeCol[j]) 
					continue;

				// If Show returns false, it means not all of the contents of the cell were able to be 
				// printed in the allotted space. In this case, keep cycling until all of the cell's contents
				// are able to be printed.
				//
				while (!cells[i][j].Show(out))
					{
					out << ncl::endl;

					// Need to spit out empty cells for all cells to the left of this one
					// before spitting out the next part of this cell's text.
					//
					for (unsigned jj = 0; jj < j; ++jj)
						{
						if (!activeCol[jj]) 
							continue;

						// If this is not the first page, need to account for the gap between left and fst
						//
						if (fst != left && jj >= left && jj < fst)
							continue;

						cells[i][jj].ShowEmpty(out);
						}
					}
				}	// j-loop over non-label cells
			}	// i-loop over non-header rows

		// Set fst and lst for start of next round
		//
		if (lst >= lastCol)
			done = true;
		else
			fst = lst;

		out << '\n' << ncl::endl;
	} // while !done
}

#if 0 //HIDE_PROTECTION_FUNCS
	NxsTable::NxsTable(const NxsTable &)
		{
		assert(0);
		throw BadCopyConstrXcp("NxsTable");
		}
	NxsTable &NxsTable::operator=(const NxsTable &)
		{
		assert(0);
		throw BadEqOpXcp("NxsTable");
		}

#endif
