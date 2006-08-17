#ifndef NCL_NXSTABLE_H
#define NCL_NXSTABLE_H

#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/output/nxs_table_cell.hpp"

class NxsSimpleStats;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates an output table. Allows tabular output to be formatted according to a specified output width, column 
|	width and numerical precision. If output width is not enough to accommodate entire table, table will be split into 
|	two or more separate tables. Each separate table copies initial columns designated row headers and initial rows 
|	designated column headers. Below is example code showing how to create a NxsTable object, fill it with tabular 
|	output, then send it to an output stream object. The two-dimensional array of doubles called data has already been
|	initialized at this point with the values to be printed. Values in the data array considered invalid are indicated
|	by the value 9999, which is detected in the code below, resulting in the NxsTableCell::Bad object being put in the 
|	table instead of the actual array value. The AddHyphens method causes a NxsTableCell::Hyphens object to occupy the 
|	next cell, and the default behavior is for hyphens to fill up the entire cell width except for the extreme left
|	character. To specify that 5 hyphens are to be used rather than enough to span the column, use AddHyphens(5).
|>
|	NxsTable t;
|	t.AddString("Locus");
|	t.SetLeftMargin();
|	t.AddString("f");
|	t.AddString("F");
|	t.AddString("Theta");
|	t.AddString("S3");
|	t.AddString("S2");
|	t.AddString("S1");
|	t.AddString("C");
|	t.SetRightMargin();
|	for (unsigned i = 0; i < 8; ++i)
|	  t.AddHyphens();
|	t.SetTopMargin();
|	for (i = 0; i < 5; ++i) 
|	   {
|	   std::string name = "locus-";
|	   name << (i+1);
|	   t.AddString(name);
|	   for (j = 0; j < 7; ++j) 
|	      {
|	      if (data[i][j] > 999)
|	         t.AddAsterisks();
|	      else
|	         t.AddDouble(data[i][j]);
|	      }
|	   }
|	t.SetBottomMargin();
|	t.SetPrnWidth(70);
|	t.SetWidthAndPrecision(12, 6);
|	t.Show(outStream);
|>
|	The methods SetTopMargin, SetBottomMargin, SetLeftMargin, and SetRightMargin determine the boundaries of the table
|	body. Any columns added to the table before SetLeftMargin become row headings and any rows added to the table before
|	SetTopMargin become column headings. Specifying SetRightMargin and SetBottomMargin determine the total number of 
|	body columns and rows, respectively. Here is the resulting table as it looks printed:
|>
|	       Locus           f           F       Theta          S3
|	 ----------- ----------- ----------- ----------- -----------
|	     locus-1         ***         ***         ***   88.000000
|	     locus-2    0.250514    0.302529    0.069401   44.000000
|	     locus-3   -0.034807   -0.030052    0.004595   27.333333
|	     locus-4    0.012751    0.036738    0.024297   21.500000
|	     locus-5         ***         ***         ***   88.000000
|	
|	
|	       Locus          S2          S1           C
|	 ----------- ----------- ----------- -----------
|	     locus-1   88.000000   88.000000   88.000000
|	     locus-2   42.000000   39.118056   38.409091
|	     locus-3   21.333333   16.555952   15.829268
|	     locus-4   15.000000    9.262500    8.215116
|	     locus-5   88.000000   88.000000   88.000000
|>
*/
class NxsTable
	{
	public:
		class	NxsX_InsufficientWidth {};

	public:
				NxsTable();
				NxsTable(unsigned r, unsigned c);
				~NxsTable();

		void 	FactorySettings();

		void 	Reset(unsigned nr = 5, unsigned nc = 5);
		void 	SetNominalCIPercent(unsigned p);
		void 	SetPrnWidth(unsigned w);

		void 	SetWidth(unsigned w);
		void 	SetPrecision(unsigned p);
		void 	SetWidthAndPrecision(unsigned w, unsigned p);
		void 	HideColumn(unsigned i);
		void 	CalcCI(unsigned i);
		void 	CalcStats(NxsSimpleStats &, unsigned which_col = UINT_MAX);

		void 	UpperLine();
		void 	LowerLine();
		void 	PercentLine();
		void 	HyphenLine();
		void 	SumLine();
		void 	SumSqLine();
		void 	MeanLine();
		void 	VarianceLine();
		void 	StdDevLine();
		void 	CoeffVarLine();

		long 	SampleSize();
		double 	Sum();
		double 	SumSq();
		double 	Mean();
		double 	Variance();
		double 	StdDev();
		double 	CoeffVar();

		
		void 	HideFrom();
		void 	HideTo();

		void 	SetLeftMargin();
		void 	SetRightMargin();
		void 	SetTopMargin();
		void 	SetBottomMargin();

		void 	LeftJustifyColumn();
		void 	CenterColumn();
		void 	RightJustifyColumn();

		void 	CheckExpand();
		void 	AddString(std::string s, unsigned w = NxsTableCell::def_width);
		void 	AddUInt(unsigned i, unsigned w = NxsTableCell::def_width);
		void 	AddInt(int i, unsigned w = NxsTableCell::def_width);
		void 	AddIntArr(int *d, unsigned len, unsigned w = NxsTableCell::def_width);
		void 	AddIntVec(const std::vector<int> &v,unsigned w = NxsTableCell::def_width );
		void 	AddLong(long l, unsigned w = NxsTableCell::def_width);
		void 	AddDouble(double d, unsigned w = NxsTableCell::def_width, unsigned p = NxsTableCell::def_precision);
		void 	AddDoubleArr(double *d, unsigned len, unsigned w = NxsTableCell::def_width, unsigned p = NxsTableCell::def_precision);
		void 	AddDoubleVec(const std::vector<double> &v,unsigned w = NxsTableCell::def_width , unsigned p  = NxsTableCell::def_precision );
		void 	AddRatio(NxsRatio r, unsigned w = NxsTableCell::def_width, unsigned p = NxsTableCell::def_precision);
		void 	AddAsterisks(unsigned num = NxsTableCell::def_asterisks, unsigned w = 0);
		void 	AddHyphens(unsigned num = 0, unsigned w = 0);
		void 	AddEmpty(unsigned w = 0);
		void 	AddSampleSize(unsigned w = NxsTableCell::def_width);
		void 	AddSum(unsigned w = NxsTableCell::def_width, unsigned p = NxsTableCell::def_precision);
		void 	AddSumOfSquares(unsigned w = NxsTableCell::def_width, unsigned p = NxsTableCell::def_precision);
		void 	AddMean(unsigned w = NxsTableCell::def_width, unsigned p = NxsTableCell::def_precision);
		void 	AddVariance(unsigned w = NxsTableCell::def_width, unsigned p = NxsTableCell::def_precision);
		void 	AddStdDev(unsigned w = NxsTableCell::def_width, unsigned p = NxsTableCell::def_precision);
		void 	AddCV(unsigned w = NxsTableCell::def_width, unsigned p = NxsTableCell::def_precision);
		void 	AddCIU(unsigned w = NxsTableCell::def_width, unsigned p = NxsTableCell::def_precision);
		void 	AddCIL(unsigned w = NxsTableCell::def_width, unsigned p = NxsTableCell::def_precision);
		void 	AddCIP(unsigned w = NxsTableCell::def_width, unsigned p = NxsTableCell::def_precision);

		void 	Show(NxsOutputStream &out);	//@POL Mark, wasn't sure why this was private before. Any reason?

		//friend	ostream	&operator<<(ostream &, NxsTable &);
		//friend	NxsTable &operator<<(NxsTable &, NxsTableCell &);

	public:
		STATIC_DATA unsigned	rowincr;			/* the number of rows to add when NxsTable resized */
		STATIC_DATA unsigned	colincr;			/* the number of columns to add when NxsTable resized */
		STATIC_DATA unsigned	colspacer;			/* minimum number of characters separating columns */

	
		unsigned 		left;				/* the column forming left side of table body; columns before `left' are considered headers and are repeated on each page */
		unsigned 		top;				/* the row forming top side of table body; rows above `top' are considered headers and are repeated on each page */
		unsigned 		right;				/* the column forming right side of table body; any columns added after column `right' will go in the next row */
		unsigned 		bottom;				/* the row forming bottom side of table body */
		unsigned 		hideFrom;			/* the first row in group of hidden rows */
		unsigned 		hideTo;				/* the last row in group of hidden rows */
		unsigned 		CItop;				/* first row with cell of type Lower, Upper, or Percent */
		unsigned 		nominal_CIpct;		/* requested confidence interval size (e.g. 95 for 95 percent confidence interval) */
		unsigned 		across;				/* if true, cells added from left to right; otherwise, cells are added top to bottom */
		unsigned 		prnwidth;			/* number of characters (not table columns) in a page */
		unsigned 		curr_col;			/* the column currently being added to the NxsTable */
		unsigned 		ncols;				/* the number of columns */
		unsigned 		curr_row;			/* the row currently being added to the NxsTable */
		unsigned 		nrows;				/* the number of rows */
		bool			*activeCol;			/* an entire column may be hidden by setting corresponding element of this array to false */
		NxsTableCell	**cells;			/* the array of NxsTableCell elements */

	
		void 			Clear();
		void 			AddRows(unsigned = rowincr);
		void 			AddColumns(unsigned = colincr);
		void 			NextCell();
		void 			NextDown();
		void 			NextAcross();
		void 			QSort(unsigned col, unsigned left, unsigned right); // call only via SortColumn (below)
		unsigned		SortColumn(unsigned);
		void 			CheckColumnWidth(unsigned);
		void 			SetUpper(unsigned, unsigned, unsigned);
		void 			SetLower(unsigned, unsigned);
		void 			SetPercent(unsigned, unsigned);
		unsigned		FloatStrLength(double d, unsigned p);
		
	private:
						NxsTable(const NxsTable &);				//never use	- don't define
						NxsTable &operator=(const NxsTable &);	//never use	- don't define
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the desired confidence interval percent `nominal_CIpct' to the supplied value `p'. Assumes `p' is in the range
|	0 to 100 (inclusive).
*/
inline void NxsTable::SetNominalCIPercent(unsigned p)
	{
	NXS_ASSERT(p <= 100);
	nominal_CIpct = p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the desired confidence interval percent `nominal_CIpct' to the supplied value `p'. Assumes `p' is in the range
|	0 to 100 (inclusive).
*/
inline void NxsTable::AddUInt(unsigned i, unsigned w)
	{
	AddInt((int) i, w);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes memory for the `cells' matrix and the `activeCol' array.
*/
inline NxsTable::~NxsTable()
	{
	Clear();
	}


#endif
