#ifndef NCL_NXSTABLECELL_H
#define NCL_NXSTABLECELL_H

#include "phypy/src/ncl/output/nxs_ratio.hpp"
// The macro below specifies the upper limit for the width of any column. Columns wider than this value will be wrapped.
// Used to create a STATIC_DATA character array (cell_workspace) that can be used as a workspace to use in deciding where to 
// wrap the contents of a cell. This should be set to a value large enough that only string cell types could exceed it 
// (i.e.wrapping is only attempted for string cell types).
//
#define ABSOLUTE_MAXIMUM_COLUMN_WIDTH  100

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates a single cell in a NxsTable object. NxsTableCell objects behave differently depending on their type, 
|	which is specified by the value containing in the tag data member. Several of the cell types actually calculate 
|	their values based on the cells above them in the NxsTable object into which they have been inserted. See the
|	NxsTableCellTypeTag enum for the possible NxsTableCell types.
|	
|	Some of the cell types require additional information:
|	
|	The ratio cell type (R) stores the numerator and denominator of a ratio separately (although the quotient is what 
|	is shown in the NxsTable when it is output). The main reason one would use an R cell type over, say, a D type is
|	that the M cell type averages R cell types by summing the numerators of the ratios separately from the denominators,
|	forming the mean by dividing the numerator sum by the denominator sum. This type of averaging behaves better for 
|	ratios in many cases.
|	
|	The hyphens cell type (H) upon output produces a line of hyphens. The number of hyphens output by default equals 
|	the cell width minus one. The line of hyphens is right-justified in the cell. This cell type is useful for 
|	producing lines under column headers that adjust automatically to the width of the cells in that column.
|	
|	The B cell type (stands for Bad value) can be used to output asterisks in cells representing undefined values. The
|	number of asterisks output by default equals the value stored in the data member asterisks. If, however, this value 
|	is larger than the cell width, the number of asterisks will be equal to one less than the cell width.
|
|	See the documentation for NxsTable for an example of how to create and display a NxsTable object.
*/
class NxsTableCell
	{
	friend class NxsTable;

	public:
						NxsTableCell();
						NxsTableCell(const std::string& ss);
						NxsTableCell(const double dd) ;
						NxsTableCell(const int ii);
						NxsTableCell(const long ll);
						NxsTableCell(const NxsRatio rr);

		NxsTableCell	&operator =(const std::string &ss);
		NxsTableCell	&operator =(const int ii);
		NxsTableCell	&operator =(const long ll);
		NxsTableCell	&operator =(const double dd);
		NxsTableCell	&operator =(const NxsRatio rr);
		bool			operator <(const NxsTableCell &);

		void			FactorySettings();
		std::string		&FillWith(unsigned n, char ch = ' ');

		bool			Show(NxsOutputStream &out);
		void			ShowEmpty(NxsOutputStream &out) const;

		enum TableColumnJustification	/* Values determine justification of a column of cells when displayed */
			{
			nojustify,		/* user did not specify justification */
			leftjustify,	/* table column left justified */
			centerjustify,	/* table column center */
			rightjustify	/* table column right justified */
			};

		enum TableCellTypeTag	/* Values determine the type of a particular NxsTableCell object */
			{
			S,		/* string-type cell */
			I,		/* integer-type cell */
			D,		/* double-type cell */
			L,		/* long-type cell */
			R,		/* ratio-type cell */
			B,		/* bad-type cell */
			H,		/* hyphen-type cell */
			E,		/* empty-type cell */
			N,		/* sample-size-type cell */
			SUM,	/* sum-type cell */
			SS,		/* sum-of-squares-type cell */
			M,		/* mean-type cell */
			V,		/* variance-type cell */
			SD,		/* standard-deviation-type cell */
			CV, 	/* coefficient-of-variation-type cell */
			CIU,	/* confidence-interval-upper-bound-type cell */
			CIL,	/* confidence-interval-lower-bound-type cell */
			CIP		/* confidence-interval-percent-type cell */
		};
	
		STATIC_DATA unsigned		def_width;		/* default value for `width' */
		STATIC_DATA unsigned 	def_precision;	/* default value for `precision' */
		STATIC_DATA unsigned 	def_asterisks;	/* default value for `asterisks' */
		unsigned 			width;			/* width of cell; if 0, column will expand to accommodate the widest cell */
		unsigned 			precision;		/* precision for floating point types */
		unsigned 			asterisks;		/* number of asterisks for cells that contain asterisks to indicate undefined values */

	private:
		// Only one of these variables actually used for any given NxsTableCell object (the one used depends on its tag)
		//
		std::string					s;				/* holds contents if cell is type S */
		NxsRatio					r;				/* holds contents if cell is type R */
		double						d;				/* holds contents if cell is type D */
		int							i;				/* holds contents if cell is type I */
		long						l;				/* holds contents if cell is type L */

		TableCellTypeTag			tag;			/* the type of this cell; set to one of the Tag enum values */
		STATIC_DATA char					cell_workspace[ABSOLUTE_MAXIMUM_COLUMN_WIDTH+1];
		unsigned					resume_at;		/* holds position in `s' where output will continue if text needs to be wrapped within column */
		TableColumnJustification	justification;	/* holds one of the values of the TableColumnJustification enum */
	};

#endif
