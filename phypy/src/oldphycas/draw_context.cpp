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

//#include "phycas/force_include.h"
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/oldphycas/draw_context.hpp"
#include "phypy/src/oldphycas/tree.hpp"
#include "phypy/src/oldphycas/tree_node.hpp"
using std::string;
const unsigned char ASCIIDrawContext::codeToAscii[]
	= { ' ', '|', '-', '\\', '|', '|', '/', '+', '-', '/', '-', '\\', '+', '+', '+', '+'}; 
const unsigned char ASCIIDrawContext::codeToPaupMonaco[]
	= { ' ', 222, 220, 222, 223, 221, 223, 225, 220, '9', 220, 227, 'C', 224, 'E', 219, 226, 208 };


void GUIDrawContext::DrawTree(const TreeNode *, bool)
	{
	//canvas->GetSize();
	}
	
ASCIIDrawContext::ASCIIDrawContext(NxsOutputStream &os)
	: compressed(false), 
	usePaupMonaco(false),
	heightOfNextGraphic(0U),
	outStream(os)
	{
	GetPlotWidth();
	}

#if defined (SUPPORT_TREE_OUTPUT)  && defined(NCL_USE_NXS_CONSOLE_OUTPUT)
	 NxsStdOutputStream  &operator<<(NxsStdOutputStream & outStr, const Tree &t)
		{
		ASCIIDrawContext adc(outStr);
		t.Draw(adc);
		return outStr;
		}
#endif


double ASCIIDrawContext::GetHeightScaler(double max_y)
	{ 
	double m =  (compressed ? 1.0 : 2.0);
	heightOfNextGraphic  = (unsigned) max_y*(unsigned) m;
	return m;
	}
	
double ASCIIDrawContext::GetPlotWidth() const
	{
#   if defined (NCL_USE_NXS_CONSOLE_OUTPUT)
		windowWidth = outStream.GetWidth() - 1; //@POL adding the subtraction of 1 to see if that cures the blank line problem
#   else
		windowWidth = 80; //@temp hack will be replaced as we modularize output
#   endif
	return windowWidth;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Translates codes in line prepared by FillLine for output.
*/
void ASCIIDrawContext::TranslateLine(string &lineBuffer)
	{
	const unsigned char *translation = (usePaupMonaco ? codeToPaupMonaco : codeToAscii);
	unsigned lastDarkChar = 0;

	for (unsigned col = 0; col < windowWidth; ++col)
		{
		char ch = lineBuffer[col];
		if (ch > 0)
			lastDarkChar = col;
		if (ch < 16)
			{
			#if 1//
			lineBuffer[col] = (char)translation[ch];
			#else//temp debugging
			if (ch == 0)
				ch = ' ';
			else if (ch < 10)
				ch = '0' + ch;
			else
				ch = 'a' + ch - 10;
			lineBuffer[col] = ch;
			#endif
			}
		}

	lineBuffer[lastDarkChar+1] = '\0';
	}
