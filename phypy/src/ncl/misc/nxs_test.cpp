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
#include "phypy/src/ncl/misc/nxs_test_impl.hpp"
using std::string;
/*--------------------------------------------------------------------------------------------------------------------------
|	Walks through the vector of tests, checking each one.  if returns false, then errM is set to the first
|	message from test that failed
*/
bool PerformAllTests(
  std::vector<NxsTestShPtr>	& reqVec, std::string * errM)
  	{
  	for (std::vector<NxsTestShPtr>::iterator rIt = reqVec.begin(); rIt != reqVec.end(); ++rIt)
  		{
  		if (! (**rIt)())
  			{
  			if (errM != NULL)
				*errM = (*rIt)->GetMessage();
  			return false;
  			}
  		}
  	return true;
  	}

string NxsNamedOperandMsgSource::GenerateMessage(NxsTest::ComparisonOperator o, const string &leftVal, const string &rightVal, NxsTest::msgContext context)
	{
	string s;
	NxsTest::ComparisonOperator oppositeOp;
	switch(o)
		{
		case (NxsTest::kLessThan) :
			s << GetOperandString(leftOperandName, leftVal) << " is ";
			if (context == NxsTest::explain_failure)
				s << "not ";
			s << "less than " << GetOperandString(rightOperandName, rightVal);
			return s;
		case (NxsTest::kLessOrEq) :
			s << GetOperandString(leftOperandName, leftVal) << " is ";
			if (context == NxsTest::explain_success)
				s << "not ";
			s << "greater than " << GetOperandString(rightOperandName, rightVal);
			return s;
		case (NxsTest::kEquals) :
			s << GetOperandString(leftOperandName, leftVal) << " is ";
			if (context == NxsTest::explain_failure)
				s << "not ";
			s << "equal to " << GetOperandString(rightOperandName, rightVal);
			return s;
		case (NxsTest::kGreaterOrEq) : oppositeOp = NxsTest::kLessThan;
								break;
		case (NxsTest::kGreaterThan) : oppositeOp =NxsTest:: kLessOrEq;
								break;
		default	: oppositeOp = NxsTest::kEquals;
		}
	context = (context == NxsTest::explain_success ? NxsTest::explain_failure : NxsTest::explain_success);
	return GenerateMessage(oppositeOp, leftVal, rightVal, context);
	}



