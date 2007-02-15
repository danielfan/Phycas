/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
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

#ifndef NCL_NXS_TEST_H
#define NCL_NXS_TEST_H

#include "phycas/src/ncl/nxs_defs.hpp"
class NxsTest;
typedef boost::shared_ptr<NxsTest>			NxsTestShPtr;//NOTE:this typedef is also in nxscommandcommon.h change it in both places

/*----------------------------------------------------------------------------------------------------------------------
|	NxsTest is an interface for a callback functions that perform run-time checks.
|	Derived classes must overload the function operator (taking void and returning a bool) and the GetMessgae function
|	that supplies a phrase to the user explaining why the test evaluated failed or succeeded
|	Templated classes override Nxstest to make it easy tailor the test, but the interface is used to cut down on code
|	bloat for classes that use the tests (they don't need to be templated).
*/
class NxsTest
	{
	public :
		// NxsCompareTestOperator and NxsBinaryTestOperator must be kept sequential (some functionality depends on them behaving
		//	as if they were one enum, but they are split for better type-checking).
		enum ComparisonOperator
			{
			kLessThan = 0,	
			kLessOrEq = 1,
			kEquals = 2,
			kGreaterOrEq = 3, 
			kGreaterThan = 4, 
			kNotEqual = 5
			};

		enum BinaryOperator
			{
			kXor = 6,	
			kIfFirstSecond = 7,
			kNotOne = 8,
			kNotBoth = 9, 
			kEither = 10,
			kBoth = 11,
			kNeither = 12,
			kFirstOnly = 13
			};

		enum msgContext 		/* Argument for GetMessage  */
			{
			explain_failure, 	/* GetMessage should return a string to complete the phrase "test failed because ..." */
			explain_success 	/* GetMessage should return a string to complete the phrase "test failed because ..." */
			};
			
		virtual bool operator()() = 0;
		virtual std::string 		GetMessage(msgContext context = explain_failure) = 0;	
		
		virtual ~NxsTest()	{} 
		
	protected:
			std::string		failureMessage;
			std::string		successMessage;
		
	};
	
bool PerformAllTests(std::vector<NxsTestShPtr> & reqVec, std::string * errM);	



#endif
