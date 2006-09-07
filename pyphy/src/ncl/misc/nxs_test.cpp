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



