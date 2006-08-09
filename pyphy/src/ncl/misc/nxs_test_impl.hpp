#if !defined NCL_NXS_TEST_IMPLEMENTATION_H
#define NCL_NXS_TEST_IMPLEMENTATION_H

#include "ncl/misc/nxs_test.hpp"
#include <boost/function.hpp>
#include "ncl/misc/nxs_copy.hpp"
#include "ncl/misc/utilities.hpp"
#include "ncl/misc/string_extensions.hpp"
/*----------------------------------------------------------------------------------------------------------------------
|	To avoid code-bloat the NxsTest class is abstract and not templated.  The implemenation of the tests is templated.
|
|	There are 3 general types of tests: NxsCompareTest (for <, <=, ==, >=, >, or != tests), NxsBoolTest (for calling a
|	function that returns bool or checking a bool variable), and NxsBinaryTest (for tests that are logical combinations 
|	of other tests, e.g. OR, XOR, IfFirstThenSecond...)
|
|	Error messages when test fail are supplied by a MessageSource class. Only 2 are currently implemented StoredMsgSource
|	and DummyMsgSource (the latter can be used for internal test that will never be required to explain themselves to the user.
|
|	NxsCompareTest has 3 template parameters: type being compared, LeftOperandPolicy, RightOperandPolicy.
|	The operator is determined by picking a value from the NxsTest::ComparisonOperator enum and specifying it in the constructor.
|	NxsConstValueWrapper, NxsVariableWrapper, and NxsFuncWrapper are 3 very simple inlined classes that can be used as Left or Right
|	operands that return a value of the appropriate type that is a constant supplied at construction (NxsConstValueWrapper),
|	the current value held in an address specified at construction (NxsVariableWrapper), or the result of function call 
|	NxsFuncWrapper.
|	To avoid code bloat, constants should be the right operand.  Typedefs at the end of this file should be used when
|	declaring or defining new test.
*/



/// interface for providing messages to the user if a comparison test fails
class NxsComparsionTestMessageSource
	{
	public:
		virtual ~NxsComparsionTestMessageSource(){}
		virtual NxsComparsionTestMessageSource *Clone() const = 0;
		virtual std::string GenerateMessage(NxsTest::ComparisonOperator o, const std::string &leftVal, const std::string &rightVal, NxsTest::msgContext context) = 0;
	};

/// interface for providing messages to the user if a comparison test fails
class NxsBinaryTestMessageSource
	{
	public:
		virtual ~NxsBinaryTestMessageSource(){}
		virtual NxsBinaryTestMessageSource *CloneBinaryTestMsgSource() const = 0;
		virtual std::string GenerateBinaryTestMessage(NxsTest::BinaryOperator , NxsTest::msgContext context) = 0;
	};


/// wrapper so that compile-time constants can be used as operands in tests
template <typename T> class NxsConstValueWrapper
	{
		const T val;
	public:
		NxsConstValueWrapper(const T & v)
		  :val(v)
		  	{
		  	}
		T operator()() const 
			{
			return val;
			}
	};
	
/// wrapper so that variables can be used as operands in tests
template <typename T> class NxsVariableWrapper
	{
		const T * valueAddress;
	public:
		NxsVariableWrapper(const T *v)
		  :valueAddress(v)
		  	{
		  	assert(v != NULL);
		  	}
		T operator()() const 
			{
			return *valueAddress;
			}
	};
	
/// wrapper for function callback to allow them to act as operands in tests
/// \todeprecate If the Test interface is stable, we can use boost::function0<T> instead of NxsFuncWrapper
template <typename T> class NxsFuncWrapper
	{
		boost::function0<T> functionPtr;
	public:
		NxsFuncWrapper(boost::function0<T> v)
		  :functionPtr(v)
		  	{
		  	assert(v);
		  	}
		T operator()() const 
			{
			return (functionPtr)();
			}
	};

	
/// implementation of NxsTest interface for tests that involve numerical comparisons between operands
template<typename ComparedType, 
		template<typename> class LeftOperandTemplate,
		template<typename> class RightOperandTemplate > 
	class NxsCompareTest : public NxsTest
	{
	public :
	
		//typedef boost::function0< ComparedType >  ValueSource;
			
		typedef LeftOperandTemplate<ComparedType> LeftOperand;
		typedef RightOperandTemplate<ComparedType> RightOperand;
		
		NxsCompareTest(const LeftOperand &leftOp, NxsTest::ComparisonOperator testOperator, const RightOperand &rightOp, const NxsComparsionTestMessageSource &m)
			:leftOperand(leftOp),
			op(testOperator),
			rightOperand(rightOp),
			msgSource(m.Clone())
			{
			}
		~NxsCompareTest()
			{
			delete msgSource;
			}

		bool 		operator()()
			{
			if(op == NxsTest::ComparisonOperator(kLessThan))
				return leftOperand() < rightOperand();
			if(op == NxsTest::ComparisonOperator(kLessOrEq))
				return leftOperand() <= rightOperand();
			if(op == NxsTest::ComparisonOperator(kEquals))
				return leftOperand() == rightOperand();
			if(op == NxsTest::ComparisonOperator(kGreaterOrEq))
				return leftOperand() >= rightOperand();
			if(op == NxsTest::ComparisonOperator(kGreaterThan))
				return leftOperand() > rightOperand();
			assert(op == NxsTest::ComparisonOperator(kNotEqual));
			return leftOperand() != rightOperand();
			}
		std::string 	GetMessage(NxsTest::msgContext context = explain_failure)
			{
			std::string lStr, rStr;
			lStr << leftOperand();
			rStr << rightOperand();
			return msgSource->GenerateMessage(op, lStr, rStr, context);
			}



	protected :
		LeftOperand			leftOperand;
		NxsTest::ComparisonOperator op;
		RightOperand		rightOperand;
		NxsComparsionTestMessageSource *msgSource;
	};	

#if 0
	template<typename ComparedType, 
			template<typename> class LeftOperandTemplate,
			template<typename> class RightOperandTemplate > 
	inline NxsCompareTest<ComparedType, LeftOperandTemplate, RightOperandTemplate>::NxsCompareTest(
	  const LeftOperand &leftOp,
	  NxsTest::ComparisonOperator testOperator,
	  const RightOperand &rightOp,
	  const NxsComparsionTestMessageSource &m)
		:leftOperand(leftOp),
		op(testOperator),
		rightOperand(rightOp),
		msgSource(m.Clone())
		{
		}
		
	template<typename ComparedType, 
			template<typename> class LeftOperandTemplate,
			template<typename> class RightOperandTemplate > 
	bool NxsCompareTest<ComparedType, LeftOperandTemplate, RightOperandTemplate>::operator()()
		{
		switch(op)
			{
			case(NxsTest::ComparisonOperator(kLessThan)) : return leftOperand() < rightOperand();
			case(NxsTest::ComparisonOperator(kLessOrEq)) : return leftOperand() <= rightOperand();
			case(NxsTest::ComparisonOperator(kEquals)) 	: return leftOperand() == rightOperand();
			case(NxsTest::ComparisonOperator(kGreaterOrEq)) : return leftOperand() >= rightOperand();
			case(NxsTest::ComparisonOperator(kGreaterThan)) : return leftOperand() > rightOperand();
			}
		assert(op == NxsTest::ComparisonOperator(kNotEqual));
		return leftOperand() != rightOperand();
		}

	template<typename ComparedType, 
			template<typename> class LeftOperandTemplate,
			template<typename> class RightOperandTemplate > 
	std::string NxsCompareTest<ComparedType, LeftOperandTemplate, RightOperandTemplate>::GetMessage(
	  NxsTest::msgContext context)
		{
		std::string lStr, rStr;
		lStr << leftOperand();
		rStr << rightOperand();
		return msgSource->GenerateMessage(op, lStr, rStr, context);
		}

	template<typename ComparedType, 
			template<typename> class LeftOperandTemplate,
			template<typename> class RightOperandTemplate > 
	NxsCompareTest<ComparedType, LeftOperandTemplate, RightOperandTemplate>::~NxsCompareTest()
		{
		delete msgSource;
		}
#endif		

/// implementation of NxsTest interface for tests that involve logical operations between operands
template<NxsTest::BinaryOperator binOperator> 
class NxsBinaryTest : public NxsTest
	{
	public :
		NxsBinaryTest(NxsTestShPtr leftTest, NxsTestShPtr rightTest, const NxsBinaryTestMessageSource & m)
			:leftOperand(leftTest),
			rightOperand(rightTest),
			msgSource(m.CloneBinaryTestMsgSource())
			{
			}
		
		bool operator()()
			{
			return RunTest(Int2Type<binOperator>());
			}

		
		bool RunTest(Int2Type<NxsTest::BinaryOperator(kXor)>)
			{
			if ((*leftOperand)())
				return !(*rightOperand)(); 
			else
				return (*rightOperand)() ;
			}
			
		bool RunTest(Int2Type<NxsTest::BinaryOperator(kBoth)>)
			{
			if ((*leftOperand)())
				return (*rightOperand)(); 
			return false;
			}
			
		bool RunTest(Int2Type<NxsTest::BinaryOperator(kNeither)>)
			{
			if (!(*leftOperand)())
				return !(*rightOperand)(); 
			return false;
			}
			
		bool RunTest(Int2Type<NxsTest::BinaryOperator(kIfFirstSecond)>)
			{
			const bool firstResult = (*leftOperand)();
			return  ( !firstResult ||  (*rightOperand)());
			}
			
		bool RunTest(Int2Type<NxsTest::BinaryOperator(kNotOne)>)
			{
			const bool firstResult = (*leftOperand)();
			const bool secondResult = (*rightOperand)();
			return ((firstResult && secondResult) || (!firstResult && !secondResult));
			}
			
		bool RunTest(Int2Type<NxsTest::BinaryOperator(kNotBoth)>)
			{
			const bool firstResult = (*leftOperand)();
			return ( !firstResult ||  ! ((*rightOperand)()));
			}
			
		bool RunTest(Int2Type<NxsTest::BinaryOperator(kEither)>)
			{
			return ((*leftOperand)() || (*rightOperand)());
			}
		
		std::string GetMessage(NxsTest::msgContext context = explain_failure)
			{
			return msgSource->GenerateBinaryTestMessage(binOperator, context);
			};
		
		~NxsBinaryTest()
			{
			delete msgSource;
			DeletePtr(leftOperand);
			DeletePtr(rightOperand);
			}
	protected :
		NxsTestShPtr leftOperand;
		NxsTestShPtr rightOperand;
		NxsBinaryTestMessageSource * msgSource;
	};

class NxsNamedOperandMsgSource : public NxsComparsionTestMessageSource
	{
	protected:
		const std::string leftOperandName;
		const std::string rightOperandName;
		std::string GetOperandString(const std::string & opName, const std::string &opVal)
			{
			if (opName.empty())
				return opVal;
			std::string s = opName ;
			s << " (currently " << opVal << ')';
			return s;
			}		
	public:
		NxsNamedOperandMsgSource(const std::string &leftOpName, const std::string &rightOpName)
			:leftOperandName(leftOpName),
			rightOperandName(rightOpName)
			{
			}
			
		std::string GenerateMessage(NxsTest::ComparisonOperator o, const std::string &leftVal, const std::string &rightVal, NxsTest::msgContext context);
	};
	
class StoredMsgSource : public NxsComparsionTestMessageSource, public NxsBinaryTestMessageSource
	{
	protected:
		const std::string failMsg;
		const std::string successMsg;
	
	public:
		StoredMsgSource(const std::string &failureMessage, const std::string &successMessage)
			:failMsg(failureMessage),
			successMsg(successMessage)
			{
			}
			
		std::string GenerateMessage(NxsTest::ComparisonOperator , const std::string &, const std::string &, NxsTest::msgContext context)
			{
			if (context == NxsTest::explain_failure)
				return failMsg;
			return  successMsg;
			}
		std::string GenerateBinaryTestMessage(NxsTest::BinaryOperator , NxsTest::msgContext context)
			{
			if (context == NxsTest::explain_failure)
				return failMsg;
			return  successMsg;
			}
		
		NxsComparsionTestMessageSource * Clone() const 
			{
			return new StoredMsgSource(failMsg, successMsg);
			}
		NxsBinaryTestMessageSource * CloneBinaryTestMsgSource() const 
			{
			return new StoredMsgSource(failMsg, successMsg);
			}
	
	};

class DummyMsgSource : public NxsComparsionTestMessageSource , public NxsBinaryTestMessageSource
	{
	public:
		std::string GenerateMessage(NxsTest::ComparisonOperator , const std::string &, const std::string &, NxsTest::msgContext)
			{
			return  std::string();
			}
		std::string GenerateBinaryTestMessage(NxsTest::BinaryOperator , NxsTest::msgContext)
			{
			return  std::string();
			}
		NxsComparsionTestMessageSource *Clone() const {return new DummyMsgSource();}
		NxsBinaryTestMessageSource *CloneBinaryTestMsgSource() const {return new DummyMsgSource();}
		
	};

template<class BoolSource> 
	class NxsBoolTest : public NxsTest
	{
	public :
	
			
		NxsBoolTest(BoolSource boolSource,const StoredMsgSource &m)
			:valueSource(boolSource),
			msgSource(m.CloneBinaryTestMsgSource())
			{
			}
		NxsBoolTest(BoolSource boolSource,const DummyMsgSource &m)
			:valueSource(boolSource),
			msgSource(m.CloneBinaryTestMsgSource())
			{
			}
		
		bool operator()()
			{
			return valueSource();
			}
		std::string GetMessage(NxsTest::msgContext context = explain_failure)
			{
			return msgSource->GenerateBinaryTestMessage(NxsTest::BinaryOperator(kFirstOnly), context);
			};
	
		~NxsBoolTest()
			{
			delete msgSource;
			}
			
	protected :
		BoolSource  	valueSource;
		NxsBinaryTestMessageSource *msgSource;
	};	


//	Lots of typedefs of the named form Type_LeftOp_RightOpt_Test 
//	
typedef NxsFuncWrapper<UInt>	NxsUIntCallbackWrapper;
typedef NxsFuncWrapper<int>		NxsIntCallbackWrapper;
typedef NxsFuncWrapper<double>  NxsDblCallbackWrapper;
typedef NxsFuncWrapper<bool>	NxsBoolCallbackWrapper;

typedef NxsVariableWrapper<UInt>	NxsUIntVarWrapper;
typedef NxsVariableWrapper<int>		NxsIntVarWrapper;
typedef NxsVariableWrapper<double>  NxsDblVarWrapper;
typedef NxsVariableWrapper<bool>	NxsBoolVarWrapper;

typedef NxsConstValueWrapper<UInt>		NxsUIntValueWrapper;
typedef NxsConstValueWrapper<int>		NxsIntValueWrapper;
typedef NxsConstValueWrapper<double>	NxsDblValueWrapper;
typedef NxsConstValueWrapper<bool>		NxsBoolValueWrapper;

typedef NxsCompareTest<unsigned, NxsFuncWrapper, NxsConstValueWrapper > 		UIntFuncToConstTest;
typedef NxsCompareTest<unsigned, NxsFuncWrapper, NxsFuncWrapper> 				UIntFuncToFuncTest;
typedef NxsCompareTest<unsigned, NxsFuncWrapper, NxsVariableWrapper>			UIntFuncToVarTest;
typedef NxsCompareTest<unsigned, NxsVariableWrapper, NxsConstValueWrapper> 		UIntVarToConstTest;
typedef NxsCompareTest<unsigned, NxsVariableWrapper, NxsFuncWrapper>			UIntVarToFuncTest;
typedef NxsCompareTest<unsigned, NxsVariableWrapper, NxsVariableWrapper> 		UIntVarToVarTest;

typedef NxsCompareTest<int, NxsFuncWrapper, NxsConstValueWrapper> 		IntFuncToConstTest;
typedef NxsCompareTest<int, NxsFuncWrapper, NxsFuncWrapper> 			IntFuncToFuncTest;
typedef NxsCompareTest<int, NxsFuncWrapper, NxsVariableWrapper>			IntFuncToVarTest;
typedef NxsCompareTest<int, NxsVariableWrapper, NxsConstValueWrapper>	IntVarToConstTest;
typedef NxsCompareTest<int, NxsVariableWrapper, NxsFuncWrapper>			IntVarToFuncTest;
typedef NxsCompareTest<int, NxsVariableWrapper, NxsVariableWrapper> 	IntVarToVarTest;

typedef NxsCompareTest<double, NxsFuncWrapper, NxsConstValueWrapper> 		DblFuncToConstTest;
typedef NxsCompareTest<double, NxsFuncWrapper, NxsFuncWrapper> 				DblFuncToFuncTest;
typedef NxsCompareTest<double, NxsFuncWrapper, NxsVariableWrapper>			DblFuncToVarTest;
typedef NxsCompareTest<double, NxsVariableWrapper, NxsConstValueWrapper> 	DblVarToConstTest;
typedef NxsCompareTest<double, NxsVariableWrapper, NxsFuncWrapper>			DblVarToFuncTest;
typedef NxsCompareTest<double, NxsVariableWrapper, NxsVariableWrapper> 		DblVarToVarTest;

typedef NxsBinaryTest<NxsTest::kXor>			Nxs_Xor_Test;
typedef NxsBinaryTest<NxsTest::kIfFirstSecond> 	Nxs_If1st2nd_Test;
typedef NxsBinaryTest<NxsTest::kNotOne>			Nxs_NotOne_Test;
typedef NxsBinaryTest<NxsTest::kNotBoth>		Nxs_NotBoth_Test;
typedef NxsBinaryTest<NxsTest::kEither>			Nxs_Or_Test;		//standard || test
typedef NxsBinaryTest<NxsTest::kBoth> 			Nxs_And_Test;		//standard && test
typedef NxsBinaryTest<NxsTest::kNeither>		Nxs_Neither_Test;		

typedef NxsBoolTest< NxsFuncWrapper<bool> > 					BoolFuncTest;
typedef NxsBoolTest< NxsVariableWrapper<bool> > 				BoolVarTest;
#endif

