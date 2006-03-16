#// dummytest.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <string>
#if defined (_WIN32)
#	include "stdafx.h"
#endif
//#include "ncl/output/nxs_output.hpp"
using std::ostream;
using std::string;
using std::cout;
using std::cin;
using std::endl;


	// enumeration of some basic output styles
enum NxsObjectOutStyle
	{
	kUnspecifiedOutStyle,
	kVerboseOutStyle,
	kConciseOutStyle,
	kGUITabularOutStyle,
	kTabSeparatedTabularOutStyle,
	kCommandOutStyle
	};

	// dummy test class
class A
	{
	public:
		string name;
		A(): name("A"){}
	private:
		A(const A &);//noncopyable
	};
	
	// dummy test class
class B
	{
	public:
		string name;
		B(): name("B"){}
	private:
		B(const B &);//noncopyable
	};
class C: public B
	{};
	
#if defined (_MSC_VER)
		//@POL my version of vc.net does not provide an const string output operator<<
	ostream & operator<<(ostream & out, const string & str);

	inline ostream & operator<<(ostream & out, const std::string & str)
		{
		out << str.c_str();
		return out;
		}
#endif
	
	// 	generic class that is the provides hooks for dispatching output based on 
	//	formatting code, output stream characteristics, and object type
	// 	output should occur in the constructor of specialized versions of 
	//	GenericPrinterClass
template<int FORMAT_HINT, class OBJ, class OUT_STREAM>
class GenericPrinterClass
	{
	public:
		GenericPrinterClass(OUT_STREAM & out, const OBJ & obj);
	};

	
	// 	function that provides front end to GenericPrinterClass, and a
	//	returns the stream again
template<int FORMAT_HINT, class OBJ, class OUT_STREAM>
OUT_STREAM & Emit(OUT_STREAM & out, const OBJ & obj);

template<int FORMAT_HINT, class OBJ, class OUT_STREAM>
inline OUT_STREAM & Emit(OUT_STREAM & out, const OBJ & obj)
	{
	typedef GenericPrinterClass<FORMAT_HINT, OBJ, OUT_STREAM> GPC;
	GPC gpc(out, obj);
	return out;
	}

	// 	2-template parameter version of the emit function fills in 
	//	kUnspecifiedOutStyle for the missing FORMAT_HINT 
template<class OBJ, class OUT_STREAM>
OUT_STREAM & Emit(OUT_STREAM & out, const OBJ & obj);

template<class OBJ, class OUT_STREAM>
OUT_STREAM & Emit(OUT_STREAM & out, const OBJ & obj)
	{
	typedef GenericPrinterClass<kUnspecifiedOutStyle, OBJ, OUT_STREAM> GPC;
	GPC gpc(out, obj);
	return out;
	}


	// 	example of a minimal specialization for any object of class A
	//	this overload will suck up any output calls that do not get overloaded
	//	with a more specific form
template<int FORMAT_HINT, class OUT_STREAM>
class GenericPrinterClass<FORMAT_HINT, A, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & out, const A & obj)
			{
			out << "name is " << obj.name << '\n';
			}
	};

	// 	specialization for verbose output of A to any output stream
template<class OUT_STREAM>
class GenericPrinterClass<kVerboseOutStyle, A, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & out, const A & obj)
			{
			out << "I'm quite sure that the name is " << obj.name.c_str() << '\n';
			}
	};

	// 	specialization that relays kUnspecifiedOutStyle to the kVerboseOutStyle spec.
template<class OUT_STREAM>
class GenericPrinterClass<kUnspecifiedOutStyle, A, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & out, const A & obj)
			{
			GenericPrinterClass<kVerboseOutStyle, A, OUT_STREAM>(out, obj);
			}
	};

	// 	specialization for kUnspecifiedOutStyle of B
template<class OUT_STREAM>
class GenericPrinterClass<kUnspecifiedOutStyle, B, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & out, const B & obj)
			{
			out << "the name (" << obj.name.c_str() << ") should be B.\n";
			}
	};

/// broad templated support of the << operator seems prone to interfere with other 
///	uses of this operator:
///	template<class STREAM, class OBJ>
///	STREAM & operator<<(STREAM &, const OBJ &);
///	would override operators that use inheritance (I'm not sure about this, but I think 
/// that it would get hairy).  Additionally, our use of << to concatenate to strings
///	might then get funnelled through this operator.
///	Instead of generically templating << , we can provide syntactic sugar for outputting:
///		cout << objectFormatter<kVerboseOutStyle>(a);
/// 

	///	Utility class used to implement the objectFormatter output through the << operator
	///	type that encapsulates the object to be printed, its type and the format code.
template<int OUT_STYLE, typename T>
class _ObjectFormatterImpl
	{
	public:
		explicit _ObjectFormatterImpl(const T & outObj): obj(outObj) {}
		const T & obj;
	};

	/// declaration and definitoin of the templated output operator for _ObjectFormatterImpl objects
template<int OUT_STYLE, typename T, class OUT_STREAM>
OUT_STREAM & operator<<(OUT_STREAM & outStream,  const _ObjectFormatterImpl<OUT_STYLE, T> &ofi);
template<int OUT_STYLE, typename T, class OUT_STREAM>
inline OUT_STREAM & operator<<(OUT_STREAM & outStream,  const _ObjectFormatterImpl<OUT_STYLE, T> &ofi)
	{
	return Emit<OUT_STYLE, T, OUT_STREAM>(outStream, ofi.obj);
	}
	
	///	Templated function that bind stream type, object type and object reference
	/// Function can be used instead of the class definitions because compilers will 
	/// infer the full template argument signature from:
	///		objectFormatter<kVerboseOutputFormat>(object);
	/// but would require
	///		_ObjectFormatterImpl<kVerboseOutputFormat, OBJECT_TYPE>(object);
template<int OUT_STYLE, typename T>
_ObjectFormatterImpl<OUT_STYLE, T> objectFormatter(const T & obj);

template<int OUT_STYLE, typename T>
inline _ObjectFormatterImpl<OUT_STYLE, T> objectFormatter(const T & obj)
	{
	return _ObjectFormatterImpl<OUT_STYLE, T>(obj);
	}

/*	/// deprecated, non-templated version ofth 
template<int OUT_STYLE, typename T>
NxsWritableStream & operator<<(NxsWritableStream & outStream, const _ObjectFormatterImpl<OUT_STYLE, T> & of);

template<int OUT_STYLE, typename T>
std::ostream & operator<<(std::ostream & outStream, const _ObjectFormatterImpl<OUT_STYLE, T> & of);

template<int OUT_STYLE, typename T>
NxsWritableStream & operator<<(NxsWritableStream & outStream, const _ObjectFormatterImpl<OUT_STYLE, T> & of)
	{
	return Emit<OUT_STYLE, T, NxsWritableStream>(outStream, of.obj);
	}
	
template<int OUT_STYLE, typename T>
std::ostream & operator<<(std::ostream & outStream, const _ObjectFormatterImpl<OUT_STYLE, T> & of)
	{
	return Emit<OUT_STYLE, T, std::ostream>(outStream, of.obj);
	}
*/

const char * expectedOutput = "Verbose output of a:\n"
"	I\'m quite sure that the name is A\n"
" << operator style verbose output of a:\n"
"	I\'m quite sure that the name is A\n"
"Concise output of a via low priority templated form:\n"
"	name is A\n"
"kUnspecifiedOutStyle output of a:\n"
"	I\'m quite sure that the name is A\n"
"kUnspecifiedOutStyle output of b:\n"
"	the name (B) should be B.\n"
"kUnspecifiedOutStyle output of c via inheritance:\n"
"	the name (B) should be B.\n";

#if defined (_WIN32)
int _tmain(int argc, _TCHAR* argv[])
	{
#else
int main(int, char *[])
	{
#endif
	cout << "Expected output:";
	cout << "\n========================================\n";
	cout << expectedOutput;
	cout << "\n========================================\n";
	A a;
	B b;

	cout << "Verbose output of a:\n\t";
	Emit<kVerboseOutStyle>(cout, a);

	cout << " << operator style verbose output of a:\n\t";
	cout << objectFormatter<kVerboseOutStyle>(a);
		// bare output via << is not supported (unless you define an operator for A)
		// this - generates compile-time error
		//	cout << a;
	cout << "Concise output of a via low priority templated form:\n\t";
	Emit<kConciseOutStyle>(cout, a);

	cout << "kUnspecifiedOutStyle output of a:\n\t";
	Emit(cout, a);


	cout << "kUnspecifiedOutStyle output of b:\n\t";
	Emit(cout, b);
	//	Use of object without low priority definition or specialization.
	//	this results in a link error.	
	//	Emit<kVerboseOutStyle>(cout, b);
	C c;
	cout << "kUnspecifiedOutStyle output of c via inheritance:\n\t";
	Emit<B>(cout, c);
	
	
	//	Lack of specializaiton to handle kVerboseOutStyle applies to C as well
		//	Link error
		//Emit<kVerboseOutStyle, B>(cout, c);
	
	//	Specification of the incorrect base class is compile-time error
		//Emit<kVerboseOutStyle, A>(cout, c);
	
	
#	if defined (_WIN32)
		// prompt for a key so that the console does not disappear
		char c;
		cout << "Enter a key to quit" << endl;
		cin >> c;
#	endif
	return 0;
	}

