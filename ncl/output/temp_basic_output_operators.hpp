#if !defined (TEMP_BASIC_OUTPUT_OPERATORS_HPP)
#define TEMP_BASIC_OUTPUT_OPERATORS_HPP
#include "ncl/misc/generic_type_mapping.hpp"
#include "ncl/misc/compile_assert.hpp"
#include "ncl/output/nxs_output_stream_wrapper.hpp"

	// declare primitive output operators
NxsWritableStream & operator<<(NxsWritableStream & outS, const std::string & s);
NxsWritableStream & operator<<(NxsWritableStream & outS, const char  s);
NxsWritableStream & operator<<(NxsWritableStream & outS, const double s);
NxsWritableStream & operator<<(NxsWritableStream & outS, const int s);
NxsWritableStream & operator<<(NxsWritableStream & outS, const unsigned s);
typedef NxsWritableStream & (* WritableStreamModFunc)(NxsWritableStream &);
NxsWritableStream & operator<<(NxsWritableStream &, WritableStreamModFunc);

/**
 *  The following enumerations are formatting hints for the output of objects using the Emit command
 *  This is for formatting decisions known at compile-time.
 *  NxsObjectOutStyle apply to a wide variety of objects
 *  different directories of the ncl tree can have specialized formatting hints for their types using 
 *  enums such as NxsMiscDirOutStyle for the ncl/misc/ directory
 *  Uniqueness of these numbers is guaranteed by using directory-specific starting indices listed 
 *  in NxsOutStyleOffsets (we can add more as needed)
*/
enum NxsOutStyleOffsets
	{
	kNCLMiscOutStyleOffset	= 0x0100,
	kNCLTaxaOutStyleOffset	= 0x0200,
	kNCLTreesOutStyleOffset	= 0x0300
	};
	// enumeration of some basic output styles
enum NxsObjectOutStyle
	{
	kUnspecifiedOutStyle = 0x00,	/// (default) no formatting was specified (there is only one way to print this object, or the caller has no preference)
	kVerboseOutStyle = 0x01,		/// used to produce output shown to the user that explains itself as fully as possible and is formatted to be easy to be read by humans
	kConciseOutStyle = 0x02,		/// print a brief summary
	kGUITabularOutStyle = 0x03,		/// print an ob
	kGUITabularLabelsOutStyle = 0x04,
	kTabSeparatedTabularOutStyle = 0x05,
	kTabSeparatedTabularLabelsOutStyle = 0x06,
	kCommandOutSyle = 0x07,
	kCommandStateOutSyle = 0x08,
	};
	
enum NxsMiscDirOutStyle
	{
	kVerboseSetOutStyle				= kVerboseOutStyle,					/// prints a set using user-defined labels (instead of numbers): Homo_sapiens Pongo Drosophila
	kConciseSetOutStyle				= kConciseOutStyle,					/// prints the set definition with numbers and ranges:  1 4 7-22
	k01BitSetOutStyle				= kNCLMiscOutStyleOffset + 0x01,	/// prints set membership using a bit representation where 1 indicates a member: 0100100001
	kDotAsteriskBitSetOutStyle		= kNCLMiscOutStyleOffset + 0x02,	/// similar to k01BitSetOutStyle but using . and * (like PAUP's bipartition table): .**.****..
	kHexBitSetOutStyle				= kNCLMiscOutStyleOffset + 0x04,	/// similar to k01BitSetOutStyle, but print a character string 4 bits at a time.  The string can be interpretted as a hex string: 1A7B5
	kRawBitSetOutStyle				= kNCLMiscOutStyleOffset + 0x05,		/// similar to k01BitSetOutStyle, produces a string of characters which should treated as integers (may contain unprintable characters when treated as ASCII)
	kCommandStateKnownSetOutStyle   = kNCLMiscOutStyleOffset + 0x06
	};

enum NxsTaxaDirOutStyle
	{
	kNxsTreeTranslateCmdOutStyle	= kNCLTaxaOutStyleOffset + 0x01		/// prints a TRANSLATE command that can go in the trees block
	};
	
enum NxsTreesDirOutStyle
	{
	kNewickTreeOutStyle				= kNCLTreesOutStyleOffset + 0x01, /// print just the newick description
	kNexusTreeCmdOutStyle			= kNCLTreesOutStyleOffset + 0x02, /// print "tree {NAME} = {NEWICK};"
	kGraphicalTreeOutStyle			= kNCLTreesOutStyleOffset + 0x03, /// uses ASCII tree to draw the tree
	kSplitSummaryTable				= kNCLTreesOutStyleOffset + 0x04,
	kVerboseModelParameterOutStyle	= kNCLTreesOutStyleOffset + 0x50
	};
	
	// 	generic class that is the provides hooks for dispatching output based on 
	//	formatting code, output stream characteristics, and object type
	// 	output should occur in the constructor of specialized versions of 
	//	GenericPrinterClass
template<int FORMAT_HINT, class OBJ, class OUT_STREAM>
class GenericPrinterClass
	{
	public:
		GenericPrinterClass(OUT_STREAM & out, const OBJ & obj)
			{
			COMPILE_TIME_ASSERT(false, NoOverloadedOutputClassForEmitCall);
			}
	};

	
	// 	function that provides front end to GenericPrinterClass, and a
	//	returns the stream again
template<int FORMAT_HINT, class OBJ, class OUT_STREAM>
OUT_STREAM & Emit(OUT_STREAM & out, const OBJ & obj);
template<int FORMAT_HINT, class OBJ>
inline NxsOutputStreamWrapper & EmitGeneric(NxsOutputStreamWrapper & out, const OBJ & obj);
template<class OBJ>
inline NxsOutputStreamWrapper & EmitGenericPlotData(NxsOutputStreamWrapper & out, const OBJ & obj);
template<class OBJ>
inline NxsOutputStreamWrapper & EmitGenericPlotLabels(NxsOutputStreamWrapper & out, const OBJ & obj);

template<int FORMAT_HINT, class OBJ, class OUT_STREAM>
inline OUT_STREAM & Emit(OUT_STREAM & out, const OBJ & obj)
	{
	typedef GenericPrinterClass<FORMAT_HINT, OBJ, OUT_STREAM> GPC;
	GPC gpc(out, obj); //@POL 15-Nov-2005 gcc complains about gpc being unused - I was afraid to touch this however
	return out;
	}

template<int FORMAT_HINT, class OBJ>
inline NxsOutputStreamWrapper & EmitGeneric(NxsOutputStreamWrapper & out, const OBJ & obj)
	{
	if (out.GetFilePtr())
		Emit<FORMAT_HINT, OBJ, std::ofstream>(*out.GetFilePtr(), obj);
	const unsigned int np = out.GetNumNCLStreamPtr();
	for (unsigned int i = 0; i < np; ++i)
		Emit<FORMAT_HINT, OBJ, NxsWritableStream>(*out.GetNCLStreamPtr(i), obj);
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

#if !defined(NCL_USE_NXS_CONSOLE_OUTPUT) && ! defined(NCL_USE_STD_OUTPUT)
#   include "ncl/output/nxs_sax_output_wrapper.hpp" 
#endif
template<class OBJ>
inline NxsOutputStreamWrapper & EmitGenericPlotLabels(NxsOutputStreamWrapper & out, const OBJ & obj)
	{
	EmitGeneric<kTabSeparatedTabularLabelsOutStyle, OBJ>(out, obj);
#   if !defined(NCL_USE_NXS_CONSOLE_OUTPUT) && ! defined(NCL_USE_STD_OUTPUT)
		// in GUI mode we might have plot stream to deal with 
		if (out.GetPlotStreamPtr())
			{
			NxsSaxOutputWrapper nsow(*out.GetPlotStreamPtr(), "label"); //@temp hard coding the xml here isn't optimal
			Emit<kTabSeparatedTabularLabelsOutStyle, OBJ, NxsWritableStream>(*out.GetPlotStreamPtr(), obj);
			}
#   endif
	return out;
	}

template<class OBJ>
inline NxsOutputStreamWrapper & EmitGenericPlotData(NxsOutputStreamWrapper & out, const OBJ & obj)
	{
	EmitGeneric<kTabSeparatedTabularOutStyle, OBJ>(out, obj);
#   if !defined(NCL_USE_NXS_CONSOLE_OUTPUT) && ! defined(NCL_USE_STD_OUTPUT)
		// in GUI mode we might have plot stream to deal with 
		if (out.GetPlotStreamPtr())
			{
			NxsSaxOutputWrapper nsow(*out.GetPlotStreamPtr(), "entry"); //@temp hard coding the xml here isn't optimal
			Emit<kTabSeparatedTabularOutStyle, OBJ, NxsWritableStream>(*out.GetPlotStreamPtr(), obj);
			}
#   endif
	return out;
	}
#include "ncl/output/nxs_output_stream_wrapper.hpp"
template<class TABLE_OUT, typename DATA>
TABLE_OUT & EmitFirstColumn(TABLE_OUT & out, const DATA & d);
	
NxsWritableStream & NextCol(NxsWritableStream &);
NxsWritableStream & EndFirstRow(NxsWritableStream &);
NxsWritableStream & EndRow(NxsWritableStream &);
NxsOutputStreamWrapper & NextCol(NxsOutputStreamWrapper &);
NxsOutputStreamWrapper & EndFirstRow(NxsOutputStreamWrapper &);
NxsOutputStreamWrapper & EndRow(NxsOutputStreamWrapper &);
std::ostream & NextCol(std::ostream &);
std::ostream & EndFirstRow(std::ostream &);
std::ostream & EndRow(std::ostream &);
inline NxsWritableStream & NextCol(NxsWritableStream &out)
	{
	return out << '\t';
	}
inline NxsWritableStream & EndFirstRow(NxsWritableStream &out)
	{
	return out << '\n';
	}
inline NxsWritableStream & EndRow(NxsWritableStream &out)
	{
	return out << '\n';
	}
inline NxsOutputStreamWrapper & NextCol(NxsOutputStreamWrapper &out)
	{
	return out << '\t';
	}
inline NxsOutputStreamWrapper & EndFirstRow(NxsOutputStreamWrapper &out)
	{
	return out << '\n';
	}
inline NxsOutputStreamWrapper & EndRow(NxsOutputStreamWrapper &out)
	{
	return out << '\n';
	}

inline  std::ostream & NextCol(std::ostream &out)
	{
	return out << '\t';
	}
	
inline  std::ostream & EndFirstRow(std::ostream &out)
	{
	return out << '\n';
	}
	
inline std::ostream & EndRow(std::ostream &out)
	{
	return out << '\n';
	}

template<class TABLE_OUT, typename DATA>
inline TABLE_OUT & EmitFirstColumn(TABLE_OUT & out, const DATA & d)
	{
	out << d;
	return out;
	}

#endif

