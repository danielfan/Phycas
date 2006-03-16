#if !defined(PHYCAS_TREE_OUTPUT_HPP)
#define PHYCAS_TREE_OUTPUT_HPP
#include "ncl/output/temp_basic_output_operators.hpp"
#include "phycas/trees/tree.hpp"
#include "phycas/trees/draw_context.hpp"

template<>
class GenericPrinterClass<kGraphicalTreeOutStyle, Tree, NxsOutputStream> //@ temp! use of typedef NxsOutputStream thwarts any rational dispatching based on output stream type
	{
	public:
		GenericPrinterClass(NxsOutputStream & outStream, const Tree & tree)
			{
			ASCIIDrawContext ascii(outStream);
			tree.Draw(ascii);
			}
	};

template<class OUT_STREAM>
class GenericPrinterClass<kNewickTreeOutStyle, Tree, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const Tree & tree)
			{
			std::string s;
			tree.AppendNewickRepresentation(s);
			outStream << s;
			}
	};

#endif
