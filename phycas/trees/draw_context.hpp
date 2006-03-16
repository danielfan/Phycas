#ifndef PHO_DRAWCONTEXT_H
#define PHO_DRAWCONTEXT_H

#include "ncl/output/nxs_user_query.hpp"
#include "ncl/output/nxs_output.hpp"

class PhoFont
	{
	};

class PhoPoint
	{
	double x, y;
	};

class PhoCanvas
	{
	public:
		//PhoPoint GetSize();
		void DrawLine(PhoPoint fromP, PhoPoint toP);
		void WriteText(const char *, PhoPoint StartPoint, PhoFont f);
	};
	
class PhoWindow	// any window
	{
	public:
		PhoCanvas 	*GetCanvas();
	};

class CommandItem;
class PhoControl
	{
	CommandItem *cmdReader;
	};
	
class PhoDialog
	{
	std::vector<PhoControl> controls;
	};

class TreeNode;
class DrawContext	
	{
		PhoFont			currentFont;
	public :
		// Services
		//
		virtual ~DrawContext(){}
		virtual unsigned 	GetStringWidth(const std::string &) const = 0;
		virtual void 		DrawTree(const TreeNode *ingroupSubRoot, bool showAsRooted) = 0;
		virtual double 		GetHeightScaler(double max_y) = 0;	// also tells DrawContext the height
		virtual double		GetPlotWidth() const = 0;
	};

class GUIDrawContext: public DrawContext
	{
	PhoCanvas	*canvas; 
		virtual ~GUIDrawContext(){}
		unsigned 	GetStringWidth(const std::string &) const;
		void 		DrawTree(const TreeNode *ingroupSubRoot, bool showAsRooted);
		double 		GetHeightScaler(double max_y);
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Handles output to portable versions 
|	into line-drawing characters for purposes of outputting ASCII representations of trees. The string representing a
|	line in such a representation is coded so that a vertical line is represented as (Nub_north | Nub_south), a 
|	horizontal line as (Nub_East | Nub_West), the 'left shoulder' of a tree with tips on the right (where down is
|	toward the root) would be (Nub_North | Nub_East), etc.
*/
class ASCIIDrawContext: public DrawContext
	{
	public:	// temporary data should be private
		enum { Nub_North = 0x01, Nub_East = 0x02, Nub_South = 0x04, Nub_West = 0x08 };

		STATIC_DATA const unsigned char codeToAscii[];
		STATIC_DATA const unsigned char codeToPaupMonaco[];
		
		virtual ~ASCIIDrawContext(){}
		void 			TranslateLine(std::string &lineBuffer);
		bool 			compressed;
		bool			usePaupMonaco;
		mutable unsigned windowWidth;
		unsigned		heightOfNextGraphic;
		
						ASCIIDrawContext(NxsOutputStream &os);

		unsigned 		GetStringWidth(const std::string &s) const {return (unsigned)s.length();};
		void			DrawTree(const TreeNode *ingroupSubRoot, bool showAsRooted);
		double 			GetHeightScaler(double max_y);
		double			GetPlotWidth() const;

	private:
		NxsOutputStream &outStream;
	};

#endif
