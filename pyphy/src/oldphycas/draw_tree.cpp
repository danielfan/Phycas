//#include "phycas/force_include.h"
#include "phycas/trees/tree_iterator.hpp"
#include "phycas/trees/draw_context.hpp"
#include "phycas/misc/multiline_bound.hpp"
#include "ncl/misc/string_extensions.hpp"
using std::string;

// 		this file uses several functions to alter a node's plot data that are only likely to be used within this file
//		rather than clutter TreeNode's class declaration with these functions, they have been placed here.
//
void 		CalcWidthScaler(TreeNode *focalNd, const DrawContext *dc, double plotWidth, double &scaler);
void 		CalcWidthScalerBounds(TreeNode *focalNd, const DrawContext *dc, double plotWidth, MultilineBound &scaler, double zeroAdditionalX, double oneAdditionalX);
unsigned 	CreateDefaultCoordsRecursively(TreeNode *focalNd, const DrawContext &dc, double *maxNameLen, double &max_y);
void 		SetUnscaledCoordsRecursively(TreeNode *focalNd, double x, double &max_y);
void 		RescaleDefaultCoords(TreeNode *focalNd, const double widthToPrint, const double yscaler);
void 		FillLine(const TreeNode *focalNd, const unsigned y, const unsigned width, string &lineBuffer);

/*----------------------------------------------------------------------------------------------------------------------
|	This is the old Tree::Draw routine from when we were rooting at a polytomy.  Now Draw makes a copy of the tree
|	roots it at the appropriate place and calls this function.
|
*/
void Tree::DrawTreeWithTempRoot(DrawContext &dc)
	{
	TreeNode *tempRoot = GetRoot();
	tip_iterator a(tempRoot);
	++a;
	NXS_ASSERT(tempRoot->GetLeftChild() != NULL);
	pre_iterator beginIt = pre_begin();
	const const_pre_iterator endIt = pre_end();
	const tip_iterator endOfTips = leaves_end();
	max_y = 0.0;
	double	widthScaler = DBL_MAX;
	tempRoot->SetFltEdgeLen(0.0);
	if (HasEdgeLens())
		{
		SetUnscaledCoordsRecursively(tempRoot, 0.0, max_y);
		for (beginIt = pre_begin(); beginIt != endIt; ++beginIt)
			CalcWidthScaler(&*beginIt, &dc, dc.GetPlotWidth(), widthScaler);	//drawer provides name width 
		if (widthScaler == DBL_MAX)
			throw XBadTree("\nCan't draw tree with all-zero branch lengths");// temp show error msg
		double heightScaler = dc.GetHeightScaler(max_y);
		for (beginIt = pre_begin(); beginIt != endIt; ++beginIt)
			beginIt->ScaleCoords(widthScaler, heightScaler);
		}
	else
		{
		tempRoot->RefreshPlotData();
		tempRoot->SetX(0.0);
		double maxNameLen = 0.0;
		CreateDefaultCoordsRecursively(tempRoot->GetLeftChild(), dc, &maxNameLen, max_y); 
		if (maxNameLen >= dc.GetPlotWidth()) 
			throw XBadTree("\nCan't draw tree because the width is too low");// temp show error msg
		double heightScaler = dc.GetHeightScaler(max_y);
		for (beginIt = pre_begin(); beginIt != endIt; ++beginIt)
			RescaleDefaultCoords(&*beginIt, (double) dc.GetPlotWidth() - maxNameLen, heightScaler);
		}
	tempRoot->SetX(0.0);
	tempRoot->GetLeftChild()->SetX(0.0);
	dc.DrawTree(tempRoot, IsRooted());
	}
	

void ASCIIDrawContext::DrawTree(
  const TreeNode *treeRoot, 	/* pointer to the root of the ingroup or outgroup sub clade which ever is the child node of the other */
  bool showAsRooted)
	{
	const_pre_iterator nd = treeRoot->clade_pre_begin();
	const const_pre_iterator endIt = treeRoot->clade_pre_end();
	if (!showAsRooted)
		++nd;
	const const_pre_iterator beginNdIt = nd;
	
	string lineBuffer;
	lineBuffer.reserve(windowWidth);
	for (unsigned y = 0; y < heightOfNextGraphic; ++y)
		{
		lineBuffer.clear();
		lineBuffer.append(windowWidth, '\0');
		for (nd = beginNdIt; nd != endIt; ++nd)
			FillLine(&*nd, y, windowWidth, lineBuffer);
		TranslateLine(lineBuffer);
		outStream << lineBuffer.c_str() << '\n';
		}
	outStream << ncl::endl;
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a tree using ASCII-style graphics. Assumes root is non-NULL and has at least one child
*/
void Tree::Draw(
  DrawContext &dc)	const	/* output stream to which tree should be sent */
	{
	//	To accommodate the new rooting scheme, we will copy the tree and root it at a fake tip so t
	//
	assert(&dc != NULL);
	Tree treeToDraw(*this);
	TreeNode *tempRoot = NULL; //will be allocated independently, but will then be attached to treeToDraw, and destroyed when it goes out of scope
	assert(treeToDraw.GetRoot() != NULL);
	if (treeToDraw.IsRooted())
		{
		NXS_ASSERT(0);//@need to add support for rooted trees
		}
	else
		{
		const bool userSpecifiedOutGroup = false;	//@need to add support for displaying outgroup taxa one on side of the tree
		if (!userSpecifiedOutGroup)
			{
			//	default out group is the lowest numbered taxon. which should be the root of the tree
			//
			TreeNode *oldRoot = treeToDraw.GetRoot();
			if (oldRoot == NULL)
				{
				assert(0);
				return;//@@@ should throw some exception
				}
			
			TreeNode *ingroupMRCA = oldRoot->GetLeftChild();
			if (ingroupMRCA == NULL)
				{
				assert(0);
				return;//@@@ should throw some exception
				}
			assert(ingroupMRCA->GetRightSib() == NULL);
			ingroupMRCA = ingroupMRCA->FindRightChild();
			tempRoot = treeToDraw.CreateTreeNode(UINT_MAX, 0.0, false);
			ingroupMRCA->AddRightSib(tempRoot);
			}
		else
			{
			//const bool showOutgroupAsMonophyletic = true;	//@ need to add support for displaying outgroup in polytomy or as monophyletic /
		
			NXS_ASSERT(0);
			// 	here we need to find the most inclusive node that is not in the outgroup and attach tempRoot
			}
		}
	NXS_ASSERT(tempRoot != NULL);
	
	treeToDraw.Reroot(tempRoot);
	treeToDraw.DrawTreeWithTempRoot(dc);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes x and y coordinates for each node in preparation for a tree plot, creating the plotData objects as
|	necessary.
|	This is a version of SetUnscaledCoords to be used when the tree doesn't have edge lengths (or the edge lengths are
|	being ignored).  
|	RescaleDefaultCoords, must be called 
|	Leaves are set to (x = 1.0, y = the # of tips before this node (when nodes are traversed in preorder)
|	Internals are set to (x = propor of remaining width on this branch, y = midpoint of childern's y)
*/
unsigned CreateDefaultCoordsRecursively(
  TreeNode *focalNd, 
  const DrawContext &dc,	/* provides lengths of leave labels */
  double *maxNameLen, /* at the end of the recursion, this will be the amount of space needed to leave room for all of the taxa names */
  double &max_y)		/* keeps track of current maximum y-coordinate, which is incremented by 1 for each leaf node encountered */
	{
	focalNd->RefreshPlotData();
	
	if (focalNd->IsShootTip())
		{
		focalNd->SetX(1.0);
		focalNd->SetY(max_y++);
		string n = "  "; 
		n << focalNd->GetName();
		unsigned labelWidth = dc.GetStringWidth(n);
		if (labelWidth > *maxNameLen)
			*maxNameLen = (double) labelWidth;
		return 1;
		}
	
	//	Recurse through the tree and see what the shortest returned branch lenghth is
	//
	unsigned maxGuess = 0;
	TreeNode *rightmost_child = focalNd->GetLeftChild();
	for (;;)
		{
		unsigned tempX = CreateDefaultCoordsRecursively(rightmost_child, dc, maxNameLen, max_y);
		if (tempX > maxGuess)
			maxGuess = tempX;
		if (rightmost_child->GetRightSib() == NULL)
			break;
		rightmost_child = rightmost_child->GetRightSib();
		}
	focalNd->SetY(0.5*(focalNd->GetLeftChild()->GetY() + rightmost_child->GetY()));
	
	if (focalNd->GetParent() == NULL || focalNd->GetParent()->GetParent() == NULL)
		focalNd->SetX(0.0);
	else 
		{
		++maxGuess;
		focalNd->SetX(1.0/(double) maxGuess);
		}
	
	return maxGuess;
	}

void RescaleDefaultCoords(TreeNode *focalNd, const double widthToPrint, const double yscaler)
	{
	if (focalNd->IsShootTip())
		focalNd->SetX(widthToPrint);
	else
		{
		const double parX = (focalNd->GetParent() == NULL ? 0.0 : focalNd->GetParent()->GetX());
		const double newX = parX + (widthToPrint-parX)*(focalNd->GetX());
		focalNd->SetX(newX);
		}
	focalNd->SetY(focalNd->GetY() * yscaler);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes x and y coordinates for each node in preparation for a tree plot, creating the plotData objects as
|	necessary.
*/
void SetUnscaledCoordsRecursively(
  TreeNode *focalNd, 
  double x,			/* the x-coordinate of the parent of this node; x for this node is set to this value plus this node's edge length */
  double &max_y)	/* keeps track of current maximum y-coordinate, which is incremented by 1 for each leaf node encountered */
	{
	focalNd->RefreshPlotData();
	
	x += focalNd->GetFltEdgeLen();
	focalNd->SetX(x);
	
	if (focalNd->IsShootTip())
		focalNd->SetY(max_y++);
	else
		{
		TreeNode *rightmost_child = focalNd->GetLeftChild();
		for (;;)
			{
			SetUnscaledCoordsRecursively(rightmost_child, x, max_y);
			if (rightmost_child->GetRightSib() == NULL)
				break;
			rightmost_child = rightmost_child->GetRightSib();
			}
		focalNd->SetY(0.5*(focalNd->GetLeftChild()->GetY() + rightmost_child->GetY()));
		}
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Scales x and y coordinates of plotData structure and then continues the recursion by visiting lChild and rSib in 
|	that order. Assumes plotData is non-NULL.
*/
void TreeNode::AddToX(
  double additionalX)	/* scaling factor to use for y-coordinates */
	{
	assert(plotData != NULL);
	SetX(GetX() + additionalX);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Scales x and y coordinates of plotData structure. Assumes plotData is non-NULL.
*/
void TreeNode::ScaleCoords(
  double xscaler,	/* scaling factor to use for x-coordinates */
  double yscaler)	/* scaling factor to use for y-coordinates */
	{
	assert(plotData != NULL);
	SetX(GetX() * xscaler);
	SetY(GetY() * yscaler);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Prepares to output line of tree by storing codes for edge-segments and overlaying taxon names. Used for drawing 
|	trees on an ASCII canvas. Continues recursion by visiting GetLeftChild() and rSib in that order.
|	Requires that plotData->(x and y) have been initialized in this node and its parent, children and siblings.
*/
void FillLine(
  const TreeNode *focalNd, 
  const unsigned y,				/* current line */
  const unsigned ,			/* number of characters in lineBuffer */
  string &lineBuffer)	/* the buffer in which drawing occurs */
	{
	unsigned	x = (unsigned)focalNd->GetX();
	
	if ((y == (unsigned)focalNd->GetY()) && !focalNd->IsRoot())
		{
		// Fill in my branch
		if (x > 0)
			{
			unsigned	x_par = (unsigned)focalNd->GetParent()->GetX();

			if (lineBuffer[x_par] < 16)
				lineBuffer[x_par] |= ASCIIDrawContext::Nub_East;
			for (unsigned col = x_par + 1; col < x; ++col)
				{
				if (lineBuffer[col] < 16)
					lineBuffer[col] |= (ASCIIDrawContext::Nub_East | ASCIIDrawContext::Nub_West);
				}
			if (lineBuffer[x] < 16)
				lineBuffer[x] |= ASCIIDrawContext::Nub_West;
			}
		
		if (focalNd->IsShootTip()) 
			{
			const string n(focalNd->GetName());
			lineBuffer.replace(x + 2, n.length(), n);		// overlay taxon name
			}
#if 1
		//POL added this section Dec. 27, 2003, in order to show posterior clade probabilities 
		// for interior edges if these probabilities have been specified. Also added the 
		// 'if (lineBuffer[x] < 16)' checks to ensure that percentages laid down below are 
		// not modified by ORing with ASCIIDrawContext::Nub_* constants.
		//
		else 
			{
			if (focalNd->IsPlotData() && focalNd->GetPct() >= 0.0)
				{
				const string pctstr = MakeStrPrintF("%d", (unsigned)(0.5 + focalNd->GetPct()));
				lineBuffer.replace(x + 2, pctstr.length(), pctstr);		// overlay information about internal node
				}
			}
#endif

		// Uncomment two lines below to always show internal nodes
		//else 
		//	lineBuffer[x] |= (ASCIIDrawContext::Nub_North | ASCIIDrawContext::Nub_South);
		}
		
	if ((lineBuffer[x] < 16) && !(focalNd->IsShootTip()))
		{
		unsigned y1 = (unsigned) focalNd->GetLeftChild()->GetY();
		unsigned y2 = (unsigned) focalNd->FindRightChild()->GetY();
		if ((y >= y1) && (y <= y2))
			{
			if (y > y1)
				lineBuffer[x] |= ASCIIDrawContext::Nub_North;
			if (y < y2)
				lineBuffer[x] |= ASCIIDrawContext::Nub_South;
			}
		}
	}

	
/*----------------------------------------------------------------------------------------------------------------------
|	For leaves, computes the scaling factor that will just place the right edge of the leaf label at the right edge of 
|	the plot area. If the scaling factor for any given leaf is found to be smaller than scaler, then scaler is set to
|	this leaf's scaling factor. After traversing the tree, scaler will be the correct scaling factor to use in scaling
|	the edge lengths of the tree so that there is no wasted space on the right side (the tree is always drawn so that
|	it lies on its side with taxon labels to the right). If this scaling were not done, wasted space shows up on the 
|	right side if the longest leaf label is not attached to the longest edge.
*/
void CalcWidthScaler(
  TreeNode *focalNd, 
  const DrawContext *dc,	/* provides lengths of leave labels */
  double plotWidth,			/* width of plot area */
  double &scaler) 			/* overall scaling factor for edge lengths */			
	{
	assert(focalNd->plotData != NULL);
	if (focalNd->IsShootTip() && focalNd->GetX() > 0.0)
		{
		// Allow for space between name and branch tip (2 spaces because cast 
		// from double to unsigned bumps tips up to next "pixel")
		//
		string n = "  "; 
		n << focalNd->GetName();
		unsigned labelWidth = dc->GetStringWidth(n);
		double x = (plotWidth - labelWidth)/focalNd->GetX();
		if (x < scaler)
			scaler = x;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	We need to find the max scaling factor that lets all of the labels fit on the tree.
|	This is made difficult by the fact that we don't know how much of the will be added to this side of the tree.
|
|	(widthToPrint - widthOfName) = ScalingFactor * (CurrentNodesX + proporBasalBranch*basalBrLen)
|	So there is a linear relationship between the inverse of the scaling factor and the amount of the basal branch
|	added to this half of the tree:
|	
|	InvScaleFactor = (CurrentNodesX + proporBasalBranch) / (widthToPrint - widthOfName)
|
|	So we want to find the lowest possible InvScaleFactor
|	Each Node adds its lowerbound on InvScaleFactor and we take the calling function (Tree::Draw takes the lowest
|	point that is consistent with all of the bounds.
*/
void CalcWidthScalerBounds(
  TreeNode *focalNd, 
  const DrawContext *dc,	/* provides lengths of leave labels */
  double plotWidth,			/* width of plot area */
  MultilineBound &allBounds, 			/* holds all of the constraints on the scaling factor */
  double minAdded, 				/* amount added to this node's x if 0% of the the basal branch is added to the outgroup side of the root */
  double deltaAddedAtOne)		/* minAdded + deltaAddedAtOne is amount added to this node's x if 100% of the the basal branch is added to the outgroup side of the root */			
	{
	if (focalNd->IsLeafOrRoot())
		{
		assert(focalNd->plotData != NULL);
		// Allow for space between name and branch tip (2 spaces because cast 
		// from double to unsigned bumps tips up to next "pixel")
		//
		string n = "  "; 
		n << focalNd->GetName();
		unsigned labelWidth = dc->GetStringWidth(n);

		double y0 = (minAdded + focalNd->GetX())/(plotWidth - labelWidth); 
		double y1 = (minAdded + deltaAddedAtOne + focalNd->GetX())/(plotWidth - labelWidth);
		Line thisNodesBound(Line::XYPoint(0.0, y0),Line::XYPoint(1.0, y1));
		allBounds.AddLowerBound(thisNodesBound);
		}
	}

