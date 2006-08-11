//#include "phycas/force_include.h"
#include "ncl/misc/algorithm_extensions.hpp"
#include "phycas/misc/multiline_bound.hpp"
using std::pair;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::fabs;
#endif

MultilineBound::MultilineBound(double min_x, double max_x)
	:minX(min_x), 
	maxX(max_x), 
	maxPt(-DBL_MAX, -DBL_MAX), 
	minPt(-DBL_MAX, -DBL_MAX),
	typeOfBound(kLowerBoundType) //need to set new defaults for maxPt and minPt if we want something other than upperbounds
	{
	assert(minX > -DBL_MAX && maxX < DBL_MAX);
	}

LineSegment::LineSegment(const Line &l, double min, double max)
	:
	unboundedLine(l), 
	minX(min), 
	maxX(max)
	{
	minY = AtLeft();
	maxY = AtRight();
	pho_sort(minY, maxY);
	}
		
	
bool MultilineBound::AddUpperBound(const Line &lineToAdd) //NOT DEBUGGED
	{
	assert(0);//NOT DEBUGGED YET!!!
	assert(typeOfBound == kUnknownBoundType || typeOfBound == kUpperBoundType);
	typeOfBound = kUpperBoundType;
	LineSegment lineSeg(lineToAdd, minX, maxX);
	return AddCroppedUpperBound(lineSeg);
	}
	
bool MultilineBound::AddCroppedUpperBound(const LineSegment &lineSeg)
	{
	assert(0);//NOT DEBUGGED YET!!!
	double newLineMin = lineSeg.AtLeft();
	double newLineMax = lineSeg.AtRight();
	pho_sort(newLineMin, newLineMax);
	if (newLineMin > GetMaxY()) //	if the new line is implies a higher upper bound that the current upper bound over the entire range of x, then we don't need to add this line.
		return false;
	if (newLineMax < GetMinY()) //	if the new line is implies a lower upper bound that the current lower bound over the entire range of x, then we can replace all of the boundSegments
		{
		boundSegments.clear();
		boundSegments.push_back(lineSeg);
		maxPt = lineSeg.GetMaxPoint();
		minPt = lineSeg.GetMinPoint();
		return true;
		}
	bool newLineWasAdded = false;
	bool previousLineIsTheNewOne = false;
	LineIter prevLine = boundSegments.begin();
	// remember to call FindMinAndMax before returning?
	
	for (LineIter bsIt = boundSegments.begin(); bsIt != boundSegments.end();)
		{
		if (lineSeg.GetMaxX() < bsIt->GetMinX())	//because boundSegments is sorted, we are done.
			{
			FindMinAndMax();
			return newLineWasAdded;
			}
		if (lineSeg.GetMinX() < bsIt->GetMaxX())
			{
			double currMaxX = bsIt->GetMaxX();
			double currMinX = bsIt->GetMinX();
			bool lowerAtLeft = (lineSeg(currMinX) < bsIt->AtLeft());
			bool lowerAtRight = (lineSeg(currMaxX) < bsIt->AtRight());
			if (lowerAtLeft)
				{
				if (lowerAtRight)
					bsIt = boundSegments.erase(bsIt);
				else
					{
					std::pair<double, double> intersection = bsIt->unboundedLine.FindIntersectionPoint(lineSeg.unboundedLine);
					currMaxX = intersection.first;
					}
				if (!previousLineIsTheNewOne)
					{
					newLineWasAdded = true;
					NXS_ASSERT(currMinX >= minX && currMinX <= maxX);
					NXS_ASSERT(currMaxX >= minX && currMaxX <= maxX);
					LineSegment lineToAdd(lineSeg.unboundedLine, currMinX, currMaxX);
					//@@@ check does insert return iterator to the added element, or one after?
					prevLine = boundSegments.insert(bsIt, lineToAdd);
					}
				else
					{
					NXS_ASSERT(currMaxX >= minX && currMaxX <= maxX);
					prevLine->SetMaxX(currMaxX);
					}
				if (!lowerAtRight)
					{
					NXS_ASSERT(currMaxX >= minX && currMaxX <= maxX);
					bsIt->SetMinX(currMaxX);
					++bsIt;
					}
				}
			else if (lowerAtRight)
				{
				std::pair<double, double> intersection = bsIt->unboundedLine.FindIntersectionPoint(lineSeg.unboundedLine);
				currMinX = intersection.first;
				bsIt->SetMaxX(currMinX);
				++bsIt;
				NXS_ASSERT(currMinX >= minX && currMinX <= maxX);
				NXS_ASSERT(currMaxX >= minX && currMaxX <= maxX);
					
				LineSegment lineToAdd(lineSeg.unboundedLine, currMinX, currMaxX);
				prevLine = boundSegments.insert(bsIt, lineToAdd);
				bsIt = prevLine;
				++bsIt;
				newLineWasAdded = true;
				}
			else
				++bsIt;
			previousLineIsTheNewOne = lowerAtRight;
			}
		}
	FindMinAndMax();
	return newLineWasAdded;
	}

bool MultilineBound::AddLowerBound(const Line &lineToAdd) //NOT DEBUGGED
	{
	assert(typeOfBound == kUnknownBoundType || typeOfBound == kLowerBoundType);
	typeOfBound = kLowerBoundType;
	LineSegment lineSeg(lineToAdd, minX, maxX);
	return AddCroppedLowerBound(lineSeg);
	}

double MultilineBound::operator()(double x) const
	{
	for (ListOfLines::const_iterator bsIt = boundSegments.begin(); bsIt != boundSegments.end(); ++bsIt)
		{
		if (x >=  bsIt->GetMinX() && x <=  bsIt->GetMaxX())
			return  (*bsIt)(x);
		}
	return DBL_MAX;
	}

	
bool MultilineBound::AddCroppedLowerBound(const LineSegment &lineSeg)
	{
	double newLineMin = lineSeg.AtLeft();
	double newLineMax = lineSeg.AtRight();
	pho_sort(newLineMin, newLineMax);
	if (newLineMax < GetMinY()) //	if the new line is implies a lower lower bound that the current lower bound over the entire range of x, then we don't need to add this line.
		return false;
	if (newLineMin > GetMaxY()) //	if the new line is implies a higher upper bound that the bounds over the entire range of x, then we can replace all of the boundSegments
		{
		boundSegments.clear();
		boundSegments.push_back(lineSeg);
		maxPt = lineSeg.GetMaxPoint();
		minPt = lineSeg.GetMinPoint();
		return true;
		}
	bool newLineWasAdded = false;
	bool previousLineIsTheNewOne = false;
	LineIter prevLine = boundSegments.begin();
	// remember to call FindMinAndMax before returning?
	
	for (LineIter bsIt = boundSegments.begin(); bsIt != boundSegments.end();)
		{
		if (lineSeg.GetMaxX() < bsIt->GetMinX())	//because boundSegments is sorted, we are done.
			{
			if (newLineWasAdded)
				FindMinAndMax();
			return newLineWasAdded;
			}
		if (lineSeg.GetMinX() < bsIt->GetMaxX())
			{
			double currMaxX = bsIt->GetMaxX();
			double currMinX = bsIt->GetMinX();
			bool higherAtLeft = (lineSeg(currMinX) > bsIt->AtLeft());
			bool higherAtRight = (lineSeg(currMaxX) > bsIt->AtRight());
			if (higherAtLeft)
				{
				if (higherAtRight)
					bsIt = boundSegments.erase(bsIt);
				else
					{
					pair<double, double> intersection = bsIt->unboundedLine.FindIntersectionPoint(lineSeg.unboundedLine);
					currMaxX = intersection.first;
					}
				if (!previousLineIsTheNewOne)
					{
					newLineWasAdded = true;
					NXS_ASSERT(currMinX >= minX && currMinX <= maxX);
					NXS_ASSERT(currMaxX >= minX && currMaxX <= maxX);
					LineSegment lineToAdd(lineSeg.unboundedLine, currMinX, currMaxX);
					//@@@ check does insert return iterator to the added element, or one after?
					prevLine = boundSegments.insert(bsIt, lineToAdd);
					}
				else
					{
					NXS_ASSERT(currMaxX >= minX && currMaxX <= maxX);
					prevLine->SetMaxX(currMaxX);
					}
				if (!higherAtRight)
					{
					NXS_ASSERT(currMaxX >= minX && currMaxX <= maxX);
					bsIt->SetMinX(currMaxX);
					++bsIt;
					}
				}
			else if (higherAtRight)
				{
				pair<double, double> intersection = bsIt->unboundedLine.FindIntersectionPoint(lineSeg.unboundedLine);
				currMinX = intersection.first;
				bsIt->SetMaxX(currMinX);
				++bsIt;
				NXS_ASSERT(currMinX >= minX && currMinX <= maxX);
				NXS_ASSERT(currMaxX >= minX && currMaxX <= maxX);
				LineSegment lineToAdd(lineSeg.unboundedLine, currMinX, currMaxX);
				prevLine = boundSegments.insert(bsIt, lineToAdd);
				bsIt = prevLine;
				++bsIt;
				newLineWasAdded = true;
				}
			else
				++bsIt;
			previousLineIsTheNewOne = higherAtRight;
			}
		}
	if (newLineWasAdded)
		FindMinAndMax();
	return newLineWasAdded;
	}

void LineSegment::SetMaxX(double newMaxX)
	{
	maxX = newMaxX;
	if (unboundedLine.slope > 0)
		maxY = AtRight();
	else
		minY = AtRight();
	}

void LineSegment::SetMinX(double newMinX)
	{
	minX = newMinX;
	if (unboundedLine.slope > 0)
		minY = AtLeft();
	else
		maxY = AtLeft();
	}


/*----------------------------------------------------------------------------------------------------------------------
|	returns pair(x, y) that decribes the point that this and `otherLine' intersect.
|	returns pair(DBL_MAX, DBL_MAX) if the lines are (essentially) parallel
*/ 
pair<double, double> Line::FindIntersectionPoint(const Line &otherLine) const
	{
	XYPoint retVal = XYPoint(DBL_MAX, DBL_MAX);
	double denom = slope - otherLine.slope;
	if (fabs(denom) > DBL_EPSILON)
		{
		retVal.first = (otherLine.yIntercept - yIntercept) / denom;
		retVal.second = (*this)(retVal.first);
		}
	return retVal;
	}

Line::Line(const XYPoint &onePt, const XYPoint &anotherPt)
	{
	if (fabs(onePt.first-anotherPt.first) < DBL_EPSILON)
		{
		assert(0);
		throw XVerticalLine();
		}
	slope = (onePt.second-anotherPt.second)/(onePt.first-anotherPt.first);
	yIntercept = onePt.second - onePt.first * slope;
	}

void MultilineBound::FindMinAndMax()
	{
	maxPt = Line::XYPoint(-DBL_MAX,-DBL_MAX); 
	minPt =  Line::XYPoint(DBL_MAX, DBL_MAX);
	for (ListOfLines::const_iterator bsIt = boundSegments.begin(); bsIt != boundSegments.end(); ++bsIt)
		{
		Line::XYPoint currSegPt = bsIt->GetMaxPoint();
		if (currSegPt.second > maxPt.second)
			maxPt = currSegPt;
		currSegPt = bsIt->GetMinPoint();
		if (currSegPt.second < minPt.second)
			minPt = currSegPt;
		}
	}

