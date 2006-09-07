#ifndef PHO_MULTILINE_BOUND_H
#define PHO_MULTILINE_BOUND_H
#include <list>

class Line
	{
	public:
		class XVerticalLine{};
		typedef std::pair<double, double> XYPoint;
		double 	yIntercept;
		double 	slope;
		void   *userInfo;
		
		Line(double linesYIntercept, double linesSlope)
			:yIntercept(linesYIntercept),
			slope(linesSlope)
			{}
			
		Line(const XYPoint &onePt, const XYPoint &anotherPt);

		double operator()(double x) const
			{
			return yIntercept + slope*x;
			}
		
		XYPoint FindIntersectionPoint(const Line &) const;
	};

class LineSegment
	{
	public:
		Line unboundedLine;
		double minX, maxX;
		double minY, maxY;
		
		LineSegment(const Line &l, double min, double max);
		double operator()(double x) const
			{
			NXS_ASSERT(minX <= x && x <= maxX);
			return unboundedLine(x);
			}
		
		double AtLeft() const
			{
			return (*this)(minX);
			}
		double AtRight() const
			{
			return (*this)(maxX);
			}
		double GetMinY() const {return minY;}
		double GetMaxY() const {return maxY;}
		double GetMinX() const {return minX;}
		double GetMaxX() const {return maxX;}
		void 	SetMaxX(double newMaxX);
		void 	SetMinX(double newMinX);
		std::pair<double, double> GetMaxPoint() const
			{
			return std::pair<double, double>((unboundedLine.slope > 0.0 ? GetMaxX() : GetMinX()), GetMaxY());
			}
		std::pair<double, double> GetMinPoint() const
			{
			return std::pair<double, double>((unboundedLine.slope > 0.0 ? GetMinX() : GetMaxX()), GetMinY());
			}
		
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Class that holds line segments used to represent where the y value of the line represents a bound for some quantity
|	(either an upper or lower bound).  Useful when a bound depends on x, and there are several conditions (differnent
|	linear relationships) that determine the bound.
|	Because the bounds are linear this class only makes sense if there is a range of x values that are legal (hence min and max in the constructor)
|
|	If we are using AddLowerBound() this function helps us keep track of the highest lower bound.
*/
class MultilineBound
	{
	public:
		MultilineBound(double min_x, double max_x);
		
		bool AddUpperBound(const Line &lineToAdd);
		bool AddLowerBound(const Line &lineToAdd);
		std::pair<double, double> GetMaxPoint() const
			{
			return maxPt;
			}
		std::pair<double, double> GetMinPoint() const
			{
			return minPt;
			}
		double GetMaxY() const
			{
			return maxPt.second;
			}
		double GetMinY() const
			{
			return minPt.second;
			}
		double operator()(double x) const;
		
	private:
		enum TypeOfBound 
			{
			kUnknownBoundType,
			kUpperBoundType,
			kLowerBoundType
			};	
		double minX, maxX;
		std::pair<double, double> maxPt; 
		std::pair<double, double> minPt; 
		TypeOfBound typeOfBound; //should really make separate classes for lower and upper (ideally template bool parameter determine which type
		typedef std::list<LineSegment> ListOfLines;
		typedef ListOfLines::iterator LineIter;
		// list of lines sorted by the range of x that they apply to 
		ListOfLines boundSegments;
		
		bool 	AddCroppedUpperBound(const LineSegment &lineSeg);
		bool 	AddCroppedLowerBound(const LineSegment &lineSeg);
		void	FindMinAndMax();
	};

#endif
