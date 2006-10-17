/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

#ifndef NCL_SIMPLESTATS_H
#define NCL_SIMPLESTATS_H

/*----------------------------------------------------------------------------------------------------------------------
|	Implements the NxsSimpleStats class, which allows easy calculation of simple summary statistics.
|>
|  Example:
|		NxsSimpleStats stats;
|		stats.Add(1.0);
|		stats.Add(2.0);
|		stats.Add(3.0);
|		stats.Add(4.0);
|		stats.Add(5.0);
|		long   n         = stats.N();        // n = 5
|		double mean      = stats.Mean();     // mean = 3.000
|		double variance  = stats.Variance(); // variance = 2.500
|		double stddev    = stats.StdDev();   // stddev = 1.581
|		double cv        = stats.CV();       // cv = 0.527
|>
*/
class NxsSimpleStats
	{
	double minimum;
	double maximum;
	double sum;
	double lsum; // lower sum used only for ratios
	double sumSq;
	double lsumSq; // lower sum of squares used only for ratios
	double mean;
	double variance;
	double stddev;
	double cv;
	long n;
	bool dirty;
	void CalcMean();
	void CalcVariance();
	void CalcStdDev() { stddev = (variance <= 0.0 ? 0.0 : std::sqrt(variance)); }
	void CalcCV() { cv = (mean == 0.0 ? 0.0 : stddev / mean); }
	void CalcAll();
	public:
		NxsSimpleStats();

		void Add(const double);
		void Add(const NxsRatio);
		long N() { return n; }
		double Sum() { return sum; }
		double SumSq() { return sumSq; }
		double Mean();
		double Variance();
		double StdDev();
		double CV();
		double Minimum();
		double Maximum();
		NxsSimpleStats& operator +=(const double v);
		NxsSimpleStats& operator +=(const NxsRatio r);
		friend std::ostream& operator<<(std::ostream&, NxsSimpleStats&);
	};

inline NxsSimpleStats::NxsSimpleStats()
	{
	minimum		= 0.0;
	maximum		= 0.0;
	sum			= 0.0;
	lsum		= 0.0;
	sumSq		= 0.0;
	lsumSq		= 0.0;
	mean		= 0.0;
	variance	= 0.0;
	stddev		= 0.0;
	cv			= 0.0;
	n			= 0L;
	dirty		= false;
	}

inline void NxsSimpleStats::CalcMean()
	{
	if (lsum > 0.0)
		mean = sum / lsum;
	else
		mean = sum / static_cast<double>(n);
	}

inline void NxsSimpleStats::CalcVariance()
	{
	double s = sum;
	double ss = sumSq;
	if (lsum > 0.0) {
		s /= lsum;
		s *= static_cast<double>(n);
		ss /= lsumSq;
		ss *= static_cast<double>(n);
	}
	variance = (ss - s * s / static_cast<double>(n)) / 
             static_cast<double>(n - 1L);
	}

inline void NxsSimpleStats::CalcAll()
	{
	if (dirty && n > 0L)
		CalcMean();
	if (dirty && n > 1L) {
		CalcMean();
		CalcVariance();
		CalcStdDev();
		CalcCV();
	}
	dirty = false;
	}

inline double NxsSimpleStats::Mean()
	{
	CalcAll();
	return mean;
	}

inline double NxsSimpleStats::Variance()
	{
	CalcAll();
	return variance;
	}

inline double NxsSimpleStats::StdDev()
	{
	CalcAll();
	return stddev;
	}

inline double NxsSimpleStats::CV()
	{
	CalcAll();
	return cv;
	}

inline double NxsSimpleStats::Minimum()
	{
	return minimum;
	}

inline double NxsSimpleStats::Maximum()
	{
	return maximum;
	}

inline void NxsSimpleStats::Add(const double v)
	{
	if (n == 0L)
		{
		minimum = v;
		maximum = v;
		}
	else
		{
		if (v < minimum)
			minimum = v;
		if (v > maximum)
			maximum = v;
		}
	sum += v;
	sumSq += v*v;
	++n;
	dirty = true;
	}

inline void NxsSimpleStats::Add(const NxsRatio r)
	{
	if (n == 0L)
		{
		minimum = r.top / r.bottom;
		maximum = r.top / r.bottom;
		}
	else
		{
		double rr = r.top / r.bottom;
		if (rr < minimum)
			minimum = rr;
		if (rr > maximum)
			maximum = rr;
		}
	sum		+= r.top;
	lsum	+= r.bottom;
	sumSq	+= r.top * r.top;
	lsumSq	+= r.bottom * r.bottom;
	++n;
	dirty = true;
	}

inline NxsSimpleStats& NxsSimpleStats::operator +=(const double v)
	{
	Add(v);
	return *this;
	}

inline NxsSimpleStats& NxsSimpleStats::operator +=(const NxsRatio r)
	{
	Add(r);
	return *this;
	}

inline std::ostream& operator<<(std::ostream& out, NxsSimpleStats& stat)
	{
	out << "\n\tN         = " << stat.N();
	out << "\n\tmean      = " << stat.Mean();
	out << "\n\tvariance  = " << stat.Variance();
	out << "\n\tstd. dev. = " << stat.StdDev();
	out << "\n\tCV        = " << stat.CV();
	return out;
	}

#endif
