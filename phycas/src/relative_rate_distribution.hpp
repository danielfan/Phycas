/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2010 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

#if !defined(RELATIVE_RATE_DISTRIBUTION_HPP)
#define RELATIVE_RATE_DISTRIBUTION_HPP

#if defined(_MSC_VER)
#	pragma warning(disable: 4267)	// warning about loss of data when converting size_t to int
#endif

#include <cmath>
#include "ncl/nxsdefs.h"

#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>

#include "phycas/src/probability_distribution.hpp"

namespace phycas
{

	
/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	This distribution is essentially a Dirichlet distribution modified so that the parameters have mean 1 rather than sum 1. It is useful as a prior for relative
|	rates (e.g. subset relative rates in a partitioned model).
*/
class RelativeRateDistribution : public DirichletDistribution
	{
	public:
											RelativeRateDistribution();
											RelativeRateDistribution(const std::vector<double> & params);
                        					RelativeRateDistribution(const RelativeRateDistribution & other);
											~RelativeRateDistribution() {}
	
        RelativeRateDistribution * 			cloneAndSetLot(Lot * other) const;
        RelativeRateDistribution * 			Clone() const;

		std::string 						GetDistributionName() const;
		std::string 						GetDistributionDescription() const;
		std::string 						GetDescriptionForPython() const;
		std::vector<double>					GetMean() const;
		std::vector<double> 				GetVar() const;
		std::vector<double> 				GetStdDev() const;
		std::vector<double>					Sample() const;
		double								ApproxCDF(const std::vector<double> &x, unsigned nsamples = 10000) const;
		double								GetLnPDF(const std::vector<double> &x) const;
		double								GetRelativeLnPDF(const std::vector<double> &x) const;
		void 								SetMeanAndVariance(const std::vector<double> &m, const std::vector<double> &v);
		unsigned							GetNParams() const;
#		if defined(PYTHON_ONLY)
		//void								AltSetMeanAndVariance(std::vector<double> m, std::vector<double> v);
		double_vect_t						GetVarCovarMatrix();
#		endif

    private:
	
		unsigned							dim;			/**< The dimension, which equals the number of parameters and is used to transform the mean, variance and density from Dirichlet to values appropriate for a relative rate distribution */
		double								sum_params;		/**< The sum of the dirichlet parameters stored in the data member `dirParams', which is provided by the base class */
	};
	
} // namespace phycas

#endif
