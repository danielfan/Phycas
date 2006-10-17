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

#if defined(USING_NUMARRAY)
#	define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#	define NO_IMPORT_ARRAY
#endif

#include "probability_distribution.hpp"
#include "ncl/nxs_exception.hpp"

namespace phycas
{

DirichletDistribution::DirichletDistribution(const std::vector<double> &params) 
	{
	unsigned params_size = (unsigned)params.size();
	//std::cerr << "params_size = " << params_size << std::endl;
	PHYCAS_ASSERT(params_size > 1);

	// Reserve necessary space to save copying
	//
	dirParams.reserve(params_size);
	scratchSpace.reserve(params_size);
	paramDistributions.reserve(params_size);

#	if defined NCL_NXS_THROW_UNDEFINED
		if (params.size() < 2)
			throw NxsX_UndefinedException("Illegal Dirichlet", __FILE__, __LINE__);
#	endif

	std::vector<double>::const_iterator iter;
	for (iter = params.begin(); iter != params.end(); ++iter)
		{
		double c_i = *iter;
		dirParams.push_back(c_i);
		scratchSpace.push_back(0.0);
		GammaDistribution gamma_i = GammaDistribution(c_i, 1.0);
		gamma_i.SetLot(&myLot);	// ensure that they all use same random number generator
		paramDistributions.push_back(gamma_i);
		}
	}

} // namespace phycas
