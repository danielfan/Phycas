/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
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

#ifndef PHO_CHARACTERS_MANAGER_H
#define PHO_CHARACTERS_MANAGER_H

#include <list>
#include "phycas/src/ncl/characters/nxs_characters_manager.hpp"
#include "phycas/src/ncl/characters/data_pattern.hpp"
class NxsDataType;
class PhoCharactersManager;
class ConstSiteInfo;
class OrdCodedArrs;

/*----------------------------------------------------------------------------------------------------------------------
|	Stores Data Patterns (columns in the data matrix) along with their weights (sum of counts*weights)
*/
class CompressedMatrix
	{
	public:
		~CompressedMatrix();
		
		void				FillInTaxonsCondLikeArray(double *clArray, unsigned taxInd, const unsigned timesToCopy) const;
		void				FillInConstantSiteInfo(ConstSiteInfo *) const;
		const NxsDataType & GetDataType() const 
			{
			return dataType;
			}
		unsigned 			GetNumTaxa() const 
			{
			return nTaxa;
			}
		unsigned 			GetNumPatterns() const 
			{
			return (unsigned)patterns.size();
			}
		std::vector<double> GetPatternWeights() const;
		const DataPatternInfo &GetPatternInfo(unsigned i) const;
		void				GetTaxonsOrdCodedStates(OrdCodedArrs *oca, unsigned rowIndex) const;
		
	protected:
		CompressedMatrix(const PhoCharactersManager *, const NxsDataType &desiredType, const NxsIndexSet &taxa, const NxsIndexSet &characters);
		
		typedef std::pair<DataPattern, DataPatternInfo>	PatternAndInfo;
		typedef std::vector<PatternAndInfo>				VecPattern;
		typedef VecPattern::const_iterator				PatConstIter;
		typedef VecPattern::iterator					PatIter;
		
		NxsDataType			dataType;
		unsigned 			nTaxa;
		VecPattern			patterns;		/**< sorted collection of the pattern (which is owned by the CompressedMatrix) and info about the pattern */
		
		friend class PhoCharactersManager;
	};
typedef boost::shared_ptr<const CompressedMatrix> CompressedMatrixShPtr;
typedef std::pair<CompressedMatrixShPtr, NxsIndexSet> MatrixAndTaxa;
  			
/*----------------------------------------------------------------------------------------------------------------------
|	adds interface for obtaining compressed data matrices to the NxsCharactersManager.
*/
class PhoCharactersManager : public NxsCharactersManager
	{
	public:
		CompressedMatrixShPtr 	GetCompressedMatrix(const NxsDataType &desiredType, const NxsIndexSet &taxaIndices);
  		PhoCharactersManager(PhoTaxaManager & taxaMgr);
  	protected:
  		void  					CreateAndInsertDataPatterns(SortedPatterns *, const NxsDataType &desiredType, const NxsIndexSet &taxaIndices, const NxsIndexSet &charIndices) const;
  		
  		typedef std::vector<MatrixAndTaxa> CompressedMatrices;
  		MatrixAndTaxa activeMatrixCompressed; /* shared pointers to matrices that have been compressed */
  		
	private:
		friend class CompressedMatrix;
	};

#endif

