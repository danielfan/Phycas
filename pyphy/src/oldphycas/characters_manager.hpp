#ifndef PHO_CHARACTERS_MANAGER_H
#define PHO_CHARACTERS_MANAGER_H

#include <list>
#include "pyphy/src/ncl/characters/nxs_characters_manager.hpp"
#include "pyphy/src/ncl/characters/data_pattern.hpp"
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

