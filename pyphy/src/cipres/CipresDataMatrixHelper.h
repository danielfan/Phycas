#if !defined(CIPRES_DATA_MATRIX_HELPERS_H)
#define CIPRES_DATA_MATRIX_HELPERS_H

#include <cassert>
#include <string>
#include <vector>
#include "pyphy/src/cipres/CipresNativeC.h"
#include "pyphy/src/cipres/AllocateMatrix.hpp"

#ifdef __cplusplus
#if !defined(NO_IDL_TYPES)
#	include "pyphy/src/cipres/CipresHelper.h"
#endif
	namespace CipresNative 
	{
#	if !defined(NO_IDL_TYPES)
		CIPR_Datatypes corbaToNativeDatatype(const CipresIDL_api1::DiscreteDatatypes dt);
		CipresIDL_api1::DiscreteDatatypes nativeToCorbaDatatype(const CIPR_Datatypes dt);
#	endif

		/**
		 * Matrix representation that serves as a bridge between CORBA communication (CipresIDL_api1::DataMatrix) and 
		 * use of the matrix in a C application using (CIPR_Matrix).
		 */
	class DiscreteMatrix 
		{
		public:
#			if !defined(NO_IDL_TYPES)
				DiscreteMatrix(const CipresIDL_api1::DataMatrix &);
#			endif
			DiscreteMatrix(const CIPR_Matrix & );
			
			const CIPR_Matrix & getConstNativeC() const
				{
				return nativeCMatrix;
				}
			CIPR_Matrix & getNativeC()
				{
				return nativeCMatrix;
				}
			unsigned	getNChar() const
				{
				return nativeCMatrix.nChar;
				}
			unsigned	getNTax() const
				{
				return nativeCMatrix.nTax;
				}
			unsigned	getNStates() const
				{
				return nativeCMatrix.nStates;
				}
			const char *	getSymbolsList() const   //POL added 15-Nov-2005 
				{
				return nativeCMatrix.symbolsList;
				}
			const std::vector<int8_t> &getStateList() const
				{
				return stateListAlias;
				}
			const std::vector<unsigned> &getStateListPos() const
				{
				return stateListPosAlias;
				}
			const CIPR_StateSet_t *getRow(unsigned i) const 
				{
				assert(i < nativeCMatrix.nTax);
				return nativeCMatrix.matrix[i];
				}
			const CIPR_StateSet_t * const * getMatrix() const 
				{
				return nativeCMatrix.matrix;
				}
#			if !defined(NO_IDL_TYPES)
				::CipresIDL_api1::DataMatrix * _retn() const;
#			endif
		private:
			typedef ScopedTwoDMatrix<CIPR_StateSet_t> ScopedStateSetTwoDMatrix;
			
			CIPR_Matrix					nativeCMatrix; 		/** taxa x characters matrix in a C struct*/
			std::string					symbolsStringAlias;	/** memory management alias to symbols field of nativeCMatrix */
			ScopedStateSetTwoDMatrix	matrixAlias;		/** memory management alias to matrix field of nativeCMatrix */
			std::vector<CIPR_State_t>	stateListAlias;		/** memory management alias to ambigList field of nativeCMatrix */
			std::vector<unsigned>		stateListPosAlias;		/** memory management alias to symbolsList field of nativeCMatrix */
			
			DiscreteMatrix(const DiscreteMatrix &); /** don't define, not copyable*/
			DiscreteMatrix & operator=(const DiscreteMatrix &); /** don't define, not copyable*/
		};
	} // CipresNative
#endif //__cplusplus

#endif  // CIPRES_DATA_MATRIX_HELPERS_H
