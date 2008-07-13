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

#if !defined(CIPRES_DATA_MATRIX_HELPERS_H)
#define CIPRES_DATA_MATRIX_HELPERS_H

#include <string>
#include <vector>
#include "phycas/src/cipres/CipresNativeC.h"
#include "phycas/src/cipres/AllocateMatrix.hpp"
#include "ncl/nxscharactersblock.h"

#ifdef __cplusplus

	namespace CipresNative 
	{

		/**
		 * Matrix representation that serves as a bridge between CORBA communication (CipresIDL_api1::DataMatrix) and 
		 * use of the matrix in a C application using (CIPR_Matrix).
		 */
	class DiscreteMatrix 
		{
		public:
			DiscreteMatrix(const CIPR_Matrix & );
			DiscreteMatrix(const NxsCharactersBlock & cb, bool convertGapsToMissing);

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
				PHYCAS_ASSERT(i < nativeCMatrix.nTax);
				return nativeCMatrix.matrix[i];
				}
			const std::vector<int8_t> getRowAsVector(unsigned i) const 
				{
				PHYCAS_ASSERT(i < nativeCMatrix.nTax);
				std::vector<int8_t> v;
				for (unsigned j = 0; j < nativeCMatrix.nChar; j++)
					{
					v.push_back(nativeCMatrix.matrix[i][j]);
					}
				return v;
				}
			const CIPR_StateSet_t * const * getMatrix() const 
				{
				return nativeCMatrix.matrix;
				}
			const int getDatatype() const
				{
				return (int)nativeCMatrix.datatype;
				}
			bool hasWeights() const
				{
				return hasIntWeights() || hasDblWeights();
				}

			bool hasIntWeights() const
				{
				return !(intWts.empty());
				}
			
			bool hasDblWeights() const
				{
				return !(dblWts.empty());
				}
			std::vector<int> & getIntWeights()
				{
				return intWts;
				}
			std::vector<double> & getDblWeights()
				{
				return dblWts;
				}
			const std::vector<int> & getIntWeightsConst() const
				{
				return intWts;
				}
			const std::vector<double> & getDblWeightsConst() const
				{
				return dblWts;
				}
		private:
			typedef ScopedTwoDMatrix<CIPR_StateSet_t> ScopedStateSetTwoDMatrix;
			
			CIPR_Matrix					nativeCMatrix; 		/** taxa x characters matrix in a C struct*/
			std::string					symbolsStringAlias;	/** memory management alias to symbols field of nativeCMatrix */
			ScopedStateSetTwoDMatrix	matrixAlias;		/** memory management alias to matrix field of nativeCMatrix */
			std::vector<CIPR_State_t>	stateListAlias;		/** memory management alias to ambigList field of nativeCMatrix */
			std::vector<unsigned>		stateListPosAlias;		/** memory management alias to symbolsList field of nativeCMatrix */
			std::vector<int>			intWts;
			std::vector<double>			dblWts;
			DiscreteMatrix(const DiscreteMatrix &); /** don't define, not copyable*/
			DiscreteMatrix & operator=(const DiscreteMatrix &); /** don't define, not copyable*/
		};
	} // CipresNative
#endif //__cplusplus

#endif  // CIPRES_DATA_MATRIX_HELPERS_H
