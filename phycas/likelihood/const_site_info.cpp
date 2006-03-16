#include "phycas/force_include.h"
#include "phycas/likelihood/const_site_info.hpp"
#include "ncl/misc/nxs_data_type.hpp"
using std::vector;

OrdCodedArrs::OrdCodedArrs(
  const vector<DataStorageType *> &bitCode,  
  unsigned wordLen,
  unsigned numStates)
  //POL-121803 : Compressed2DArray<NStateInt>(bitCode.size()),
  //POL-121803 nStates(numStates)
  : Compressed2DArray<NStateInt>((unsigned)bitCode.size(), numStates)
  	{
	typedef vector<DataStorageType *>::const_iterator BitCodePtrIter;
  	for (BitCodePtrIter bIt = bitCode.begin(); bIt != bitCode.end(); ++bIt)
  		{
  		vector<unsigned> ordCode = ConvertBitArrToVecUInt(*bIt, wordLen);
  		Append(ordCode);
  		}
  	}

