#if ! defined(STATES_PATTERNS_HPP)
#define STATES_PATTERNS_HPP

#include "pyphy/src/cipres/ConfigDependentHeaders.h"	// for int8_t typedef

typedef std::vector<int8_t>				VecStateList;
typedef std::vector<unsigned>			VecStateListPos;

typedef float										PatternCountType;
typedef	std::map<VecStateList, PatternCountType>	PatternMapType;
typedef	std::vector<PatternCountType>				CountVectorType;

#endif
