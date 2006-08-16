#ifndef CIPRES_NATIVE_C_H
#define CIPRES_NATIVE_C_H

#include "pyphy/src/cipres/ConfigDependentHeaders.h"

#ifdef __cplusplus
extern "C" 
{
#endif

typedef enum { CIPR_noScore, CIPR_intScore, CIPR_doubleScore } CIPR_ScoreType;
typedef struct
	{
	double			doubleScore;
	long			intScore;
	CIPR_ScoreType	type;
	} CIPR_TreeScore;
	
#if defined(OLD_CIPRES_MATRIX)

	typedef char	CIPR_MatrixElement;
	typedef			CIPR_MatrixElement **CIPR_Matrix;
	extern void showMatrix(const CIPR_Matrix matrix, int m, int n);
	
#else
	
	typedef int8_t CIPR_State_t; /** type used to enumerate possible states.	
									-1 is used for gaps, other negative flags may be added later.
									This size limits the maximum number of states allowed. */
	typedef int8_t CIPR_StateSet_t; /** type used to refer to unique combinations of states (the "fundamental" states and ambiguity codes)
									-1 is used for gaps.  To handle all possible data sets, this must be large enough to hold
									2^(nStates + 1) values if the datatype allows gaps.  Thus using int8_t limits us to 8 states */
	
	typedef enum {CIPR_DNA_Datatype, CIPR_RNA_Datatype, CIPR_AA_Datatype, CIPR_Codon_Datatype, CIPR_Generic_Datatype} CIPR_Datatypes;
	typedef struct CIPRES_Matrix
		{
		CIPR_State_t 	  *stateList; 		/** Flattened array of array of observed states.  If more than one state was observed, then the first element is the number of states observed.  
												Exceptions: -1 is for gaps, nStates is for missing. */
		unsigned 		   *stateListPos;  	/** Maps a state set code (the elements of the matrix) to the index in ambigList where the states are listed */
		CIPR_StateSet_t	  **matrix;			/** taxa x characters matrix of indices of state sets */
		const char		   *symbolsList;	/** array of the characters used to stand for each state ("ACGT?NRY" for example) //@temp paup depends on all symbols being unique (even ambiguity codes)*/
		unsigned			nStates;
		unsigned			nChar;
		unsigned			nTax;
		unsigned			nObservedStateSets; /* the length of stateListPos */ 
		CIPR_Datatypes		datatype;
		} CIPR_Matrix;

	extern void showMatrix(const CIPR_Matrix matrix);
	extern const char * datatypeEnumToString(const CIPR_Datatypes d);
#endif

extern int CIPR_AllocMatrix(void *pA, size_t elSize, unsigned nrows, unsigned ncols);

#define EXC_MSG_SIZE  256
struct AppException
{
	char errorMessage[EXC_MSG_SIZE];
	int errorCode;
};

#ifdef __cplusplus
}	
#endif

#endif /* __TREEINFER_HELPER_H */
