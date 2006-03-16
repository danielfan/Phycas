/*	linalg.h
|
|	Prototypes for matrix-inversion and eigensystem functions
|
|	Copyright (c) 2002 by David L. Swofford, Florida State University.
|	All rights reserved.
*/

#ifdef __cplusplus
extern "C" {
#endif

#if defined(PAUP)
//@ below defined TEMP
#	define NEED_NONSYM_EIG		/* leave undefined for PAUP* -- doesn't need nonsymmetric eigenstuff */
#	define NEED_MATINV			/* leave undefined for PAUP* -- doesn't need matrix inversion */
#endif

#if defined(NEED_NONSYM_EIG)
#	define RC_COMPLEX_EVAL 2	/* return code that complex eigenvalue obtained */
	extern int	EigenRealGeneral(int, double **, double *, double *, double **, int *, double *);
#endif

#if defined(NEED_MATINV)
	extern int	InvertMatrix(double **, int, double *, int *, double **);
#endif

extern int	LUDecompose(double **, int, double *, int *, double *);
extern int	EigenRealSymmetric(int, double **, double *, double **, double *);

#ifdef __cplusplus
}
#endif
