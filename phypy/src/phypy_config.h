#define PYTHON_ONLY	// for code that only makes sense when exported to Python (most code can be used to construct a pure C++ program)
#define POL_PHYCAS	// for POL-specific code related to PHYCAS (this one needs to go)
//#define POL_PYPHY	// for POL-specific code related to boost python (this one needs to go too)
#define NO_IDL_TYPES
#define POLPY_NEWWAY  1
#define POLPY_OLDWAY  !(POLPY_NEWWAY)
#include "phypy/src/phycas_config.h"
