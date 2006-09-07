#ifndef NCL_NXS_COMMAND_COMMON_H
#define NCL_NXS_COMMAND_COMMON_H

#include "phypy/src/ncl/nxs_defs.hpp"

class NxsCmdOption;
class NxsCommand;
class NxsTest;
typedef boost::shared_ptr<NxsTest>		NxsTestShPtr;//NOTE:this typedef is also in nxstest.h change it in both places
typedef boost::shared_ptr<NxsCommand>	NxsCommandShPtr;
typedef std::vector<NxsCommandShPtr> 	VecNxsCommandShPtr;	
typedef boost::shared_ptr<NxsCmdOption> NxsCmdOptionShPtr;
typedef std::vector<NxsCmdOptionShPtr> 	VecNxsCmdOptionShPtr;

typedef std::vector<const NxsCmdOption *> 	VecConstNxsCmdOptionPtr;


#endif
