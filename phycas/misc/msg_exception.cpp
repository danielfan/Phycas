#include "phycas/force_include.h"
#include "phycas/misc/msg_exception.hpp"
#include "ncl/misc/string_extensions.hpp"
using std::string;
MsgException::MsgException(
  const char * s, 
  const char *file, int line) 
  { 
  msg << s << " (file " << file << ", line" << line << ')'; 
  }

MsgException::MsgException(
  const string &s, 
  const char *file, 
  int line) 
  {
  msg << s << " (file " << file << ", line" << line << ')'; 
  }
