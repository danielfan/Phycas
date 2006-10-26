"""fakecompiler

Copied from the mwerkscompiler.
This "compiler" just copies prebuilt files so that extensions with complex build 
    systems (e.g. those that use the Boost python library and bjam) can fake
    compilation so that they fit into the setup.py system.
Written so that we could use bdist_mpkg for the Mac installer for phycas.
"""

#This is a hack.
#I had to add:
#                 'fake':     ('fakecompiler', 'FakeCompiler',
#                               "Dummy for a compiler that clones already built sources"),
#   to distutils.ccompiler line 1112

__revision__ = "$Id: /mirror/python24-fat/Lib/distutils/fakecompiler.py 3 2006-02-12T02:26:15.597018Z mholder@scs.fsu.edu  $"

import sys, os, string
#from types import *
from distutils.ccompiler import CCompiler
import distutils.util
import distutils.dir_util
from distutils import log

class FakeCompiler (CCompiler) :
    """Concrete class that implements an interface to MetroWerks CodeWarrior,
       as defined by the CCompiler abstract class."""

    compiler_type = 'fake'

    # Just set this so CCompiler's constructor doesn't barf.  We currently
    # don't use the 'set_executables()' bureaucracy provided by CCompiler,
    # as it really isn't necessary for this sort of single-compiler class.
    # Would be nice to have a consistent interface with UnixCCompiler,
    # though, so it's worth thinking about.
    executables = {}

    # Private class data (need to distinguish C from C++ source for compiler)
    _c_extensions = ['.c']
    _cpp_extensions = ['.cc', '.cpp', '.cxx']
    _rc_extensions = ['.r']
    _exp_extension = '.exp'

    # Needed for the filename generation methods provided by the
    # base class, CCompiler.
    src_extensions = (_c_extensions + _cpp_extensions +
                      _rc_extensions)
    res_extension = '.rsrc'
    obj_extension = '.obj' # Not used, really
    static_lib_extension = '.lib'
    shared_lib_extension = '.slb'
    static_lib_format = shared_lib_format = '%s%s'
    exe_extension = ''
    def __init__ (self,
                  verbose=0,
                  dry_run=0,
                  force=0):
        CCompiler.__init__ (self, verbose, dry_run, force)
        self.verbose = verbose or 1


    def compile (self,
                 sources,
                 output_dir=None,
                 macros=None,
                 include_dirs=None,
                 debug=0,
                 extra_preargs=None,
                 extra_postargs=None,
                 depends=None):
        print output_dir
        return [] #don't need to compile


    def link (self,
              target_desc,
              objects,
              output_filename,
              output_dir=None,
              libraries=None,
              library_dirs=None,
              runtime_library_dirs=None,
              export_symbols=None,
              debug=0,
              extra_preargs=None,
              extra_postargs=None,
              build_temp=None,
              target_lang=None):
        print 'LINK:'
        destParent, file = os.path.split(output_filename)
        if not os.path.exists(destParent):
            os.makedirs(destParent)
        exPostDict = {}
        for i in extra_postargs:
            s = i.split('=')
            if len(s) == 1:
                exPostDict[s[0]] = None,
            else:
                exPostDict[s[0]] = '='.join(s[1:])
        if '--skip' in exPostDict:
            return
        if not '--built-under' in exPostDict:
            raise ValueError, '"--built-under=<path>" must be used the "extra_link_args" argument to the Extensions __init__ to use the fakecompiler'
        assert '--path-from-package' in exPostDict
        if not '--path-from-package' in exPostDict:
            raise ValueError, '"--path-from-package=<path>" must be used the "extra_link_args" argument to the Extensions __init__ to use the fakecompiler'
        sourceParent = os.path.join(exPostDict['--built-under'], exPostDict['--path-from-package'])
        if not os.path.exists(sourceParent):
            raise ValueError, '%s does not exist' % sourceParent
        source = os.path.join(sourceParent, file)
        self._copyFile(source, output_filename)
        if '--colocate-lib' in exPostDict:
            source = os.path.join(sourceParent, exPostDict['--colocate-lib'])
            dest = os.path.join(destParent, exPostDict['--colocate-lib'])
            self._copyFile(source, dest)

    def _copyFile(self, source, dest):
        import shutil
        if self.verbose or 1:
            print 'Copying %s to %s' %(source, dest)
        shutil.copy(source, dest)
        
    def library_dir_option (self, dir):
        """Return the compiler option to add 'dir' to the list of
        directories searched for libraries.
        """
        return

    def runtime_library_dir_option (self, dir):
        """Return the compiler option to add 'dir' to the list of
        directories searched for runtime libraries.
        """
        return

    def library_option (self, lib):
        """Return the compiler option to add 'dir' to the list of libraries
        linked into the shared library or executable.
        """
        return

    def find_library_file (self, dirs, lib, debug=0):
        """Search the specified list of directories for a static or shared
        library file 'lib' and return the full path to that file.  If
        'debug' true, look for a debugging version (if that makes sense on
        the current platform).  Return None if 'lib' wasn't found in any of
        the specified directories.
        """
        return 0
