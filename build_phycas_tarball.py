#/usr/bin/env python

# This script creates a tar archive (tarfile) inside the tardir directory. The tarfile
# contains the rootdir directory and everything inside it. Normally, rootdir is defined
# to be the phycas directory. This makes it possible to build phycas in one's home directory
# and then easily copy the entire phycas package to the site-packages directory of the 
# Python installation. The src subdirectory and all .svn subdirectories, which are not
# needed for a functioning phycas package, are excluded.  The list of file name extensions 
# that are matched is given in the variable validset. When the program finishes, it spits
# out a list of file extensions that were ignored because they were not in validset.
#
# Note: after you unpack the tarfile, you will probably need to run "ldconfig" or specify
# a "LD_LIBRARY_PATH" in order for the copied .so files to be loadable. For example, if
# you unpack tarfile at "/usr/lib/python2.5/site-packages", here is what you need to do:
#
# /sbin/ldconfig /usr/lib/python2.5/site-packages/phycas/Conversions
# 
# or, alternatively, define the LD_LIBRARY_PATH environmental variable like this (in 
# your .tcshrc file):
#
# setenv LD_LIBRARY_PATH /usr/lib/python2.5/site-packages/phycas/Conversions
#
# or like this (in your .bashrc file):
# 
# export LD_LIBRARY_PATH="/usr/lib/python2.5/site-packages/phycas/Conversions"

import os
import sets

# Specify the name of the tar file to be created, and the directory in which it should be created
tardir = '.'
tarfile = 'phycas_distr.tar'

# Specify the root directory (everything in this dir will be tarred except those files and dirs specifically excluded)
rootdir = './phycas'

# Spedify directories to be ignored
rmlist = ['.svn','src']

# Spedify valid file name extensions (files with other extensions will be skipped)
validset = sets.Set(['.py','.so','.afm','.p','.t','.tre','.nex','.pdf','.html','.txt'])

# This set will hold file name extensions encountered that were not in validset
invalid = sets.Set()

# This function handles each directory encounteed in os.path.walk
def visit(outf, dirname, names):
        # walk lets you delete directories from names, which prevents these from being visited
        for r in rmlist:
                if r in names:
                        names.remove(r)
        for fn in names:
                fp = os.path.join(dirname,fn)
                if not os.path.isdir(fp):
                        x = os.path.splitext(fp)[1]
                        if x in validset:
                                os.system('tar rf %s %s' % (tarpath, fp))
                                outf.write('%s\n' % fp)
                        else:
				a,b = os.path.splitext(fn)
				if a.find('libboost_python') > -1:
					os.system('tar rf %s %s' % (tarpath, fp))
					outf.write('%s\n' % fp)
				else:
                                	invalid.add(x)

# Open the manifest file, which will contain the path of every file included in the tar file
outf = open('distr_manifest.txt', 'w')

# If the tarball already exists, delete it
tarpath = os.path.join(tardir, tarfile)
if os.path.exists(tarpath):
        os.remove(tarpath)

# Conduct the directory walk
os.path.walk(rootdir, visit, outf)

# Finish up by closing the manifest file and printing a list of the invalid extensions encountered
outf.close()
print 'These extensions ignored:'
for x in invalid:
        print x
