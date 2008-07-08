# This script creates a tar file containing the phycas directory and everything inside it
# needed for running Phycas. This makes it possible to build phycas in one's home directory
# and then easily copy the relevant parts of the phycas directory to the site-packages
# directory of the Python installation on a linux machine, for example. 
# The src subdirectory and all .svn subdirectories are completely excluded, and the list
# of file name extensions that are matched is given in the variable validset. The only
# thing this script does not do is move the libboost*.so file and its symbolic links,
# which you will have to do manually.

import os
import sets

# Specify the name of the tar file to be created
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
                                invalid.add(x)

# Open the manifest file, which will contain the path of every file included in the tar file
outf = open('distr_manifest.txt', 'w')

# If the tarball already exists, delete it
tarpath = os.path.join(rootdir, tarfile)
if os.path.exists(tarpath):
        os.remove(tarpath)

# Conduct the directory walk
os.path.walk(rootdir, visit, outf)

# Finish up by closing the manifest file and printing a list of the invalid extensions encountered
outf.close()
print 'These extensions ignored:'
for x in invalid:
        print x
