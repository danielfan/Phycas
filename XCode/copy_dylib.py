#!/usr/bin/python

import shutil,os
product=os.environ['PRODUCT_NAME']
phycas_root= os.environ['PHYCAS_ROOT']
build_dir= os.environ['TARGET_BUILD_DIR']
fromfile=os.path.join(build_dir,product)+'.dylib'
# drop the leading "_" and  trailing "Ext", and rename *.dylib to *.so
tofile=os.path.join(phycas_root,'phycas',product[1:-3], product+'.so')
print "copying fromfile=",fromfile," to tofile=",tofile
shutil.copyfile(fromfile, tofile)
