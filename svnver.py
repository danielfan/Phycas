#!/bin/env python

import subprocess, sys
try:
	svn_revision = subprocess.Popen('svnversion', shell=False, stdout=subprocess.PIPE).communicate()[0].strip()
except OSError:
	svn_revision = '?'
open('phycas/svnver.txt', 'w').write('%s\n' % (svn_revision,))
