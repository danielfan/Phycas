import subprocess, sys, os

if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
    phycas_version = open(sys.argv[1],'r').read().strip()
else:
    phycas_version = '?'
    
try:
    svn_revision = subprocess.Popen('svnversion', shell=False, stdout=subprocess.PIPE).communicate()[0].strip()
except OSError:
    svn_revision = '?'
    
open('phycas/svnver.txt', 'w').write('%s-r%s\n' % (phycas_version,svn_revision))
