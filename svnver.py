import subprocess, sys, os

if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
    phycas_version = open(sys.argv[1],'r').read().strip()
else:
    phycas_version = '?'
    
try:
    p = subprocess.Popen(['git', 'show'], shell=False, stdout=subprocess.PIPE)
    t = p.communicate()
    svn_revision = t[0].split()[1]
except OSError:
    svn_revision = '?'
    
open('phycas/svnver.txt', 'w').write('%s git SHA1 %s\n' % (phycas_version,svn_revision))
