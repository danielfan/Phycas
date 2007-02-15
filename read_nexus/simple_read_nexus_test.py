from phycas.ReadNexus import *
from phycas.DataMatrix import *
import sys, os

defaultFileName = 'nyldna4.nex'

def askokcancel(title,msg, automaticMode= True):
    print title
    if automaticMode:
        print msg
        return True
    else:
        resp = raw_input(msg)
        return bool(resp)
	    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        fn = defaultFileName
    else:
        fn = sys.argv[1]
    r = NexusReader(4)
    try:
        r.ReadFile(fn)
    except ValueError, e:
        askokcancel('Error:', '%s' % e)
    except int:
    	sys.exit('Unknown exception!')
    else:
        askokcancel('Success!', 'Number of characters read = %d' % r.GetNChar())
        discMatrix = getDiscreteMatrix(r, 0)
        print 'According to the discrete matrix, nchar =', discMatrix.getNChar()
