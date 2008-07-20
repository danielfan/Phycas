import sys

class CommonFunctions(object):
    def __init__(self):
        self.logf = None
        self.quiet = False
        
    def phycassert(self, assumption, msg):
        if not assumption:
            if phycassertRaisesException:
                raise AssertionError(msg)
            sys.exit('Error: ' + msg)

    def output(self, msg = None):
        if not self.quiet:
            if msg:
                print msg
            else:
                print
        if self.logf is not None:
            if not msg:
                msg = ''
            self.logf.write('%s\n' % (msg))
            self.logf.flush()
            
    def abort(self, msg):
        s = '\n***** Fatal error: %s' % msg
        sys.exit(s)
        
    def warning(self, msg):
        s = '\n***** Warning: %s' % msg
        self.output(s)
        
