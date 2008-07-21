import sys
from phycas import OutputFilter, getDefaultOutFilter

output_stream = None
outputter = None

def getDefaultOutputContainer():
    return getDefaultOutputter()

def writeMessageToStdOut(m):
    "Writes a message to standard out and adds a trailing newline"
    sys.stdout.write("%s\n" % m)

def getDefaultOutputStream():
    return writeMessageToStdOut

def getDefaultOutputter():
    global outputter
    if outputter is None:
        outputter = OutputFilter(getDefaultOutFilter(), getDefaultOutputStream())
    return outputter

class CommonFunctions(object):
    def __init__(self, opts):
        self.log_file_spec = None
        self.quiet = False
        self.opts = opts
        try:
            o = self.opts.out
            self.optsout = o
            self.stdout = o.getStdOutputter()
            self.open()
        except:
            self.optsout = None
            self.stdout = getDefaultOutputter()

    def __del__(self):
        self.close()

    def open(self):
        if (self.log_file_spec is None):
            try:
                self.log_file_spec = self.optsout.log
            except:
                pass
            if self.log_file_spec:
                self.log_file_spec.openAsLog(self.stdout)
    def close(self):
        if self.log_file_spec:
            self.log_file_spec.close()
            self.log_file_spec = None
        
    def phycassert(self, assumption, msg):
        self.stdout.phycassert(assumption, msg)

    def output(self, msg=""):
        self.stdout.info(msg)
    info = output
            
    def abort(self, msg):
        self.stdout.abort(m)
    def warning(self, msg):
        self.stdout.warning(m)
    def error(self, msg):
        self.stdout.error(m)
    def add_mirror(self, m):
        self.stdout.add_mirror(m)
    def remove_mirror(self, m):
        self.stdout.remove_mirror(m)



