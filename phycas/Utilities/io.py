import os, sys, subprocess

_phycas_dir = None

def getPhycasDir():
    "Returns the absolute path to the directory that is the top of the phycas package"
    global _phycas_dir
    if _phycas_dir is None:
        _phycas_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    return _phycas_dir
    
def getPhycasTestDir():
    return os.path.join(getPhycasDir(), "Tests")

def getPhycasTestDataDir():
    return os.path.join(getPhycasTestDir(), "Data")

def getPhycasTestData(filen):
    """Takes a string `filen` that represents the name of one of the files in
    Phycas' test suite.
    The function returns the full path to the test file."""
    return os.path.join(getPhycasTestDataDir(), filen)

def _runRegressionTests(out):
    d = getPhycasTestDir()
    r = os.path.join(d, "runall.py")
    spawnPython(r)
    r = os.path.join(d, "doctestall.py")
    spawnPython(r)
    out.info("\nAll tests passed.")
    
        
def spawnPython(f):
    if not os.path.exists(f):
        raise RuntimeError('The python script "%s" does not exist' % f)
    retcode = subprocess.call([sys.executable, f])
    if retcode < 0:
        raise RuntimeError("python execution of %s was terminated by signal %s" % (f, str(-retcode)))
    elif retcode > 0:
        raise RuntimeError("python execution of %s was failed with code %s" % (f, str(retcode)))
