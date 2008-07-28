import os, sys, subprocess
from phycas.ReadNexus._NexusReader import FileFormats
from phycas.Utilities.CommonFunctions import getDefaultOutputter


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

def _readFileSanityCheck(filepath, format=FileFormats.NEXUS, out=None):
    if not format == FileFormats.NEXUS:
        if out is None:
            out = getDefaultOutputter()
        out.phycassert(format == FileFormats.NEXUS, "Currently only the NEXUS format is supported")
    if not os.path.exists(filepath):
        raise ValueError('The file "%s" does not exist' % filepath)
    
def readData(filepath, format=FileFormats.NEXUS, out=None):
    """Returns a data matrix (or None if there is no data) from `filepath`
    
    Currently only supports NEXUS and only returns the last data matrix, but
    this will be generalized to read other formats and return the 
    supermatrix of all data matrices in the file."""
    _readFileSanityCheck(filepath, format, out)
    from phycas.ReadNexus import NexusReader
    reader = NexusReader()
    reader.readFile(filepath)
    return reader.getLastDiscreteMatrix(True)



def readTrees(filepath, format=FileFormats.NEXUS, out=None):
    """Returns a (taxon_list, list of TreeDescription objects) (or an empty list if there are no trees) from `filepath`
    
    Currently only supports NEXUS and only returns the last data matrix, but
    this will be generalized to read other formats and return the 
    supermatrix of all data matrices in the file."""

    _readFileSanityCheck(filepath, format, out)
    from phycas.ReadNexus import NexusReader
    reader = NexusReader()
    reader.readFile(filepath)
    return reader.taxa, reader.getTrees()

def readFileToSourceObjects(filepath, format=FileFormats.NEXUS, out=None):
    """Returns a (DataSource, TreeCollection) from the file `filepath`
    
    Currently only supports NEXUS and only returns the last data matrix, but
    this will be generalized to read other formats and return the 
    supermatrix of all data matrices in the file."""

    _readFileSanityCheck(filepath, format, out)
    from phycas.ReadNexus import NexusReader
    reader = NexusReader()
    reader.readFile(filepath)
    t = reader.taxa
    dm = DataSource(reader.getLastDiscreteMatrix(True), taxon_labels=t)
    tm = TreeCollection(reader.getTrees(), taxon_labels=t)
    return dm, tm

class TreeCollection(object):
    def __init__(self, arg, **kwargs):
        """`arg` can be a string or an iterable containing trees.
        If `arg` is a string, then the `format` keyword argument can be used to
        specify the file format (the default is NEXUS).
        If the trees are passed in as a list then a `taxon_labels` keyword
        argument is expected.
        """
        self.title = kwargs.get("title")
        if isinstance(arg, str):
            self.filename = arg
            self.format = kwargs.get("format", FileFormats.NEXUS)
            self.trees = None
            self.taxon_labels = None
        else:
            if arg is None:
                self.trees = None
            else:
                self.trees = list(arg)
            self.taxon_labels = kwargs.get("taxon_labels")
            self.filename = None
            self.format = None
            
    def __iter__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns an iterable over the trees (this will trigger the reading of
        the input file if the object was initialized from a string.
        Stores tree descriptions in self.stored_tree_defs and taxon labels in
        self.taxon_labels.
        """
        if self.trees is None:
            if not self.filename:
                return iter([])
            self.taxon_labels, self.trees = readTrees(self.filename, self.format)
        return iter(self.trees)

    def __str__(self):
        if self.title:
            return self.title
        s, t = "", ""
        d = self.trees
        if d:
            s = "Collection of %d trees in memory (id %d). " % (len(self.trees), id(self.trees))
        if self.filename:
            t =  "Trees from the file %s" % repr(self.filename)
        if s or t:
            return s + t
        return "None"

    def __deepcopy__(self, memo):
        #trees are expensive, so we don't make a deepcopy
        if self.trees:
            return TreeCollection(self.trees, taxon_labels=self.taxon_labels, title=self.title)
        elif self.filename:
            return TreeCollection(self.filename, format=self.format, title=self.title)
        return TreeCollection(None, taxon_labels=self.taxon_labels, title=self.title)

    def writeTree(self, tree, name="", rooted=None):
        if not self.trees:
            self.trees = []
        self.trees.append(tree)

    def finish(self):
        pass

class DataSource(object):
    def __init__(self, arg, **kwargs):
        """`arg` can be a string or an iterable containing trees.
        If `arg` is a string, then the `format` keyword argument can be used to
        specify the file format (the default is NEXUS).
        If the trees are passed in as a list then a `taxon_labels` keyword
        argument is expected.
        """
        self.title = kwargs.get("title")
        if isinstance(arg, str):
            self.filename = arg
            self.format = kwargs.get("format", FileFormats.NEXUS)
            self.data_obj = None
            self.taxon_labels = None
        else:
            self.data_obj = arg
            self.taxon_labels = kwargs.get("taxon_labels")
            self.filename = None
            self.format = None
            

    def getMatrix(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns an iterable over the trees (this will trigger the reading of
        the input file if the object was initialized from a string.
        Stores tree descriptions in self.stored_tree_defs and taxon labels in
        self.taxon_labels.
        """
        if self.data_obj is None:
            if not self.filename:
                return None
            self.data_obj = readData(self.filename, self.format)
        self.data_obj

    def __str__(self):
        if self.title:
            return self.title
        s, t = "", ""
        d = self.data_obj
        if d:
            s = "Character matrix in memory (id %d) with %d characters for %d taxa.  " % (id(d), d.n_char, d.n_tax)
        if self.filename:
            t =  "Characters from the file %s" % self.filename
        if s or t:
            return s + t
        return "None"

    def __deepcopy__(self, memo):
        #trees are expensive, so we don't make a deepcopy
        if self.data_obj:
            return DataSource(self.data_obj, taxon_labels=self.taxon_labels, title=self.title)
        elif self.filename:
            return DataSource(self.filename, format=self.format, title=self.title)
        return DataSource(None, taxon_labels=self.taxon_labels, title=self.title)

