from _ReadNexus import *
from phycas.DataMatrix._DataMatrix import DataMatrixWrapper
import re

_ROUND_TRIP_EVERY_NEXUS_READ = False

class TreeDescription(object):
    def __init__(self, nxs_full_tree_description):
        self.ftd = nxs_full_tree_description
    def getNameFromC(self):
        return self.ftd.getName()
    def getNewickFromC(self):
        return self.ftd.getNewick()
    def getIsRootedFromC(self):
        return self.ftd.isRooted()
    def writeNEXUSTreeCommand(self, out):
        c = self.rooted and "R" or "U"
        out.write("Tree %s = [&%s] %s;\n" % (escapeNexus(self.name), c, self.newick))
    def buildTree(self, tree=None):
        """Calls buildFromString or other appropriate method to construct a tree 
        in the the variable `tree`
        Returns the `tree` instance."""
        if tree is None:
            tree = Phylogeny.Tree()
        tree.buildFromString(self.newick, False) # tree descriptions from NCL 2.1 are 1-based not 0-based
        return tree
    name = property(getNameFromC)
    newick = property(getNewickFromC)
    rooted = property(getIsRootedFromC)

_tokenBreakers = re.compile(r"[(){}_\"-\[\]/\\,;:=*`+<>' \t\n]")
_quoteDemandingTokenBreakers = re.compile('[-\\]\\[/\\,;:=*`+<>\'\\t\\n_(){}\\"]')
def escapeNexus(s):
    global _tokenBreakers, _quoteDemandingTokenBreakers
    ls = len(s)
    if ls > 1:
        if len(_tokenBreakers.split(s)) > 1:
            if len(_quoteDemandingTokenBreakers.split(s)) > 1:
                return "'%s'" % "''".join(s.split("'"))
            return s.replace(" ", "_")
        return s
    if ls == 0:
        return "''"
    if (s in "'_[") or (s.isspace()):
        if s == "'":
            return "''''"
        return "'%s'" % s
    return s
    
class NexusReader(NexusReaderBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Need to write.
    
    """
    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Need to write.

        >>> from phycas import *
        >>> reader = ReadNexus.NexusReader()
        >>> reader.readFile('../Tests/Data/nyldna4.nex')
        >>> print reader.getNChar()
        3080

        """
        NexusReaderBase.__init__(self, -1)

    def readFile(self, fn):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Need to write.

        """
        NexusReaderBase.readFile(self, fn)
        if _ROUND_TRIP_EVERY_NEXUS_READ:
            import sys
            sys.stdout.write("==========================DEBUG OUTPUT FROM %s  ================\n" % (__file__))
            self.writeNEXUS(sys.stdout)
            sys.stdout.write("========================== END DEBUG OUTPUT FROM %s  ================\n" % (__file__))
        
    def getNChar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Need to write.

        """
        return NexusReaderBase.getNChar(self)
        
    def getErrorMessage(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Need to write.

        """
        return NexusReaderBase.getErrorMessage(self)
        
    def getTrees(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Need to write.

        """
        return [TreeDescription(t) for t in NexusReaderBase.getTrees(self)]
    
    def getLastDiscreteMatrix(self, gaps_to_missing=True):
        raw_mat = getLastDiscreteMatrix(self, gaps_to_missing)
        if raw_mat:
            return DataMatrixWrapper(raw_mat)
        return None


    def writeNEXUS(self, out, appending=False):
        taxa = self.getTaxLabels()
        et = [escapeNexus(t) for t in taxa]
        if et:
            if not appending:
                out.write("#NEXUS\n")
            out.write("Begin Taxa;\n Dimensions ntax = %d;\n" % len(et))
            out.write(" Taxlabels %s ;\nEnd;\n\n" % " ".join(et))
        dm = self.getLastDiscreteMatrix()
        if dm:
            out.write("Begin Characters;\n Dimensions nchar = %d;\n " % self.getNChar())
            out.write(dm.getNEXUSFormatCommand())
            out.write("\n Matrix\n")
            assert(len(et) == dm.n_tax)
            for ind, name in enumerate(et):
                out.write("%-20s "% name)
                r = dm.getRow(ind)
                for c in r:
                    out.write(dm.getCodeToSymbol(c))
                out.write("\n")
            out.write(";\nEnd;\n\n")
            
        tr = self.getTrees()
        if tr:
            out.write("Begin Trees;\n")
            for t in tr:
                t.writeNEXUSTreeCommand(out)
            out.write("End;\n\n")
            
