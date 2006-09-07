from _ReadNexus import *

class NexusReader(NexusReaderBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Need to write.
    
    """
    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Need to write.

        >>> from phypy import *
        >>> reader = ReadNexus.NexusReader()
        >>> reader.readFile('../Examples/Data/nyldna4.nex')
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
        return NexusReaderBase.getTrees(self)
        
