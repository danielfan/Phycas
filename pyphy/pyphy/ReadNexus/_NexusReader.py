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

        >>> from ReadNexus import *
        >>> reader = NexusReader()
        >>> reader.readFile('nyldna4.nex')
        >>> print reader.getNChar()
        3080

        """
        NexusReaderBase.__init__(self, -1)

    def readFile(self, fn):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Need to write.

        """
        print fn
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
        
