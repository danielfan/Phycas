from _DataMatrixBase import *

class DataMatrix(DataMatrixBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapsulates a data matrix used for phylogenetic analysis. It is not
    possible to construct an empty DataMatrix object; the way to obtain
    one is to import ReadNexus and call its getDiscreteMatrix function:

    >>> from phypy import *
    >>> r = ReadNexus.NexusReader()
    >>> r.readFile('../Tests/Data/nyldna4.nex')
    >>> m = ReadNexus.getDiscreteMatrix(r, 0)

    """
    def getNChar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of characters in the data matrix.

        >>> from phypy import *
        >>> r = ReadNexus.NexusReader()
        >>> r.readFile('../Tests/Data/nyldna4.nex')
        >>> m = ReadNexus.getDiscreteMatrix(r, 0)
        >>> m.getNChar()
        3080

        """
        DataMatrixBase.getNChar(self)

    def getNTax(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of taxa in the data matrix.

        >>> from phypy import *
        >>> r = ReadNexus.NexusReader()
        >>> r.readFile('../Tests/Data/nyldna4.nex')
        >>> m = ReadNexus.getDiscreteMatrix(r, 0)
        >>> m.getNTax()
        4
        
        """
        DataMatrixBase.getNTax(self)

    def getNStates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of states associated with the data type of the
        matrix. Note: there may well be more state codes than this, because
        ambiguities are coded as separate states.

        >>> from phypy import *
        >>> r = ReadNexus.NexusReader()
        >>> r.readFile('../Tests/Data/nyldna4.nex')
        >>> m = ReadNexus.getDiscreteMatrix(r, 0)
        >>> m.getNStates()
        4

        """
        DataMatrixBase.getNStates(self)

    def getSymbolsList(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string showing the symbols found in the data matrix. These
        original symbols or ambiguity/polymorphism specifications have been
        recoded, and getStateList() returns the codes used.

        >>> from phypy import *
        >>> r = ReadNexus.NexusReader()
        >>> r.readFile('../Tests/Data/nyldna4.nex')
        >>> m = ReadNexus.getDiscreteMatrix(r, 0)
        >>> print m.getSymbolsList()
        ACGT?N

        """
        DataMatrixBase.getSymbolsList(self)

    def getStateList(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a tuple representing the state list. The state list is a
        vector of integers containing all relevant information about each of
        the states found in the data matrix. Each state discovered is encoded
        as an integer, and each ambiguity or polymorphism found is recorded
        as if it were a separate state. The state list can be used to uncover
        the original meaning of these coded states.

        >>> from phypy import *
        >>> r = ReadNexus.NexusReader()
        >>> r.readFile('../Tests/Data/nyldna4.nex')
        >>> m = ReadNexus.getDiscreteMatrix(r, 0)
        >>> print m.getStateList()
        (1, 0, 1, 1, 1, 2, 1, 3, 5, -1, 0, 1, 2, 3, 4, 0, 1, 2, 3)

        Here is a translation of the state list in the above example. In this
        case, the states encountered in the original data file (nyldna4.nex)
        were A, C, G, T, ?, N. 

        state   state
        list    list
        index   element  meaning
        ------------------------------------------------------------
           0       1     next 1 element defines next state (A)
           1       0     first state has code 0
           2       1     next 1 element defines next state (C)
           3       1     second state has code 1
           4       1     next 1 element defines next state (G)
           5       2     third state has code 2
           6       1     next 1 element defines next state (T)
           7       3     fourth state has code 3
           8       5     next 5 elements define next state (?)
           9      -1     fifth state includes gaps
          10       0     fifth state includes state 0
          11       1     fifth state includes state 1
          12       2     fifth state includes state 2
          13       3     fifth state includes state 3
          14       4     next 4 elements define next state (N)
          15       0     sixth state includes state 0
          16       1     sixth state includes state 1
          17       2     sixth state includes state 2
          18       3     sixth state includes state 3
        ------------------------------------------------------------

        """
        DataMatrixBase.getStateList(self)

    def getStateListPos(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a tuple showing the index corresponding to each state found
        in the data matrix into the state list.

        >>> from phypy import *
        >>> r = ReadNexus.NexusReader()
        >>> r.readFile('../Tests/Data/nyldna4.nex')
        >>> m = ReadNexus.getDiscreteMatrix(r, 0)
        >>> print m.getStateListPos()
        (0, 2, 4, 6, 8, 14)

        See first column of table in the documentation for getStateList()
        function to see what the state list values mean.

        """
        DataMatrixBase.getStateListPos(self)


