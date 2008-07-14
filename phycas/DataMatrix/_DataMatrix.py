from _DataMatrixBase import *

class DataMatrix(DataMatrixBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapsulates a data matrix used for phylogenetic analysis. It is not
    possible to construct an empty DataMatrix object; the way to obtain
    one is to import ReadNexus and call its getLastDiscreteMatrix function:

    >>> from phycas import *
    >>> r = ReadNexus.NexusReader()
    >>> r.readFile('../Tests/Data/nyldna4.nex')
    >>> m = r.getLastDiscreteMatrix(True)

    """
    def getNChar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of characters in the data matrix.

        >>> from phycas import *
        >>> r = ReadNexus.NexusReader()
        >>> r.readFile('../Tests/Data/nyldna4.nex')
        >>> m = r.getLastDiscreteMatrix(True)
        >>> m.getNChar()
        3080

        """
        DataMatrixBase.getNChar(self)

    def getNTax(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of taxa in the data matrix.

        >>> from phycas import *
        >>> r = ReadNexus.NexusReader()
        >>> r.readFile('../Tests/Data/nyldna4.nex')
        >>> m = r.getLastDiscreteMatrix(True)
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

        >>> from phycas import *
        >>> r = ReadNexus.NexusReader()
        >>> r.readFile('../Tests/Data/nyldna4.nex')
        >>> m = r.getLastDiscreteMatrix(True)
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

        >>> from phycas import *
        >>> r = ReadNexus.NexusReader()
        >>> r.readFile('../Tests/Data/nyldna4.nex')
        >>> m = r.getLastDiscreteMatrix(True)
        >>> print m.getSymbolsList()
        ['A', 'C', 'G', 'T', '?', 'B', 'D', 'H', 'K', 'M', 'R', 'S', 'V', 'W', 'Y']

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

        >>> from phycas import *
        >>> r = ReadNexus.NexusReader()
        >>> r.readFile('../Tests/Data/nyldna4.nex')
        >>> m = r.getLastDiscreteMatrix(True)
        >>> print m.getStateList()
        (1, 0, 1, 1, 1, 2, 1, 3, 4, 0, 1, 2, 3, 3, 1, 2, 3, 3, 0, 2, 3, 3, 0, 1, 3, 2, 2, 3, 2, 0, 1, 2, 0, 2, 2, 1, 2, 3, 0, 1, 2, 2, 0, 3, 2, 1, 3)

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

        >>> from phycas import *
        >>> r = ReadNexus.NexusReader()
        >>> r.readFile('../Tests/Data/nyldna4.nex')
        >>> m = r.getLastDiscreteMatrix(True)
        >>> print m.getStateListPos()
        (0, 2, 4, 6, 8, 13, 17, 21, 25, 28, 31, 34, 37, 41, 44)

        See first column of table in the documentation for getStateList()
        function to see what the state list values mean.

        """
        DataMatrixBase.getStateListPos(self)


class DataMatrixWrapper(object):
    """Wraps a DataMatrixBase object and enables access to field. Note that
    this class may be expanded later to allow editing, but (given that it 
    is tedious to recode a matrix's underlying data structure via python)
    only editing of the cell (no new state combinations will be introduced)."""
    DNA_Datatype = 0
    RNA_Datatype = 1
    AA_Datatype = 2 
    Codon_Datatype = 3
    Generic_Datatype = 4
    NEXUS_DATATYPE_NAMES = ("DNA", "RNA", "Protein", "Standard", "Standard")
    def __init__(self, dataMatrixBaseObj):
        assert(dataMatrixBaseObj)
        self.mat = dataMatrixBaseObj
        self.n_states = dataMatrixBaseObj.getNStates()
        raw_symbols_list = dataMatrixBaseObj.getSymbolsList()
        assert(len(raw_symbols_list) >= self.n_states)
        self.symbols = []
        self.symbols_to_code = {}
        
        for n in range(self.n_states):
            i = raw_symbols_list[n]
            if i != " ":
                self.symbols.append(i)
                self.symbols_to_code[i] = n
        u = 0
        for n, i in enumerate(raw_symbols_list):
            if i != " ":
                if n >= self.n_states:
                    self.symbols.append(i)
                    self.symbols_to_code[i] = n
            else:
                s = getListOfStates()
                if len(s) == 1:
                    sn = str(u)
                    while sn in self.symbols_to_code:
                        u += 1
                        sn = str(u)
                    self.symbols.append(sn)
                    self.symbols_to_code[sn] = n
                else:
                    s.sort()
                    try:
                        sym = [self.symbols[j] for j in s]
                    except KeyError:
                        # if this happens then we have a multi-state code that 
                        # refers to a higher state code -- all of the 
                        # "fundamental" states are supposed to come firet
                        assert(false) 
                    sym.sort()
                    t = tuple(s)
                    self.symbols.append("{%s}"% " ".join(sym))
                    self.symbols_to_code[t] = n
        self._state_list_pos = None
        self._state_list = None
        self._n_char = None
        self._n_tax = None
        self._datatype_enum = None
        self._int_wts = None
        self._float_wts = None

    def getDatatype(self):
        if self._datatype_enum is None:
            self._datatype_enum = self.mat.getDatatype()
        return self._datatype_enum
    datatype_enum = property(getDatatype)
    
    def getNChar(self):
        if self._n_char is None:
            self._n_char = self.mat.getNChar()
        return self._n_char
    n_char = property(getNChar)
    
    def getNTax(self):
        if self._n_tax is None:
            self._n_tax = self.mat.getNTax()
        return self._n_tax
    n_tax = property(getNTax)
    
    def getNStates(self):
        return self.n_states

    def getSymbolsList(self):
        return self.symbols

    def getStateList(self):
        if self._state_list is None:
            self._state_list = self.mat.getStateList()
        return self._state_list
    state_list = property(getStateList)

    def getStateListPos(self):
        if self._state_list_pos is None:
            self._state_list_pos = self.mat.getStateListPos()
        return self._state_list_pos
    state_list_pos = property(getStateListPos)

    def getIntWts(self):
        if self._int_wts is None:
            self._int_wts = self.mat.getIntWeights()
        return self._int_wts
    int_wts = property(getIntWts)

    def getFloatWts(self):
        if self._float_wts is None:
            self._float_wts = self.mat.getFloatWeights()
        return self._float_wts
    float_wts = property(getFloatWts)

    def getCodedDataMatrix(self):
        return self.mat.getCodedDataMatrix()

    def getRow(self, index):
        return self.mat.getRow(index)

    def getListOfStates(self, intCode):
        """Converts a single integer code (the codes in a matrix) to list of 
        integer codes."""
        if intCode == -1:
            return [-1]
        s = []
        index = self.state_list_pos[intCode]
        ns = self.state_list[index]
        for i in range(ns):
            index += 1
            s.append(self.state_list[index])
        return s

    def getSymbolToCode(self, sym):
        try:
            return self.symbols_to_code[sym]
        except:
            pass
        if isinstance(sym, str):
            if len(sym) == 1:
                return self.symbols_to_code[sym]
            if sym.startswith("{") and sym.endswith("}"):
                symList = symList[1:-1].split(" ")
                if len(symList) == 1:
                     return self.symbols_to_code[symList[0]]
            else:
                return self.symbols_to_code[sym]
        else:
            symList = list(sym)
        symList.sort()
        symT = tuple(symList)
        return self.symbols_to_code[symT]

    def getCodeToSymbol(self, intCode):
        if intCode == -1:
            return '-'
        return self.symbols[intCode]

    def getNEXUSFormatCommand(self):
        return "Format datatype=%s missing = ? gap = - ;" % DataMatrixWrapper.NEXUS_DATATYPE_NAMES[self.datatype_enum]
