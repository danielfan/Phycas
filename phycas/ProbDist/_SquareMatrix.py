from _ProbDistExt import *

class SquareMatrix(SquareMatrixBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapsulates a square matrix of floating point values (underlying C++
    implementation stores these as doubles).
        
    """
    def __init__(self, dimension, value):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Create a square matrix of size dimension containing value in every 
        cell.
        
        """
        SquareMatrixBase.__init__(self, dimension, value)
        
    def identity(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Converts existing matrix to an identity matrix (1s on diagonal, 0s 
        everywhere else). Dimension of the matrix is not changed.
        
        """
        SquareMatrixBase.identity(self)
    
    def trace(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the sum of the elements on the main diagonal.
        
        """
        return SquareMatrixBase.trace(self)
    
    def inverse(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a SquareMatrix that is the inverse of this matrix.
        
        """
        return SquareMatrixBase.inverse(self)
    
    def getDimension(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns an integer representing the number of rows of the matrix. The
        number of columns is the same value because this is a square matrix.
        
        """
        return SquareMatrixBase.getDimension(self)
        
    def getElement(self, i, j):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns (i,j)th element of the square matrix.
        
        """
        return SquareMatrixBase.getElement(self, i, j)
        
    def setElement(self, i, j, v):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets (i,j)th element of the square matrix to value v.
        
        """
        SquareMatrixBase.setElement(self, i, j, v)
        
    def getMatrix(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns square matrix in the form of a two-dimensional list.
        
        """
        dim = self.getDimension()
        v = SquareMatrixBase.getMatrix(self)
        m = []
        k = 0
        for i in range(dim):
            tmp = []
            for j in range(dim):
                tmp.append(v[k])
                k += 1
            m.append(tmp)
        return m
        
    def setMatrixFromFlattenedList(self, dim, v):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Replaces existing or creates a new square matrix using the supplied
        unidimensional list or tuple v. The supplied list v is expected to 
        have length equal to the square of dim, the number of elements in a
        single row or column of the matrix.
        
        """
        SquareMatrixBase.setMatrix(self, dim, v)
        
    def setMatrix(self, m):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Replaces existing or creates a new square matrix using the supplied
        two-dimensional list or tuple m.
        
        """
        dim = len(m[0])
        v = []
        for row in m:
            for col in row:
                v.append(col)
        SquareMatrixBase.setMatrix(self, dim, v)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Represents matrix as string.
        
        """
        s = SquareMatrixBase.__repr__(self)
        return s
        