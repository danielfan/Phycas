#!/usr/bin/python
import sys, re
from omniORB import CORBA
import CosNaming, CipresIDL
from corba_util import *

def makeCorbaCharacters(matRow, symbols):
	return [symbols.index(i) for i in matRow]

def makeCorbaRawMatrix(mat, symbols):
	return [makeCorbaCharacters(row, symbols) for row in mat]

class CipresMatrix:
	def __init__(self, mat, symbols = 'ACGT'):
		self.matrix = mat
		self.symbols = symbols
	def getNChars(self):
		return len(self.matrix)> 0 and len(self.matrix[0]) or 0
	def getCorbaMatrix(self):
		rawMat = makeCorbaRawMatrix(self.matrix, self.symbols)
		return CipresIDL.DataMatrix(self.symbols, len(self.symbols), 1, self.nChars , rawMat)
	def __str__(self): return str(self.matrix)
	nChars = property(getNChars)
	