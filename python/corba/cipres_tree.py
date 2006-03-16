#!/usr/bin/python
import sys, re
from omniORB import CORBA
import CosNaming, CipresIDL
from corba_util import *
from tree import Tree

numericTaxaInTreeRE = re.compile(r'\s*[(,]\s*(\d+)')
class CipresTree(Tree):
	def makeCorbaTree(self):
		global numericTaxaInTreeRE
		newick = self.newick + ';'
		return CipresIDL.Tree(newick, -1 , [int(i) for i in numericTaxaInTreeRE.findall(newick)], 'one')

