import sys
from phycas import *

if __name__ == '__main__':
	filename = getPhycasTestData('green.nex')
	blob = readFile(filename)
	like.data_source = blob.characters
	like.tree_source =  TreeCollection(filename='green.tre')
	for td in like.tree_source:
		n = td.getNewickFromC()
		print n
		t = Phylogeny.Tree(newick=n)
		t.debugListTree()
	like.starting_edgelen_dist = None
	model.type = 'jc'
	lnL = like()
	print "lnL = ", lnL
	

