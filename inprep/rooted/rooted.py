import sys
from phycas import *

if __name__ == '__main__':
	#filename = getPhycasTestData('green.nex')
	#filename = 'test.nex'
	filename = 'junk.nex'
	blob = readFile(filename)
	like.data_source = blob.characters
	like.tree_source = TreeCollection(filename='test.tre')
	for t in like.tree_source:
		t.debugListTree()
	like.starting_edgelen_dist = None
	model.type = 'jc'
	#model.curr()
	#print model.state_freqs
	#like.store_site_likes = True
	lnL = like()
	#print like.site_likes
	print "lnL = ", lnL
	

