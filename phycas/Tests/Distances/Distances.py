#   #!/usr/bin/env python

import sys
from phycas import *
#from phycas.Constants import *

print sys.prefix


# def calcDistanceMatrix():
# 	matrix = []
# 	
# 	for i in range(p.ntax):
# 		print "%-20s" % p.taxon_labels[i]
# 		dataRow = p.data_matrix.getRow(i)
# 		for d in dataRow:
# 			print d,
# 		matrix.append(dataRow)
# 	
# 	distance = {}
# 	
# 	for i in range(p.ntax):
# 		for j in range(i + 1,p.ntax):
# 			diff = 0
# 			for k in range(p.nchar):
# 				xik = matrix[i][k]
# 				xjk = matrix[j][k]
# 				if xik != xjk:
# 					diff += 1
# 			distance[i,j] = diff
# 			
# 	for i,j in distance.keys():
# 		print "%5d%5d%10d" % (i, j, distance[i,j])

p = Phycas()
p.data_file_name = "../Data/green.nex"
p.default_model = 'gtr'

#p.readDataFromFile()
# p.calcDistances()
p.samc()

#calcDistanceMatrix()
