from phycas import *

model.type = 'jc'
jc = model()

model.type = 'hky'
hky = model()

model.type = 'gtr'
gtr = model()

partition.addSubset(range(1,100,3), jc, 'first')
partition.addSubset(range(2,100,3), hky, 'second')
partition.addSubset(range(3,100,3), gtr, 'third')
partition()
