from phycas.Phylogeny import *

outf = file('out.txt', 'w')

z = '-**---****-----******------*'
outf.write('pattern = %s\n' % z)

s = Split();
s.createFromPattern(z)
outf.write('pattern = %s (after Split construction)\n' % s.createPatternRepresentation())

s.setExcluded([1,3,5,7,9])
outf.write('pattern = %s (after excluding 1,3,5,7,9)\n' % s.createPatternRepresentation())

on_list = s.getOnList()
outf.write('list of bits set:')
for onbit in on_list:
    outf.write(' %d' % onbit)
outf.write('\n')

off_list = s.getOffList()
outf.write('list of bits unset:')
for offbit in off_list:
    outf.write(' %d' % offbit)
outf.write('\n')

excl_list = s.getExcludedList()
outf.write('list of excluded bits:')
for exclbit in excl_list:
    outf.write(' %d' % exclbit)
outf.write('\n')

s.setOnSymbol('^')
s.setOffSymbol('~')
s.setExcludedSymbol('!')
outf.write('pattern = %s (on = %s, off = %s, excluded = %s)\n' % (s.createPatternRepresentation(), s.getOnSymbol(), s.getOffSymbol(), s.getExcludedSymbol()))

s.setOnSymbol('*')
s.setOffSymbol('-')
s.setExcludedSymbol('x')
outf.write('pattern = %s (on = %s, off = %s, excluded = %s)\n' % (s.createPatternRepresentation(), s.getOnSymbol(), s.getOffSymbol(), s.getExcludedSymbol()))

num_on_bits = s.countOnBits()
outf.write('pattern = %s (number of set bits is %d)\n' % (s.createPatternRepresentation(), num_on_bits))
num_off_bits = s.countOffBits()
outf.write('pattern = %s (number of unset bits is %d)\n' % (s.createPatternRepresentation(), num_off_bits))
s.invertSplit()
num_on_bits_inv = s.countOnBits()
outf.write('pattern = %s (after inverting, number of set bits now %d)\n' % (s.createPatternRepresentation(), num_on_bits_inv))
complexity = s.calcComplexity()
outf.write('Complexity is %d (minimum of %d and %d)\n' % (complexity, num_on_bits, num_off_bits))

# These have yet to be tested in this file
#
#setBit(self, b):
#setBits(self, btuple):
#unsetBit(self, b):
#unsetBits(self, btuple):
#isBitSet(self, b):
#createNewickRepresentation(self):
#equals(self, other_split):
#isCompatible(self, other_split):
#subsumedIn(self, other_split):

outf.close()
