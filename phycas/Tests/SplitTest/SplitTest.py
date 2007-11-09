from phycas.Phylogeny import *

outf = file('out.txt', 'w')

z = '-**---****-----******------*'
outf.write('pattern = %s\n' % z)

s = Split();
s.createFromPattern(z)
s.setExcluded([1,3,5,7,9])

outf.write('pattern = %s\n' % s.createPatternRepresentation())

on_list = s.getOnList()
outf.write('on list:')
for onbit in on_list:
    outf.write(' %d' % onbit)
outf.write('\n')

off_list = s.getOffList()
outf.write('off list:')
for offbit in off_list:
    outf.write(' %d' % offbit)
outf.write('\n')

excl_list = s.getExcludedList()
outf.write('excluded list:')
for exclbit in excl_list:
    outf.write(' %d' % exclbit)
outf.write('\n')

outf.close()
