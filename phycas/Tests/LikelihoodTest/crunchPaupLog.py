import re,math

stuff = open('paup.log', 'r').read()

# Tree             1
# -ln L     7710.492
paup_results = re.findall('Tree\s+1\s+-ln L\s+([0-9.]+)', stuff, re.MULTILINE)

# Phycas GTR+G+I lnL = -7710.49182
phycas_results = re.findall('Phycas ([CFGHIJKTRY018+]+) lnL = ([0-9.-]+)', stuff, re.MULTILINE)

print '%12s\t%12s\t%12s\t%12s' % ('model', 'paup', 'phycas', 'fabs(diff)')
for a,(nm,b) in zip(paup_results,phycas_results):
    paup_lnL = -float(a)
    phycas_lnL = float(b)
    diff_lnL = math.fabs(paup_lnL - phycas_lnL)
    print '%12s\t%12.5f\t%12.5f\t%12.5f' % (nm, paup_lnL, phycas_lnL, diff_lnL)


