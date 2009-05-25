#!/usr/bin/env python
import sys
import subprocess
import re
lnLPat = re.compile(r"^[0-9]+\s+(-[.0-9]+)\s+.*")
f = open("paupStdOut", "w")
subprocess.call(["paup", "-n", "debug-with-paup.nex"], stdout=f)
f.close()

expectedLnLPattern = re.compile("phycas lnL = (-[.0-9]+) .*")
f = open("paupStdOut", "rU")
expectedLnL = []
for line in f:
    m = expectedLnLPattern.search(line)
    if m:
        expectedLnL.append(float(m.groups(1)[0]))
f.close()
if not expectedLnL:
    sys.exit("No expected lnL found in comments in the output of PAUP")

paupLnLPattern = re.compile("^-ln L\s+([.0-9]+)")
f = open("paupStdOut", "rU")
for line in f:
    m = paupLnLPattern.match(line)
    if m:
        paupLnL = -1.0*float(m.groups(1)[0])
        try:
            phycasLnL = expectedLnL.pop(0)
        except:
            sys.exit("Found more PAUP LnL values than Phycas LnL values")
        if abs(paupLnL - phycasLnL) > 10e-4:
            sys.exit("%f != %f" %(paupLnL, phycasLnL))
if expectedLnL:
    sys.exit("Found more Phycas LnL values than PAUP LnL values")
print "It worked!"    
sys.exit(0)
