import math, sys
try:
    fn = sys.argv[1]
except:
    fn = 'data.txt'

lines = open(fn,'r').readlines()

marglikes = {}
for line in lines:
   parts = line.split()
   tree = int(parts[0])
   lnL = float(parts[1])
   marglikes[tree] = lnL

mv = marglikes.values()
maxlnl = max(mv)
sum_diffs = 0.0
for m in marglikes.values():
   term = m - maxlnl
   sum_diffs += math.exp(term)

n = float(len(marglikes))
log_marg_like = maxlnl + math.log(sum_diffs) - math.log(n)
print 'log(marginal likelihood) =', log_marg_like
print 'based on %d distinct tree topologies' % len(marglikes)
log_total_marg_like = maxlnl + math.log(sum_diffs)
print 'posterior probs'
tot = 0.0
for i in range(len(mv)):
   t = i + 1
   if t in marglikes.keys():
       prob = math.exp(marglikes[t] - log_total_marg_like)
       tot += prob
       print '%d\t%g' % (t, prob)
   else:
       print '%d\t<not available>' % (t,)
print 'Total prob =', tot
