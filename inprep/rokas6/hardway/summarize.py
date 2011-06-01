import math

# Read file into lines list
lines = open('data_first1000jc.txt','r').readlines()

# Create dictionary with 
#   key   = tree number (1 to 105) 
#   value = log marginal likelihood
marglikes = {}
for line in lines:
    parts = line.split()
    tree = int(parts[0])
    lnL = float(parts[1])
    marglikes[tree] = lnL
    
# Find maximum log marginal likelihood, which will be 
# factored out for numerical purposes
maxlnl = max(marglikes.values())

# Sum the ratio of each value to the maximum value
sum_diffs = 0.0
for m in marglikes.values():
    term = m - maxlnl
    sum_diffs += math.exp(term)
    
n = float(len(marglikes))
log_marg_like = maxlnl + math.log(sum_diffs) - math.log(n)
print
print 'Overall log marginal likelihood (based on %d distinct tree topologies):' % len(marglikes)
print '  log(marginal likelihood) =', log_marg_like

print
print 'Individual log marginal likelihoods:'
for i in range(105):
    t = i + 1
    if t in marglikes.keys():
        print '  %d\t%g' % (t, marglikes[t])
    else:
        print '  %d\t<not available>' % (t,)
        
print
print 'Marginal posterior distribution of tree topologies:'
log_total_marg_like = maxlnl + math.log(sum_diffs)
for i in range(105):
    t = i + 1
    if t in marglikes.keys():
        prob = math.exp(marglikes[t] - log_total_marg_like)
        print '  %d\t%g' % (t, prob)
    else:
        print '  %d\t<not available>' % (t,)
