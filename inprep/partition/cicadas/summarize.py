# Log of marginal likelihood (harmonic mean method) = -94684.183203
#  -95397.43634390 Stabilized Steppingstone Sampling (SSS) method

smallest_possible_lnL = -1000000.0

import re, glob, sys, math
nschemes = 4
smap = {}
def addToSchemeMap(scheme, method, mlike):
    if (scheme,method) in smap.keys():
        smap[(scheme,method)].append(mlike)
    else:
        smap[(scheme,method)] = [mlike]

if __name__ == '__main__':
    hm_results = False
    ss_results = False
    filenames = glob.glob('results2/*.sump.txt')
    filenames += glob.glob('results3/*.sump.txt')
    for fn in filenames:
        print 'Processing %s...' % fn
        m = re.match('(results2|results3)/([0-9]+)_([0-9]+)_[hms]+_([0-9]+)', fn)
        if m:
            pid    = int(m.group(2))
            scheme = int(m.group(3))
            rnseed = int(m.group(4))
            print '--> pid=%d, scheme=%d, rnseed=%d' % (pid,scheme,rnseed)
        else:
            pid = 0
            scheme = 0
            rnseed = 0
        lines = open(fn, 'r').readlines()
        line = lines[-1].strip()
        if 'SSS' in line:
            m = re.match('([0-9.-]+) Stabilized Steppingstone Sampling \(SSS\) method', line)
            if m:
                marglike = float(m.group(1))
                addToSchemeMap(scheme, 'ss', marglike)
                ss_results = True
            else:
                print 'ERROR: no match for last line in file name = %s' % fn
        else:
            m = re.match('Log of marginal likelihood \(harmonic mean method\) = ([0-9.-]+)', line)
            if m:
                marglike = float(m.group(1))
                addToSchemeMap(scheme, 'hm', marglike)
                hm_results = True
            else:
                print 'ERROR: no match for last line in file name = %s' % fn
    
    hm_mean = [0.0]*nschemes
    hm_sd = [0.0]*nschemes
    hm_max = [smallest_possible_lnL]*nschemes
    hm_min = [0.0]*nschemes
    ss_mean = [0.0]*nschemes
    ss_sd = [0.0]*nschemes
    ss_max = [smallest_possible_lnL]*nschemes
    ss_min = [0.0]*nschemes
    for i in range(nschemes):
        index = i
        print 'Partition scheme %d:' % i
        try:
            v = smap[(i,'hm')]
            for x in v:
                if x > hm_max[index]:
                    hm_max[index] = x
                if x < hm_min[index]:
                    hm_min[index] = x
            n = len(v)
            mean = sum(v)/float(n)
            if n > 1:
                sd = math.sqrt((sum([x**2.0 for x in v]) - float(n)*mean**2.0)/float(n-1))
            else:
                sd = None
            span = max(v) - min(v)
            hm_mean[index] = mean
            hm_sd[index] = sd
            print '  hm: %.5f (n = %d, max-min = %.5f, sd = %.5f)' % (mean, n, span, sd)
        except KeyError:
            print '  (no data)'
            
        try:
            v = smap[(i,'ss')]
            for x in v:
                if x > ss_max[index]:
                    ss_max[index] = x
                if x < ss_min[index]:
                    ss_min[index] = x
            n = len(v)
            mean = sum(v)/float(n)
            if n > 1:
                sd = math.sqrt((sum([x**2.0 for x in v]) - float(n)*mean**2.0)/float(n-1))
            else:
                sd = None
            span = max(v) - min(v)
            ss_mean[index] = mean
            ss_sd[index] = sd
            print '  ss: %.5f (n = %d, max-min = %.5f, sd = %.5f)' % (mean,n,span,sd)
        except KeyError:
            print '  (no data)'

    # determine minima and maxima for purposes of setting limits on y-axis
    hm_lower = [hm-sd for hm,sd in zip(hm_mean,hm_sd)]
    ss_lower = [ss-sd for ss,sd in zip(ss_mean,ss_sd)]
    minysd = min(hm_lower) < min(ss_lower) and min(hm_lower) or min(ss_lower)

    hm_upper = [hm+sd for hm,sd in zip(hm_mean,hm_sd)]
    ss_upper = [ss+sd for ss,sd in zip(ss_mean,ss_sd)]
    maxysd = max(hm_upper) > max(ss_upper) and max(hm_upper) or max(ss_upper)

    print
    print 'R code for line plot:'
    print '------------------------------------------------------------'
    
    # needed for errbar function
    print 'library(Hmisc)'
    
    # need fake x axis values to give to errbar
    print 'x <- c(1,2,3,4)'
    
    # save hm and ss along with their upper and lower error bar locations
    print 'hm <- c(%s)' % ','.join(['%.5f' % x for x in hm_mean])
    print 'hmsd <- c(%s)' % ','.join(['%.5f' % x for x in hm_sd])
    print 'hmupper <- c(%s)' % ','.join(['%.5f' % x for x in hm_max])
    print 'hmlower <- c(%s)' % ','.join(['%.5f' % x for x in hm_min])
    print 'ss <- c(%s)' % ','.join(['%.5f' % x for x in ss_mean])
    print 'sssd <- c(%s)' % ','.join(['%.5f' % x for x in ss_sd])
    print 'ssupper <- c(%s)' % ','.join(['%.5f' % x for x in ss_max])
    print 'sslower <- c(%s)' % ','.join(['%.5f' % x for x in ss_min])
    
    # plot points and error bars
    if hm_results:
        print 'plot(x, hm, type="l", lty="dotted", lwd=2, xaxt="n", xlab="Partition scheme", ylab="Marginal likelihood", ylim=c(%g,%g))' % (minysd,maxysd)
        print 'lines(x, ss, lty="solid", lwd=2)'
        print 'axis(side=1,at=x,labels=c("None","Gene","Codon","Both"))'
        print 'errbar(x, hm, hmupper, hmlower, add=TRUE, pch=".")'
        print 'errbar(x, ss, ssupper, sslower, add=TRUE, pch=".")'
    else:
        print 'plot(x, ss, type="l", lty="solid", lwd=2, xaxt="n", xlab="Partition scheme", ylab="Marginal likelihood", ylim=c(%g,%g))' % (minysd,maxysd)
        print 'axis(side=1, at=x, labels=c("None","Gene","Codon","Both"))'
        print 'errbar(x, ss, ssupper, sslower, add=TRUE, pch=".")'
    print '------------------------------------------------------------'

    #print
    #print 'R code for bar plot:'
    #print '--------------------'
    #print 'm <- matrix(data=c(hm,ss),ncol=4)'
    #print 'barplot(m,beside=TRUE,ylim=c(%g,%g),xpd=FALSE,names.arg=c("None","Gene","Codon","Both"))' % (miny,maxy)
    
