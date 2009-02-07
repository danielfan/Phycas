import os,sys,math,random
from phycas import *
from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions

class VarianceZeroError(Exception):
    def __init__(self):
        self.message = 'Cannot calculate autocorrelation because variance is zero'
    def __str__(self):
        return self.message
    
class VarianceUndefinedError(Exception):
    def __init__(self):
        self.message = 'Cannot calculate standard deviation because sample size is less than 2'
    def __str__(self):
        return self.message
    
class InvalidNumberOfColumnsError(Exception):
    def __init__(self, nparts, nexpected, line_num):
        self.message = 'Number of values (%d) on line %d inconsistent with number of column headers (%d)' % (nparts, line_num, nexpected)
    def __str__(self):
        return self.message

class ParamSummarizer(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Summarizes parameter sample contained in the param file, which is one
    of the two files output from an mcmc analysis (the other being the 
    tree file).
    
    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes ParamSummarizer object.
        
        """
        CommonFunctions.__init__(self, opts)
        
    def interpolate(self, xx, x, y):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Given x and y, which define three reference points, (x[0],y[0]), 
        (x[1],y[1]) and (x[2],y[2]), find the point on the interpolating 
        polynomial corresponding to x-coordinate xx (scalar).
        
        """
        term0 = (xx - x[1])*(xx - x[2])/((x[0] - x[1])*(x[0] - x[2]))
        term1 = (xx - x[0])*(xx - x[2])/((x[1] - x[0])*(x[1] - x[2]))
        term2 = (xx - x[0])*(xx - x[1])/((x[2] - x[0])*(x[2] - x[1]))
        retval = term0*y[0] + term1*y[1] + term2*y[2]
        return retval
        
    def cumLagrange(self, which, x, y):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Given x and y, which define three reference points, (x[0],y[0]), 
        (x[1],y[1]) and (x[2],y[2]), find the integral under the interpolating
        polynomial for the first (which=1) or both (which=2) segments.
        
        """
        xx = x[which]
        
        x0 = x[0]
        x1 = x[1]
        x2 = x[2]
        
        y0 = y[0]
        y1 = y[1]
        y2 = y[2]
        
        psi0 = (x0 - x1)*(x0 - x2)
        psi1 = (x1 - x0)*(x1 - x2)
        psi2 = (x2 - x0)*(x2 - x1)
        
        xterm0 = (xx**3.0 - x0**3.0)/3.0
        xterm1 = (xx**2.0 - x0**2.0)/2.0
        xterm2 = xx - x0
        
        term0 = xterm0*(y0/psi0 + y1/psi1 + y2/psi2)
        term1 = xterm1*(y0*(x1 + x2)/psi0 + y1*(x0 + x2)/psi1 + y2*(x0 + x1)/psi2)
        term2 = xterm2*(y0*x1*x2/psi0 + y1*x0*x2/psi1 + y2*x0*x1/psi2)
        
        cum = term0 - term1 + term2
        return -cum
        
    def ps_simpsons(self, betas, means):
        """
        This approach uses Simpson's method to interpolate between beta values using the
        interpolation polynomial in Lagrange form. Simpson' method is described in most
        calculus textbooks, and (as used here) fits a parabola to every three consecutive
        points, using the area under the parabola as an approximation to the integral.
        This approach provides two estimates for the integral corresponding to each segment
        except for the first and last. We use the simple average of the two estimates when
        they are available.
        
        """
        nbetas = len(betas)
        if nbetas < 3:
            raise Exception("Must have at least 3 beta values to compute path sampling using Simpson's rule")
        marginal_like = 0.0
        for i in range(nbetas - 2):
            x = [betas[i], betas[i+1], betas[i+2]]
            y = [means[i], means[i+1], means[i+2]]
            # print '\ni = %d (dx = %.5f, dy = %.5f)' % (i,x[0] - x[2], y[0] - y[2])
            # print 'x%d <- c(%s)' % (i,','.join(['%.6f' % z for z in x]))
            # print 'y%d <- c(%s)' % (i,','.join(['%.6f' % z for z in y]))
            # print 'plot(x%d,y%d,type="b")' % (i,i)
            if i == 0:
                # no averaging on the first segment
                a = self.cumLagrange(1, x, y)
                marginal_like += a
                b = self.cumLagrange(2, x, y)
                marginal_like += (b - a)/2.0
            elif i == nbetas - 3:
                # no averaging on the last segment
                a = self.cumLagrange(1, x, y)
                marginal_like += a/2.0
                b = self.cumLagrange(2, x, y)
                marginal_like += (b - a)
            else:
                # average two estimates for each of the middle segments
                a = self.cumLagrange(1, x, y)
                b = self.cumLagrange(2, x, y)
                marginal_like += b/2.0
                # if i == 7 or i == 8:
                #     xxx = []
                #     yyy = []
                #     for z in range(21):
                #         xxx_value = x[0] + (x[2] - x[0])*float(z)/20.0
                #         yyy_value = self.interpolate(xxx_value, x, y)
                #         xxx.append(xxx_value)
                #         yyy.append(yyy_value)
                #     print 'xx%d <- c(%s)' % (i,','.join(['%.5f' % z for z in xxx]))
                #     print 'yy%d <- c(%s)' % (i,','.join(['%.5f' % z for z in yyy]))
                #     print 'plot(xx%d,yy%d,type="b")' % (i,i)
                #     print 'a     = %.5f' % a
                #     print 'b - a = %.5f' % (b - a)
                #     print '(x[0] - x[1])*(y[0] + y[1])/2 = ',(x[0] - x[1])*(y[0] + y[1])/2.0
                #     print '(x[1] - x[2])*(y[1] + y[2])/2 = ',(x[1] - x[2])*(y[1] + y[2])/2.0
        self.output(" %.8f Path sampling method (using Simpson's rule)" % (marginal_like))

    def ps_trapezoid(self, betas, means):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This method approximates the integral under the curve defined by betas
        (x-coordinates) and means (y-coordinates) using the trapezoid method
        (straight line interpolation). This is the method advocated by in the
        Lartillot and Phillippe (2006) paper that introduced the thermodynamic
        integration method to phylogenetics.
        
        """
        nbetas = len(betas)
        marginal_like = 0.0
        for i in range(nbetas):
            if i == 0:
                before = betas[0]
            else:
                before = betas[i-1]
            if i == nbetas - 1:
                after = betas[nbetas - 1]
            else:
                after = betas[i+1]
            diff = before - after
            if diff < 0.0:
                raise Exception('Phycas does not currently support path sampling from prior toward posterior')
            #print 'mean %d = %f (diff = %f)' % (i, means[i], diff)
            marginal_like += means[i]*diff/2.0
        self.output(' %.8f Path sampling method (using trapezoid rule)' % (marginal_like))
                
    def ss(self, betas, likes):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This method estimates the marginal likelihood using the product of 
        ratios (the stepping stones) bridging the gap between the posterior
        and the prior. Each ratio is estimated using importance sampling, with
        the importance distribution being the power posterior defined by 
        the smaller of the two beta values in the ratio. See Xie et al. (2009;
        Systematic Biology; submitted Jan. 2009) for details.
        
        """
        # Calculate marginal likelihood using steppingstone sampling method
        # betas is a list of beta values ordered from the first sampled to the last sampled
        # likes is a map: i.e. likes[b] holds list of log-likelihoods sampled for beta value b
        # Assumes that betas[0] = 1.0 and betas[-1] = 0.0
        lnR = 0.0
        seR = 0.0
        nbetas = len(betas)
        if betas[0] != 1.0 or betas[-1] != 0.0:
            raise Exception('Steppingstone sampling requires beta values to be ordered from 1.0 (first) to 0.0 (last)')
        for i in range(1,nbetas):
            blarger = betas[i - 1]
            bsmaller = betas[i]
            beta_incr = blarger - bsmaller
            loglikes = likes[bsmaller]
            tmp = 0.0
            n = len(loglikes)
            Lmax = max(loglikes)
            for lnL in loglikes:
                tmp += math.exp(beta_incr*(lnL - Lmax))
            tmp /= float(n)
            lnR += beta_incr*Lmax + math.log(tmp)
        
            # standard error calculation
            tmp1 = 0.0
            for lnL in loglikes:
                aa = math.exp(beta_incr*(lnL - Lmax))/tmp
                tmp1 += math.pow((aa - 1.0),2.0)
            tmp1 /= float(n)
            seR += tmp1/float(n)
        self.output(' %.8f Steppingstone sampling method (se = %.8f)' % (lnR, seR))
        
    def autocorr_ess(self, values):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes the lag-1 autocorrelation (r) for the supplied values list.
        Also computes the effective sample size (ess). Returns the tuple 
        (r,ess).
        
        """
        nvalues = len(values)
        n = float(nvalues)
        
        # calculate the mean
        m = sum(values)/n
        
        # calculate the variance
        ss = 0.0
        for x in values:
            ss += x**2.0
        var = (ss - n*m*m)/(n - 1.0)
        
        # calculate the covariance
        cov = 0.0
        for i in range(nvalues - 1):
            x = values[i] - m
            y = values[i+1] - m
            cov += x*y
        cov /= n - 1.0
        if var <= 0.0:
            raise VarianceZeroError()
        r = cov/var
        ess = n*(1.0 - r)/(1.0 + r)
        return (r,ess)
    
    def marginal_likelihood(self, headers, lines, burnin):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Estimates the marginal likelihood for the path sampling/steppingstone
        sampling case and outputs the autocorrelation and effective sample 
        size for each parameter/beta combination. The supplied list headers 
        holds the column headers from the param file. The parameter names are
        taken from this. The supplied lines is a list of lines from the param
        file. The supplied burnin indicates the number of initial lines to 
        ignore (usually will be 1 because the starting values of all 
        parameters is always output and should always be skipped).
        
        """
        # Create params such that, for example, params['lnL'][1.0] is a list of 
        # all log-likelihood values sampled for beta = 1.0 (in the order in which
        # they were sampled)
        betas = []
        params = {}
        for h in headers:
            params[h] = {}
        start = 2 + burnin
        curr_beta = None
        for i,line in enumerate(lines[start:]):
            parts = line.split()
            if len(parts) != len(headers):
                raise InvalidNumberOfColumnsError(len(parts), len(headers), i + start + 1)
            beta = float(parts[1])
            if curr_beta is None or curr_beta != beta:
                curr_beta = beta
                betas.append(beta)
                for h,x in zip(headers,parts):
                    params[h][beta] = [float(x)]
            else:
                for h,x in zip(headers,parts):
                    params[h][beta].append(float(x))
                    
        # Output first-order autocorrelation for each parameter (and the log-likelihood)
        # for each beta value separately
        self.output('\nAutocorrelations (lag 1):\n')
        self.output('%12s%s' % ('beta ->',' '.join(['%12.5f' % b for b in betas])))
        for h in headers[2:]:   # skip 'Gen' and 'beta' headers
            s = []
            for b in betas:
                p = params[h][b]
                try:
                    r,ess = self.autocorr_ess(p)
                except VarianceZeroError:
                    s.append('%12s' % '---')
                else:
                    s.append('%12.5f' % r)
            self.output('%12s%s' % (h,' '.join(s)))

        # Output effective sample size for each parameter (and the log-likelihood)
        # for each beta value separately
        actual_sample_size = len(params['LnL'][1.0])
        all_same = True
        self.output('\nEffective sample sizes (actual sample size = %d):\n' % actual_sample_size)
        self.output('%12s%s' % ('beta ->',' '.join(['%12.5f' % b for b in betas])))
        for h in headers[2:]:   # skip 'Gen' and 'beta' headers
            s = []
            for b in betas:
                p = params[h][b]
                n = len(p)
                if n != actual_sample_size:
                    all_same = False
                try:
                    r,ess = self.autocorr_ess(p)
                except VarianceZeroError:
                    s.append('%12s' % '---')
                else:
                    s.append('%12.1f' % ess)
            self.output('%12s%s' % (h,' '.join(s)))
        if not all_same:
            self.output('\nWarning: actual sample sizes varied across beta values and/or parameters!')

        # Compute means of log-likelihoods for each beta value (used for ps calculation)
        means = []
        for b in betas:
            p = params['LnL'][b]
            m = sum(p)/float(len(p))
            means.append(m)
        self.output('\nMean log-likelihood for each value of beta used\nfor marginal likelihood estimation:\n')
        self.output('%12s %12s' % ('beta','mean lnL'))
        for b,m in zip(betas, means):
            s = ''
            if b == 1.0:
                s = '(posterior)'
            elif b == 0.0:
                s = '(prior)'
            self.output('%12.5f %12.5f %s' % (b,m,s))

        self.output('\nMarginal likelihood estimates:')
        try:
            self.ss(betas, params['LnL'])
        except Exception,e:
            self.output(' %s' % e.message)
        try:
            self.ps_trapezoid(betas, means)
        except Exception,e:
            self.output(' %s' % e.message)
        try:
            self.ps_simpsons(betas, means)
        except Exception,e:
            self.output(' %s' % e.message)
            
    def harmonic_mean(self, v):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calculate marginal likelihood using the harmonic mean method using 
        log-likelihoods supplied in v.
        
        """
        nignored = 0
        n = len(v)
        min_lnL = min(v)
        sum_diffs = 0.0
        for lnl in v:
            diff = lnl - min_lnL
            if diff < 500.0:
                sum_diffs += math.exp(-diff)
            else:
                nignored += 1
        log_harmonic_mean = math.log(n) + min_lnL - math.log(sum_diffs)
        if nignored > 0:
            self.warning('ignoring %d sampled log-likelihoods in harmonic mean calculation' % nignored)
        self.output('Log of marginal likelihood (harmonic mean method) = %f' % log_harmonic_mean)
            
    def summary_stats(self, v, cutoff=95):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes the following summary statistics for the supplied vector v:
        first-order autocorrelation (lag=1), effective sample size, lower 
        credible interval bound, upper credible interval bound, minimum, 
        maximum, sample mean, and sample standard deviation (i.e. divide by 
        n-1). The value of cutoff is the percentage to use for the credible
        interval. If v is None, returns tuple of header strings. If v is not
        None, returns tuple of summary statistics.
        
        """
        if v is None:
            h = ('autocorr', 'ess', 'lower %d%%' % int(cutoff), 'upper %d%%' % int(cutoff), 'min', 'max', 'mean', 'stddev')
            return h
            
        if len(v) < 2:
            raise VarianceUndefinedError()
        s = []
        
        # compute autocorr and ess
        try:
            r,ess = self.autocorr_ess(v)
        except VarianceZeroError:
            raise VarianceZeroError()
        else:
            s.extend((r,ess))
            
        # compute lower and upper
        v.sort()
        n = float(len(v))
        p = float(cutoff)/100.0
        lower_at = int(math.ceil(n*(1 - p)))
        upper_at = int(math.ceil(n*p))
        lower = v[lower_at]
        upper = v[upper_at]
        s.extend((lower, upper))
        
        # compute min and max
        s.extend((v[0], v[-1]))
        
        # compute mean and stddev
        mean = sum(v)/n
        ss = sum([x**2 for x in v])
        var = (ss - n*mean**2)/(n - 1.0)
        sd = math.sqrt(var)
        s.extend((mean, sd))
        
        return tuple(s)
    
    def std_summary(self, headers, lines, burnin):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Produces a table of summary statistics for each parameter using the
        data in a param file. This is the companion to marginal_likelihood for
        the standard mcmc case. The supplied list headers holds the column 
        headers from the param file. The parameter names are taken from this.
        The supplied lines is a list of lines from the param file. The 
        supplied burnin indicates the number of initial lines to ignore 
        (usually will be 1 because the starting values of all parameters is 
        always output and should always be skipped). See documentation for
        the summary_stats function for a list of summary statistics computed.
        
        """
        # Create params such that, for example, params['lnL'] is a list of 
        # all log-likelihood values (in the order in which they were sampled)
        params = {}
        for h in headers:
            params[h] = []
        start = 2 + burnin
        for i,line in enumerate(lines[start:]):
            parts = line.split()
            if len(parts) != len(headers):
                raise InvalidNumberOfColumnsError(len(parts), len(headers), i + start + 1)
            for h,x in zip(headers,parts):
                params[h].append(float(x))
                
        # Output summary statistics for each parameter (and the log-likelihood)
        self.output('\nSummary statistics:\n')
        stats_headers = ('param','n') + self.summary_stats(None)
        sz = len(stats_headers)
        sss = '%12s'*sz + '\n'
        self.output(sss % stats_headers)
        for h in headers[1:]:   # skip 'Gen' header
            v = params[h]
            try:
                stats = (h,len(v)) + self.summary_stats(v)
                sss = '%12s' + '%12d' + '%12.5f'*(sz - 2)
                self.output(sss % stats)
            except (VarianceZeroError,VarianceUndefinedError):
                sss = '%12s'*sz + '\n'
                sub = tuple(['---']*sz)
                self.output(sss % sub)
        self.output()
        self.harmonic_mean(params['LnL'])
        
    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reads the contents of the param file and decides whether to use the 
        marginal_likelihood or std_summary functions to summarize the data
        therein. The marginal_likelihood function is called if the second
        header is "beta"; otherwise, std_summary is called.
        
        """
        fn = self.opts.file
        burnin = self.opts.burnin
        lines = open(fn, 'r').readlines()
        if len(lines) < 3 + burnin:
            self.output("File '%s' does not look like a parameter file (too few lines)")
        else:
            headers = lines[1].split()
            if headers[1] == 'beta':
                try:
                    self.marginal_likelihood(headers, lines, burnin)
                except InvalidNumberOfColumnsError, e:
                    print e
            else:
                self.std_summary(headers, lines, burnin)
