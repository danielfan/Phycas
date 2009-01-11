import os,sys,math,random
from phycas import *
from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions

class ParamSummarizer(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Summarizes parameter sample.
    
    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes ParamSummarizer object.
        
        """
        CommonFunctions.__init__(self, opts)
        
    def interpolate(self, xx, x, y):
        term0 = (xx - x[1])*(xx - x[2])/((x[0] - x[1])*(x[0] - x[2]))
        term1 = (xx - x[0])*(xx - x[2])/((x[1] - x[0])*(x[1] - x[2]))
        term2 = (xx - x[0])*(xx - x[1])/((x[2] - x[0])*(x[2] - x[1]))
        retval = term0*y[0] + term1*y[1] + term2*y[2]
        return retval
        
    def cumLagrange(self, which, x, y):
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
        # This approach uses Simpson's method to interpolate between beta values using the
        # interpolation polynomial in Lagrange form. This eliminates a big chunk of the bias
        # apparent with Lartillot and Phillippe method when the number of beta values used is small
        K = len(betas)
        marginal_like = 0.0
        ps_aborted = False
        for i in range(K - 2):
            x = [betas[i], betas[i+1], betas[i+2]]
            y = [means[i], means[i+1], means[i+2]]
            # print '\ni = %d (dx = %.5f, dy = %.5f)' % (i,x[0] - x[2], y[0] - y[2])
            # print 'x%d <- c(%s)' % (i,','.join(['%.6f' % z for z in x]))
            # print 'y%d <- c(%s)' % (i,','.join(['%.6f' % z for z in y]))
            # print 'plot(x%d,y%d,type="b")' % (i,i)
            if i == 0:
                a = self.cumLagrange(1, x, y)
                marginal_like += a
                b = self.cumLagrange(2, x, y)
                marginal_like += (b - a)/2.0
            elif i == K - 3:
                a = self.cumLagrange(1, x, y)
                marginal_like += a/2.0
                b = self.cumLagrange(2, x, y)
                marginal_like += (b - a)
            else:
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
        K = len(betas)
        marginal_like = 0.0
        ps_aborted = False
        for i in range(K):
            if i == 0:
                before = betas[0]
            else:
                before = betas[i-1]
            if i == K - 1:
                after = betas[K - 1]
            else:
                after = betas[i+1]
            diff = before - after
            if diff < 0.0:
                ps_aborted = True
                ps_msg = 'Phycas does not currently support path sampling from prior toward posterior'
                break
            #print 'mean %d = %f (diff = %f)' % (i, means[i], diff)
            marginal_like += means[i]*diff/2.0
        if ps_aborted:
            self.output(ps_msg)
        else:
            self.output(' %.8f Path sampling method (using trapezoid rule)' % (marginal_like))
                
    def ss(self, betas, likes):
        # Calculate marginal likelihood using steppingstone sampling method
        lnR = 0.0
        seR = 0.0
        nbetas = len(betas)
        for i in range(1,nbetas):
            beta_incr = betas[i-1] - betas[i]
            tmp = 0.0
            n = len(likes[i])
            Lmax = max(likes[i])
            for lnL in likes[i]:
                tmp += math.exp(beta_incr*(lnL - Lmax))
            tmp /= float(n)
            lnR += beta_incr*Lmax + math.log(tmp)
        
            # standard error calculation
            tmp1 = 0.0
            for lnL in likes[i]:
                aa = math.exp(beta_incr*(lnL - Lmax))/tmp
                tmp1 += math.pow((aa - 1.0),2.0)
            tmp1 /= float(n)
            seR += tmp1/float(n)
        self.output(' %.8f Steppingstone sampling method (se = %.8f)' % (lnR, seR))
        
    def marginal_likelihood(self, lines, burnin):
        beta_sum = {}
        beta_num = {}
        beta_likes = {}
        start = 2 + burnin
        for line in lines[start:]:
            parts = line.split()
            beta = float(parts[1])
            lnL  = float(parts[2])
            try:
                beta_sum[beta] += lnL
                beta_num[beta] += 1
                beta_likes[beta].append(lnL)
            except KeyError:
                beta_sum[beta] = lnL
                beta_num[beta] = 1
                beta_likes[beta] = [lnL]
        beta_vect = beta_sum.keys()
        beta_vect.sort()
        beta_vect.reverse()
        mean_vect = []
        likes_vect = []
        for b in beta_vect:
            m = beta_sum[b]/float(beta_num[b])
            mean_vect.append(m)
            likes_vect.append(beta_likes[b])
            
        self.output('\nMean log-likelihood for each value of beta used\nfor marginal likelihood estimation:\n')
        self.output('%12s %12s' % ('beta','mean lnL'))
        for b,m in zip(beta_vect,mean_vect):
            self.output('%12.5f %12.5f' % (b,m))
        self.output('\nMarginal likelihood estimates:')
        self.ss(beta_vect, likes_vect)
        self.ps_trapezoid(beta_vect, mean_vect)
        self.ps_simpsons(beta_vect, mean_vect)

    def std_summary(self, lines, burnin):
         print 'Standard summary is not yet implemented'
        
    def run(self):
        fn = self.opts.file
        burnin = self.opts.burnin
        lines = open(fn, 'r').readlines()
        if len(lines) < 3 + burnin:
            self.output("File '%s' does not look like a parameter file (too few lines)")
        else:
            headers = lines[1].split()
            if headers[1] == 'beta':
                self.marginal_likelihood(lines, burnin)
            else:
                self.std_summary(lines, burnin)

# The stuff below here can be deleted at any time...            
# def cumLagrangeWG(self, which, x, y):
#     xx = x[which]
#     
#     term1a = (xx - x[0])
#     term1b = 2.0*(x[1] - x[0])
#     term1c = (xx - x[0])*(y[1] - y[0])/term1b
#     term1 = term1a*(y[0] + term1c)
#     
#     term2a = 2.0*(xx**2.0) - xx*x[0] - (x[0]**2.0) + 3.0*x[0]*x[1] - 3.0*xx*x[1]
#     term2b = (y[2] - y[1])/(x[2] - x[1])
#     term2c = (y[1] - y[0])/(x[1] - x[0])
#     term2d = (xx - x[0])/(x[2] - x[0])
#     term2 = term2a*(term2b - term2c)*term2d/6.0
#     cum = term1 + term2
#     return -cum
#     
# def simpson_function(self, x0, x1, x2, y0, y1, y2):
#     y = (x2 - x0)*(y0 + 0.5*(x2 - x0)*(y1 - y0)/(x1 - x0)) + \
#         (2*x2**2 - x0*x2 - x0**2 + 3*x0*x1 - 3*x1*x2) * \
#         ((y2 - y1)/(x2 - x1) - (y1 - y0)/(x1 - x0))/6
#     # Add minus because we move from 1 to 0
#     return (-y)
# 
# def ps_simpsons(self, betas, U):
#     nbetavals = len(betas)
#     
#     # betas vector assumed to be sorted from high (1.0) to low (0.0)
#     # so that betas[i] - betas[i+1] > 0
#     betaIncvect = []
#     for i in range(nbetavals-1):
#         betaIncvect.append(betas[i] - betas[i+1])
# 
#     # Trapzoidal integration approximation
#     mlike = (U[0]*betaIncvect[0] + U[nbetavals - 1]*betaIncvect[nbetavals - 2])/2
#     for i in range(1, nbetavals - 1):
#         mlike += U[i]*(betaIncvect[i - 1] + betaIncvect[i])/2
# 
#     # Simpson integration approximation
#     tmp = 0.0
#     for i in range(nbetavals - 2):
#         x0 = betas[i]
#         x1 = betas[i+1]
#         x2 = betas[i+2]
# 
#         y0 = U[i]
#         y1 = U[i+1]
#         y2 = U[i+2]
# 
#         if i == 0:
#             #mu0 = -(x1 - x0)*(y0 + y1)/2.0 + (x1 - x0)**3*((y2 - y1)/(x2 - x1) - (y1 - y0)/(x1 - x0))/6
#             mu0  = -(x1 - x0)*(y0 + y1)/2.0 + (x1 - x0)**3*((y2 - y1)/(x2 - x1) - (y1 - y0)/(x1 - x0))/(6*(x2 - x0))
# 
#         if i == nbetavals - 3:
#             #mu1 = self.simpson_function(x0,x1,x2,y0,y1,y2) + \
#             #      (x1 - x0)*(y0 + y1)/2.0 - (x1 - x0)**3*((y2 - y1)/(x2 - x1) - (y1 - y0)/(x1 - x0))/6
#             mu1 = self.simpson_function(x0,x1,x2,y0,y1,y2) + \
#                    (x1 - x0)*(y0 + y1)/2.0 - (x1 - x0)**3*((y2 - y1)/(x2 - x1) - (y1 - y0)/(x1 - x0))/(6*(x2 - x0))
# 
#         tmp += self.simpson_function(x0,x1,x2,y0,y1,y2)
# 
#     tmp = (tmp + mu0 + mu1)/2.0
#     
#     self.output('\npsm_mlike = %.8f\n' % mlike)
#     self.output('\npsm_mlike_simpson = %.8f\n' % tmp)        
    
        
        