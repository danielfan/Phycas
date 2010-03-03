# Questions for Ming:
# 1. just setting gain = 1.0 seems to work; why does Faming use gain = t0/max(t,t0)?
# 2. when theta goes out of range, what do I do?

from Tkinter import *
from phycas.ProbDist import *
from tkMessageBox import askokcancel, askquestion
from tkSimpleDialog import askstring
from tkFileDialog import askopenfilename
from samcplotter import Plotter
import math

num_levels = 5
window_width = 0.3

# MyExpon is derived from ProbDist.Exponential. MyExpon gets its
# getLnPDF function from the base class, but must add getPDF
# so that SAMCDemo can plot the density function.
class MyExpon(Exponential):
    def __init__(self):
        self.hazard_rate = 2.0
        self.energy_levels = []
        self.top_heavy = True
        self.levels_on_right = True
        Exponential.__init__(self, self.hazard_rate)

    def getPDF(self, x):
        return math.exp(self.getLnPDF(x))
        
    def getName(self):
        return 'Exponential(%g)' % self.hazard_rate
        
# MyBeta is derived from ProbDist.Beta. MyBeta gets its
# getLnPDF function from the base class, but must add getPDF
# so that SAMCDemo can plot the density function.
class MyBeta(Beta):
    def __init__(self):
        self.shape = 1.5
        self.energy_levels = []
        self.top_heavy = False
        self.levels_on_right = False
        Beta.__init__(self, self.shape, self.shape)

    def getPDF(self, x):
        return math.exp(self.getLnPDF(x))
        
    def getName(self):
        return 'Beta(%g, %g)' % (self.shape, self.shape)
        
# AdHocDensity allows you to provide an arbitrary density function
# to SAMCDemo (not necessarily one derived from ProbDist.ProbabilityDistribution)
# Pass a function calculating the natural log of a probability density function to
# the constructor of AdHocDensity, then pass the AdHocDensity thus constructed to
# the SliceSampler constructor.
class AdHocDensity(AdHocDensityBase):
    def __init__(self, f):
        AdHocDensityBase.__init__(self) # very important to call base class constructor!
        self.func = f

    def __call__(self, x):
        return self.func(x)

# Here is an example of a case where we need to use AdHocDensity. MyBimodal
# implements a mixture of two Beta distributions, and thus is not represented
# by any of the distributions in ProbDist. When rolling your own probability
# distributions like this, it is important to return a very large negative number
# if the parameter is out of range. Here, the constructor queries ProbDist
# for a suitable number (returned by the getEffectiveLnZero function)
class MyBimodal:
    def __init__(self):
        self.b1 = Beta(2,19)
        self.b2 = Beta(19,2)
        self.lnzero = getEffectiveLnZero()
        self.energy_levels = []
        self.top_heavy = True
        self.levels_on_right = False

    def getLnPDF(self, x):
        if x > 0.0 and x < 1.0:
            p1 = math.exp(self.b1.getLnPDF(x))
            p2 = math.exp(self.b2.getLnPDF(x))
            p = (0.5*p1 + 0.5*p2)
            lnp = math.log(p)
        else:
            lnp = self.lnzero
        return lnp

    def getName(self):
        return 'Equal mixture of Beta(2,19) and Beta(19,2)'
        
    def getPDF(self, x):
        if x > 0.0 and x < 1.0:
            p1 = math.exp(self.b1.getLnPDF(x))
            p2 = math.exp(self.b2.getLnPDF(x))
            p = (0.5*p1 + 0.5*p2)
        else:
            p = 0.0
        return p
        
class SAMCDemo(Frame):
    def __init__(self, parent=None):
        Frame.__init__(self, parent, bg='yellow')
        self.pack(expand=YES, fill=BOTH)

        # number of points to sample when manySteps function called
        self.manyStepsN = 1000  
        
        # initialize the density from which to sample
        self.d = MyBeta()
        
        # initialize data members related to MCMC proposal
        # delta is the width of the proposal window centered
        # around current value
        self.delta = window_width
        self.x = 0.5    # will be reset in rescalePlot
        
        # initialize SAMC-related data members
        self.using_samc = True
        self.biased_pi = True
        self.m = None
        self.levels = None
        self.loglevels = None
        self.theta =None
        self.true_pi = None
        self.est_pi = None
        self.sample_size = None
        self.gain = 1.0
        self.t0 = 10.0
        self.limit_log_theta = 100.0*math.log(10.0)
        
        # initialize miscellaneous data members
        self._step_num_within_step = 0
        self._just_point = False
        self.plotter = None
        self.n_func_evals = 0
        
        # set up pseuodorandom number generator
        self.r = Lot()
        #self.r.setSeed(13579)
        print "seed is =", self.r.getSeed()
        
        # create a frame to hold the menu buttons
        menuf = Frame(self, bg='magenta')
        menuf.pack(expand=NO, fill=X)
        
        # create the File menu button
        self.filemb = Menubutton(menuf, text='File', relief=RAISED, anchor=W, borderwidth=0)
        self.filemb.pack(expand=NO, fill=X, side=LEFT)
        self.filemb.menu = Menu(self.filemb, tearoff=0)
        self.filemb['menu'] = self.filemb.menu
        self.filemb.menu.add_command(label='Quit', command=self.quit)
        
        # create the Sample menu button
        self.samplemb = Menubutton(menuf, text='Sample', relief=RAISED, anchor=W, borderwidth=0)
        self.samplemb.pack(expand=YES, fill=X, side=LEFT)
        self.samplemb.menu = Menu(self.samplemb, tearoff=1)
        self.samplemb['menu'] = self.samplemb.menu
        self.samplemb.menu.add_command(label='Toggle between SAMC and MCMC (s)', command=self.toggleSAMC)
        self.samplemb.menu.add_command(label='One step (n)', command=self.oneStep)
        self.samplemb.menu.add_command(label='Many steps (m)', command=self.manySteps)
        self.samplemb.menu.add_command(label='Many steps using slice sampler (M)', command=self.manyStepsSliceSampler)
        self.samplemb.menu.add_command(label='Biased SAMC (b)', command=self.biasedSAMC)
        self.samplemb.menu.add_command(label='Unbiased SAMC (u)', command=self.unbiasedSAMC)
        self.samplemb.menu.add_command(label='Reset (r)', command=self.reset)
        self.samplemb.menu.add_command(label='Use Beta (B)', command=self.switchToBeta)
        self.samplemb.menu.add_command(label='Use Bimodal (I)', command=self.switchToBimodal)
        self.samplemb.menu.add_command(label='Use Exponential (E)', command=self.switchToExpon)
        
        # bind some keys to the application (doesn't matter which widget has focus
        # when you use bind_all)
        self.bind_all("<KeyPress-s>", self.keybdToggleSAMC)
        self.bind_all("<KeyPress-n>", self.keybdOneStep)
        self.bind_all("<KeyPress-m>", self.keybdManySteps)
        self.bind_all("<Shift-KeyPress-M>", self.keybdManyStepsSliceSampler)
        self.bind_all("<KeyPress-b>", self.keybdBiasedSAMC)
        self.bind_all("<KeyPress-u>", self.keybdUnbiasedSAMC)
        self.bind_all("<KeyPress-r>", self.keybdReset)
        self.bind_all("<KeyPress-q>", self.keybdQuit)
        self.bind_all("<Shift-KeyPress-B>", self.keybdSwitchToBeta)
        self.bind_all("<Shift-KeyPress-I>", self.keybdSwitchToBimodal)
        self.bind_all("<Shift-KeyPress-E>", self.keybdSwitchToExpon)
                
        # create the canvas
        # func_to_plot should return density, not log density
        self.plotter = Plotter(parent=self, func_to_plot=self.d.getPDF)
        self.plotter.pack(side=TOP, expand=YES, fill=BOTH)

        # create a status bar
        self.status_label = Label(self, justify=LEFT, relief=SUNKEN, height=1, anchor=W, text='Ready')
        self.status_label.pack(side=TOP, expand=NO, fill=X)        

        self.switchToBimodal()
        
    def getLevelIndex(self, logf):
        """
        Suppose f = 0.6 (but note that log(f), not f, should be passed into this function), and 
        supposing also that these are the levels:
        level 4 (= m - 1)
        -------- 1.00
        level 3
        -------- 0.75
        level 2        <-- 0.6
        -------- 0.50
        level 1
        -------- 0.25
        level 0
        then this function would return 2.
        """
        for i in range(self.m - 1):
            if logf < self.loglevels[i]:
                return i
        return self.m - 1
        
    def takeStep(self, suppress_repaint=False):
        """
        This function proposes a new step and shows the results as a point on the plot.
        """
         # propose a new value for self.x
        u = self.r.uniform()
        x0 = (self.x - self.delta/2.0) + self.delta*u
        if x0 < 0.0:
            x0 = -x0
        logf0 = self.d.getLnPDF(x0)
        logf = self.d.getLnPDF(self.x)
        logR = logf0 - logf
        
        # adjust logR if using SAMC
        x_index  = self.getLevelIndex(logf)
        x0_index = self.getLevelIndex(logf0)
        if self.using_samc:
            logR += self.theta[x_index] - self.theta[x0_index]
            
        logu = math.log(self.r.uniform())
        if logu < logR:
            self.x = x0
            logf = logf0
            x_index = x0_index
            self.est_pi[x0_index] += 1
        else:
            self.est_pi[x_index] += 1
        self.sample_size += 1
            
        # update weights if using SAMC
        for i in range(self.m):
            if x_index == i:
                self.theta[i] += self.gain*(1.0 - self.true_pi[i])
            else:
                self.theta[i] += self.gain*(0.0 - self.true_pi[i])
        
        # point (x,y) is the point ultimately chosen along the horizontal slice
        x = self.x
        fx = self.d.getPDF(self.x)
        u = self.r.uniform()
        if self.using_samc:
            if (x_index == self.m - 1) or (self.levels[x_index] > fx):
                y_hi = fx
            else:
                y_hi = self.levels[x_index]
            if x_index > 0:
                y_lo = self.levels[x_index - 1]
            else:
                y_lo = 0.0
            y = y_lo + (y_hi - y_lo)*u
        else:
            y = fx*u
        self.plotter.points.append((x,y))
        if not suppress_repaint:            
            self.plotter.repaint()

    def _getSAMCLnPDF(self, x):
        """
        Computes the log of the weighted PDF. For use when using slice sampler to 
        sample from the SAMC target distribution.
        """
        logf     = self.d.getLnPDF(x)
        if self.using_samc:
            x_index  = self.getLevelIndex(logf)
            q        = self.theta[x_index]
        else:
            q        = 0.0
        return logf - q
        
    def takeStepSliceSampler(self, suppress_repaint=False):
        """
        This function uses the slice sampler to draw a new value and shows the 
        results as a point on the plot.
        """
        # propose a new value for self.x
        self.x   = self.s.sample()
        logf     = self.d.getLnPDF(self.x)
        
        # update weights and counts if using SAMC
        if self.using_samc:
            x_index  = self.getLevelIndex(logf)
    
            # update counts
            self.est_pi[x_index] += 1
            self.sample_size += 1
            
            # update weights            
            for i in range(self.m):
                if x_index == i:
                    self.theta[i] += self.gain*(1.0 - self.true_pi[i])
                else:
                    self.theta[i] += self.gain*(0.0 - self.true_pi[i])
        
        # point (x,y) is the point ultimately chosen along the horizontal slice
        fx = self.d.getPDF(self.x)
        u = self.r.uniform()
        if self.using_samc:
            if (x_index == self.m - 1) or (self.levels[x_index] > fx):
                y_hi = fx
            else:
                y_hi = self.levels[x_index]
            if x_index > 0:
                y_lo = self.levels[x_index - 1]
            else:
                y_lo = 0.0
            y = y_lo + (y_hi - y_lo)*u
        else:
            y = fx*u
        self.plotter.points.append((self.x,y))
        if not suppress_repaint:            
            self.plotter.repaint()

    def toggleSAMC(self):
        if self.using_samc:
            self.modifyStatus('Switched to MCMC')
            self.using_samc = False
        else:
            self.modifyStatus('Switched to SAMC')
            self.using_samc = True
        self.reset()
        
    def keybdToggleSAMC(self, event):
        self.toggleSAMC()

    def oneStep(self):
        self.takeStep()
        self.modifyStatus('Took 1 step (%d total steps so far)' % self.sample_size)
        
    def keybdOneStep(self, event):
        self.oneStep()

    def biasedSAMC(self):
        self.modifyStatus('Lowest energy level will now be visited twice as often as other levels')
        self.using_samc = True
        self.biased_pi = True
        self.reset()
        
    def keybdBiasedSAMC(self, event):
        self.biasedSAMC()

    def unbiasedSAMC(self):
        self.modifyStatus('Each energy level will now be visited equally')
        self.using_samc = True
        self.biased_pi = False
        self.reset()
        
    def keybdUnbiasedSAMC(self, event):
        self.unbiasedSAMC()

    def switchToBeta(self):
        self.modifyStatus('Switching to Beta(1.5) distribution')
        self.d = MyBeta()
        self.reset()

    def switchToExpon(self):
        self.modifyStatus('Switching to Exponential(2.0) distribution')
        self.d = MyExpon()
        self.reset()
        
    def switchToBimodal(self):
        self.modifyStatus('Switching to a distribution that is an equal mixture of Beta(2,19) and Beta(19,2)')
        self.d = MyBimodal()
        self.reset()
       
    def keybdSwitchToBeta(self, event):
        self.switchToBeta()

    def keybdSwitchToBimodal(self, event):
        self.switchToBimodal()

    def keybdSwitchToExpon(self, event):
        self.switchToExpon()

    def manySteps(self):
        self.modifyStatus('Patience...')
        for x in range(self.manyStepsN):
            self.takeStep(suppress_repaint=True)
        self.plotter.repaint()
        self.modifyStatus('Took %d steps (%d total steps so far)' % (self.manyStepsN,self.sample_size))
            
    def keybdManySteps(self, event):
        self.manySteps()

    def manyStepsSliceSampler(self):
        self.modifyStatus('Patience...')
        for x in range(self.manyStepsN):
            self.takeStepSliceSampler(suppress_repaint=True)
        self.plotter.repaint()
        self.modifyStatus('Took %d steps using slice sampler (%d total steps so far)' % (self.manyStepsN,self.sample_size))
            
    def keybdManyStepsSliceSampler(self, event):
        self.manyStepsSliceSampler()

    #def reset(self):
    #    self.x = starting_x
    #    self.modifyStatus('Ready')
    #    self.initLevels()
    #    if self.plotter:
    #        try:
    #            pw, ph, plotm = self.plotter.plotw, self.plotter.ploth, self.plotter.plotm
    #        except:
    #            pw, ph, plotm = 800, 600, 40
    #        self.plotter.resize(pw, ph, plotm)
    #        self.plotter.reset()
                
    def keybdReset(self, event):
        self.reset()
        
    # next two functions cause the application to terminate
    def quit(self):
        Frame.quit(self)

    def keybdQuit(self, event):
        self.quit()
        
    def modifyStatus(self, msg):
        self.status_label.config(text=msg)
        self.status_label.update_idletasks()
        
    def initLevels(self, ymin, ymax, nlevels):
        """
        If self.biased_pi is True, lowest energy level is made twice as probable
        as all other energy levels. If self.biased_pi is False, all energy levels
        have the same probability.
        """
        # set up energy levels for this density
        # boundary between level 0 and 1 is log(ymin)
        # other boundaries divided evenly between log(ymin) and log(ymax)
        log_ymax = math.log(ymax)
        log_ymin = math.log(ymin)
        yincr = (log_ymax - log_ymin)/float(num_levels - 1)
        self.d.energy_levels = []
        for i in range(num_levels - 1):
            self.d.energy_levels.append(log_ymin + yincr*float(i))

        self.loglevels = self.d.energy_levels
        self.m = len(self.loglevels) + 1
        self.levels = [math.exp(v) for v in self.loglevels]
        self.theta = [0.0]*self.m
        if self.biased_pi:
            self.true_pi = [1.0]*self.m
            self.true_pi[0] = 2.0
            denom = sum(self.true_pi)
            self.true_pi = [p/denom for p in self.true_pi]
        else:
            self.true_pi = [1.0/float(self.m)]*self.m
        self.est_pi = [0]*self.m
        self.sample_size = 0

    def reset(self):
        """
        Called after switching the probabilty density to explore.
        """
        # self.plotter.func should be function that returns a probability density,
        # not a log density
        self.plotter.func = self.d.getPDF

        # self.f and self.ff should be functions that return the natural log of a 
        # probability density, or the natural log of a function proportional to a 
        # probability density
        #self.f = AdHocDensity(self.d.getLnPDF)
        self.f = AdHocDensity(self._getSAMCLnPDF)
        
        # create a slice sampler
        self.s = SliceSampler(self.r, self.f)
        
        # create the status label
        self.modifyStatus('Ready')
        
        # grab a few samples to use in scaling the plot (number to generate is passed
        # to the rescalePlot function)
        self.rescalePlot(1000)

        # reset statistics such as average slice width, number of function evalutations, etc.
        w = self.s.getSliceUnitWidth()
        self.s.setSliceUnitWidth(w/3.0)
        self.s.resetDiagnostics()
                
    def rescalePlot(self, n):
        """
        The rescalePlot function draws n samples from the target distribution to determine
        the approximate limits of the main part of the density function. It then sets
        the plotter's x and y axis ranges accordingly 
        """
        # temporarily turn off samc for purposes of rescaling the plot
        prev_using_samc = self.using_samc
        self.using_samc = False
        
        smallest = getEffectiveLnZero()
        largest = -smallest
        ymin = 0.0              # always zero
        ymin_not_zero = None    # smallest y value seen that is greater than zero
        ymax = smallest
        xmin = largest
        xmax = smallest
        for i in range(n):
            x = self.s.sample()
            if x < xmin:
                xmin = x
            if x > xmax:
                xmax = x
            y = math.exp(self.f(x))
            if y > ymax:
                ymax = y
            if (ymin_not_zero is None) or (y > 0.0 and y < ymin_not_zero):
                ymin_not_zero = y
                
        # restore using_samc setting
        self.using_samc = prev_using_samc
                
        # start x off at a random point between xmin and xmax
        self.x = xmin + (xmax - xmin)*self.r.uniform()
              
        # store the energy levels for this distribution
        self.initLevels(ymin_not_zero, ymax, num_levels)
                
        # make sure plot is wide enough to accommodate the worst-case scenario
        # in terms of positioning the initial interval of width w
        w = self.s.getSliceUnitWidth()
        self.plotter.setXRange(xmin - w, xmax + w, (xmax - xmin + 2*w)/5.0)            
        self.plotter.setYRange(ymin, ymax, (ymax - ymin)/5.0)

        try:
           pw, ph, plotm = self.plotter.plotw, self.plotter.ploth, self.plotter.plotm
        except:
           pw, ph, plotm = 800, 600, 40
        self.plotter.resize(pw, ph, plotm)
        self.plotter.reset()
        
if __name__ == '__main__':
    app = Tk()
    app.geometry("%dx%d%+d%+d" % (800, 600, 0, 0))
    SAMCDemo(app).mainloop()
