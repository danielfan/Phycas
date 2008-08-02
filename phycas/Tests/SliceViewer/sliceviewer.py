from Tkinter import *
from phycas.ProbDist import *
from tkMessageBox import askokcancel, askquestion
from tkSimpleDialog import askstring
from tkFileDialog import askopenfilename
from sliceplotter import Plotter
import math

# MyExpon is derived from ProbDist.Exponential. MyExpon gets its
# getLnPDF function from the base class, but must add getPDF
# so that SliceViewer can plot the density function.
class MyExpon(Exponential):
    def __init__(self):
        Exponential.__init__(self, 2.0)

    def getPDF(self, x):
        return math.exp(self.getLnPDF(x))

# MyBeta is derived from ProbDist.Beta. MyBeta gets its
# getLnPDF function from the base class, but must add getPDF
# so that SliceViewer can plot the density function.
class MyBeta(Beta):
    def __init__(self, strength):
        Beta.__init__(self, strength, strength)

    def getPDF(self, x):
        return math.exp(self.getLnPDF(x))

# AdHocDensity allows you to provide an arbitrary density function
# to SliceViewer (not necessarily one derived from ProbDist.ProbabilityDistribution)
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

    def getLnPDF(self, x):
        if x > 0.0 and x < 1.0:
            p1 = math.exp(self.b1.getLnPDF(x))
            p2 = math.exp(self.b2.getLnPDF(x))
            p = (0.5*p1 + 0.5*p2)
            lnp = math.log(p)
        else:
            lnp = self.lnzero
        return lnp

    def getPDF(self, x):
        if x > 0.0 and x < 1.0:
            p1 = math.exp(self.b1.getLnPDF(x))
            p2 = math.exp(self.b2.getLnPDF(x))
            p = (0.5*p1 + 0.5*p2)
        else:
            p = 0.0
        return p

# Results of applying different adapt modes to MyBimodal
# (results in each case are the mean number of extra function evaluations
# based on a sample size of 10,000)
#
#  Default |              Simple               |                Neal                |      y-conditional (from end: 0.25)
#    NA    | 0.125 0.25  0.5   1.0   2.0   4.0 | 0.25  0.5   1.0   2.0   4.0   8.0  | 0.25  0.5   0.8   1.0   2.0   2.2   4.0 |  
#  --------|-----------------------------------|------------------------------------|-----------------------------------------
#    2.9   | 2.5  [1.9]  2.0   2.6   3.2   3.9 | 6.4   3.6   2.2  [1.9]  2.1   2.8  | 4.4   2.5   2.0  [1.8] [1.8]  1.9   2.4


class SliceViewer(Frame):
    def __init__(self, parent=None):
        Frame.__init__(self, parent, bg='yellow')
        self.pack(expand=YES, fill=BOTH)

        self.adaptMode = "simple" # 'yconditional'   # choices are: simple, neal, or yconditional
        self.adaptSimpleParam = 0.25
        self.adaptNealParam = 2.0
        self.adaptYCondParam = 1.3
        self.adaptYCondFromEnds = 0.25
        
        self.manyStepsN = 1000  # number of points to sample when manySteps function called
        #self.total_steps = 0
        
        self.n_func_evals = 0
        
        self.r = Lot()
        #self.r.setSeed(13579)
        print "seed is =", self.r.getSeed()
        #self.d = MyBeta(1.5)
        # self.d = MyExpon()
        self.d = MyBimodal() 
        self.f = AdHocDensity(self.d.getLnPDF)
        
        # self.f should be function that returns natural log of a probability density, or
        # the natural log of a function proportional to a probability density
        self.s = SliceSampler(self.r, self.f)

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
        self.samplemb.menu.add_command(label='One step', command=self.oneStep)
        self.samplemb.menu.add_command(label='Many steps', command=self.manySteps)
        self.samplemb.menu.add_command(label='Adapt', command=self.adaptSampler)
        self.samplemb.menu.add_command(label='Reset', command=self.reset)

        # create the canvas
        # func_to_plot should return density, not log density
        self.plotter = Plotter(parent=self, func_to_plot=self.d.getPDF)
        self.plotter.pack(side=TOP, expand=YES, fill=BOTH)
        
        # create the status label
        self.status_label = Label(self, justify=LEFT, relief=SUNKEN, height=1, anchor=W, text='Ready')
        self.status_label.pack(side=TOP, expand=NO, fill=X)
        
        # grab a few samples to use in scaling the plot (number to generate is passed
        # to the rescalePlot function)
        self.rescalePlot(1000)

        # reset statistics such as average slice width, number of function evalutations, etc.
        w = self.s.getSliceUnitWidth()
        self.s.setSliceUnitWidth(w/2.0)
        self.s.resetDiagnostics()

        # bind some keys to the application (doesn't matter which widget has focus
        # when you use bind_all
        self.bind_all("<KeyPress-o>", self.keybdOverrelaxedOneStep)
        self.bind_all("<KeyPress-n>", self.keybdOneStep)
        self.bind_all("<KeyPress-d>", self.keybdDetailedStep)
        self.bind_all("<KeyPress-b>", self.keybdBackDetailedStep)
        self.bind_all("<Shift-KeyPress-N>", self.keybdManySteps)
        self.bind_all("<KeyPress-m>", self.keybdManySteps)
        self.bind_all("<KeyPress-a>", self.keybdAdaptSampler)
        self.bind_all("<KeyPress-r>", self.keybdReset)
        self.bind_all("<KeyPress-q>", self.keybdQuit)

        self._step_num_within_step = 0
        self._just_point = False
        # configure event is bound only to the main frame
        #self.bind("<Configure>", self.resizing)

    # The rescalePlot function draws n samples from the target distribution to determine
    # the approximate limits of the main part of the density function. It then sets
    # the plotter's x and y axis ranges accordingly 
    def rescalePlot(self, n):
        smallest = getEffectiveLnZero()
        largest = -smallest
        ymin = 0.0 # always zero
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
        # make sure plot is wide enough to accommodate the worst-case scenario
        # in terms of positioning the initial interval of width w
        w = self.s.getSliceUnitWidth()
        self.plotter.setXRange(xmin - w, xmax + w, (xmax - xmin + 2*w)/5.0)            
        self.plotter.setYRange(ymin, ymax, (ymax - ymin)/5.0)
        print 'range of x is [%f, %f]' % (xmin, xmax)
        print 'range of y is [%f, %f]' % (ymin, ymax)

    def takeStep(self, show_slice=False, detailed=False):
        """
        If detailed is False, this function displays an entire slice sampling step.
        If detailed is True, this function shows the details of a single slice sampling step.
        That is, it first shows the vertical slice, then the point chosen for the horizontal
        slice, then the horizontal slice, etc.
        """
        target_step = 1
        if (not detailed) or self._step_num_within_step == 0:
            #self.total_steps += 1
            # Tell the slice sampler to take one step and return all the information
            # about that step
            self.curr = self.s.debugSample()
        if not detailed:
            self._step_num_within_step = 1000000
            
        # debugSample returns a tuple with the following information about the 
        # slice sampling step just completed:
        #   0: sampled x
        #   1: x-coord of vertical slice
        #   2: y-coord of top of vertical slice (y-coord of bottom of vertical slice always 0.0)
        #   3: x-coord of left edge of horizontal slice
        #   4: x-coord of right edge of horizontal slice
        #   5: y-coord of horizontal slice
        #   6: horizontal slice interval width
        #   7+: x-coord of failed sampling attempts (y-coord all equal to element 5)
        
        # point (x0,y0) is the top of the vertical slice
        x0 = self.curr[1]
        y0 = math.exp(self.curr[2])
        
        # point (x,y) is the point ultimately chosen along the horizontal slice
        x = self.curr[0]
        y = math.exp(self.curr[5])
        
        # xleft and xright are the x-coordinates of the left and right ends of the 
        # horizontal slice
        xleft = self.s.getOrigLeftEdgeOfSlice()
        xright = self.s.getOrigRightEdgeOfSlice()
        
        #slice_extent = float(self.curr[4] - self.curr[3])
        #outf = file('yw.txt', 'a+')
        #outf.write('%d\t%f\t%f\n' % (self.total_steps, y, slice_extent))
        #outf.close()                    
        
        # xincr is the slice unit width
        xincr = self.curr[6]
        
        # curr_len is the length of the list returned by the slice sampler's debugSample method
        # This information helps us tell how many failed attempts there were
        curr_len = len(self.curr)
        
        # Determine the number of horiz. slice units needed
        unit_width = self.s.getSliceUnitWidth()
        assert xincr == unit_width
        n_units = int(0.5 + (xright - xleft)/xincr)
        
        # Get number of function evaluations needed for the current slice sample
        # At minimum, four function evaluations are required:
        #   1. re-evaluate current point (density has probably changed since last update)
        #   2. check left edge to see if it is outside slice
        #   3. check right edge to see if it is outside slice
        #   4. check sampled point to see if it is valid
        prev_n_func_evals = self.n_func_evals
        self.n_func_evals = self.s.getNumFuncEvals()
        extra_evals = self.n_func_evals - prev_n_func_evals - 4
        n_samples = self.s.getNumSamples()
        efficiency = (float(self.n_func_evals)/float(n_samples)) - 4.0
        #self.status_label.config(text='unit width=%.2f, units=%d, extra evals=%d, extra evals/sample=%.1f' % (unit_width, n_units, extra_evals, efficiency))
        if show_slice:
            self.plotter.repaint()
            if not self._just_point:
                self.plotter.plotLine(x0, 0.0, x0, y0)
                self._step_num_within_step += 1
                
                # Plot the vertical slice
                self.status_label.config(text='vertical slice')
                if self._step_num_within_step == target_step:
                    return
                target_step += 1
    
                n_failed = curr_len - 7
                #print 'n_failed=%d' % (n_failed)
                #print 'self.s.getNumSamples()=', self.s.getNumSamples()
                #print 'self.s.getNumFailedSamples()=', self.s.getNumFailedSamples()
                #print 'self.s.getNumFuncEvals()=', self.s.getNumFuncEvals()
                if self.adaptMode == 'yconditional':
                    mode = self.s.getMode()
                    mode_height = math.exp(self.s.getLnDensityAtMode())
                    w0 = self.s.calcW(mode_height)
                    w1 = self.s.calcW(0.0)
                    left_x0 = mode - w0/2.0
                    left_y0 = mode_height
                    left_x1 = mode - w1/2.0
                    left_y1 = 0.0
                    right_x0 = mode + w0/2.0
                    right_y0 = mode_height
                    right_x1 = mode + w1/2.0
                    right_y1 = 0.0
                    self.plotter.plotLine(left_x0, left_y0, left_x1, left_y1, "green")
                    self.plotter.plotLine(right_x0, right_y0, right_x1, right_y1, "green")
                
                # Plot an arrow indicating where the horizontal slice will be located
                if self._step_num_within_step == target_step:
                    self.plotter.plotLeftPointingArrow(x0, y, color="white", shaft_length=10, arrowhead_length=5)
                    self.status_label.config(text='location of horizontal slice')
                    return
                target_step += 1
                
                # Plot the unit that spans the vertical slice
                n_increases = int(0.5 + (xright - xleft)/xincr) # int((xright - xleft)//xincr)
        
                growingSlice = self._step_num_within_step < target_step + n_increases
                if growingSlice:
                    nunits, n_increases, first_unit_index = self.plotter.plotOneSliceUnit(xleft, xright, xincr, x0, y, which=None, showTicks=True)

                    assert nunits == n_units, 'nunits = %d, n_units = %d' %(nunits, n_units)
                    self.status_label.config(text='first unit randomly positioned around vertical slice')
                    if self._step_num_within_step == target_step:
                        return
                    target_step += 1
                    
                    # Plot the units to the right of the first unit
                    curr_index = first_unit_index + 1 
                    while curr_index  < n_increases :
                        self.plotter.plotOneSliceUnit(xleft, xright, xincr, x0, y, curr_index) 
                        if self._step_num_within_step == target_step:
                            self.status_label.config(text='add units to right until density bracketed')
                            return
                        target_step += 1
                        curr_index += 1
                        
                    # Plot the units to the left of the first unit
                    curr_index = first_unit_index - 1 
                    while curr_index >= 0:
                        self.plotter.plotOneSliceUnit(xleft, xright, xincr, x0, y, curr_index)                
                        if self._step_num_within_step == target_step:
                            self.status_label.config(text='add units to left until density bracketed')
                            return
                        target_step += 1
                        curr_index -= 1
                else:
                    target_step += n_increases
                left_border, right_border = self.plotter.getExpandedSliceDim(xleft, xright, xincr)
                if self._step_num_within_step == target_step:
                    self.plotter.plotCurrentSliceWidth(left_border, right_border, x0, y)
                    self.status_label.config(text='Now we have the full slice that we will sample from')
                    return
                target_step += 1
                if n_failed > 0:
                    for i in range(n_failed):
                        for j in (0,1):
                            xfail = self.curr[i + 7]
                            xft = self.plotter.xtranslate(xfail)
                            if j == 1:
                                if x0 < xfail:
                                    right_border = xft
                                else:
                                    left_border = xft
                            if not self._just_point:
                                self.plotter.plotPoint(xfail, y, 'red')                                
                            if self._step_num_within_step == target_step:
                                self.plotter.plotCurrentSliceWidth(left_border, right_border, x0, y)
                                if j == 0:
                                    self.status_label.config(text='sample attempt failed')
                                elif j == 1:
                                    self.status_label.config(text='We can crop the edges of the slice')                                
                                return
                            target_step += 1
            
            # this is the last thing to plot, so we reset to _step_num_within_step to 0
            self.plotter.plotPoint(x, y, 'cyan')
            if self._just_point:
                self._step_num_within_step = 0
                self._just_point = False
                self.plotter.points.append((x,y))
            else:
                self.plotter.plotCurrentSliceWidth(left_border, right_border, x0, y)
                self._just_point = True
            self.status_label.config(text='successful sample at x=%.3f' % x)
            #print 'xleft=%f, xright=%f, xincr=%f\n' % (xleft, xright, xincr)
        else:
            self.plotter.points.append((x,y))

    def takeOverrelaxedStep(self, show_slice=False):
        self.curr = self.s.debugOverrelaxedSample()
        # 0: sampled x
        # 1: x-coord of vertical slice
        # 2: y-coord of top of vertical slice (y-coord of bottom of vertical slice always 0.0)
        # 3: x-coord of left edge of horizontal slice
        # 4: x-coord of right edge of horizontal slice
        # 5: y-coord of horizontal slice
        x = self.curr[0]
        y = math.exp(self.curr[5])
        x0 = self.curr[1]
        y0 = math.exp(self.curr[2])
        self.plotter.points.append((x,y))
        if show_slice:
            self.plotter.repaint()
            xleft = self.curr[3] # self.s.getLeftEdgeOfSlice()
            xright = self.curr[4] # self.s.getRightEdgeOfSlice()
            xincr = xright - xleft
            self.plotter.plotSlice(xleft, xright, xincr, y)
            self.plotter.plotLine(x0, 0.0, x0, y0)
            self.plotter.plotPoint(x, y, 'cyan')

    #def resizing(self, event):
    #    print "resizing: x=%f, y=%f" % (event.width, event.height)

    # next two functions result in a single sample and show details of the slice sampling process
    def oneStep(self, detailed=False):
        self.takeStep(show_slice=True, detailed=detailed)
        
    def keybdOneStep(self, event):
        self.oneStep(detailed=False)

    def keybdDetailedStep(self, event):
        self.oneStep(detailed=True)

    def keybdBackDetailedStep(self, event):
        if self._step_num_within_step < 1:
            return
        self._step_num_within_step = max(1, self._step_num_within_step - 2)
        self.oneStep(detailed=True)

    # next two functions result in a single overrelaxed sample and show details of the slice sampling process
    def overrelaxedOneStep(self):
        self.takeOverrelaxedStep(True)
        
    def keybdOverrelaxedOneStep(self, event):
        self.overrelaxedOneStep()
        
    # next two functions result in many samples but do not show details of the slice sampling process
    def manySteps(self):
        for x in range(self.manyStepsN):
            self.takeStep()
        self.plotter.repaint()
            
    def keybdManySteps(self, event):
        self.manySteps()

    # next two functions clear the canvas and start the sampling over from scratch        
    def reset(self):
        self.plotter.reset()
                
    def keybdReset(self, event):
        self.reset()
        
    # next two functions adjust the slide unit width based on previous slice widths
    def adaptSampler(self):
        if self.adaptMode == 'simple':
            # average cropped slice width is divided by self.adaptSimpleParam to get new unit width
            w = self.s.adaptSimple(self.adaptSimpleParam)
            self.status_label.config(text='Simple adapt mode: multiplier = %f, new unit width = %f' % (self.adaptSimpleParam, w))
        elif self.adaptMode == 'neal':
            # average distance between successive sampled points times self.adaptNealParam is used as new unit width
            w = self.s.adaptNeal(self.adaptNealParam)
            self.status_label.config(text='Neal adapt mode: multiplier = %f, new unit width = %f' % (self.adaptNealParam, w))
        elif self.adaptMode == 'yconditional':
            # new unit width depends on how cropped slice interval varies with y (the value chosen along
            # interval from 0 to f(x) that defines the vertical height of the slice). If this mode is in
            # force, w is changed every time a sample is drawn
            self.s.adaptYConditional(self.adaptYCondFromEnds, self.adaptYCondParam)
            self.status_label.config(text='y-conditional adapt mode: from_ends = %f, multiplier = %f' % (self.adaptYCondFromEnds, self.adaptYCondParam))
        self.s.resetDiagnostics()

    def keybdAdaptSampler(self, event):
        self.adaptSampler()

    # next two functions cause the application to terminate
    def quit(self):
        Frame.quit(self)

    def keybdQuit(self, event):
        self.quit()
        
if __name__ == '__main__':
    app = Tk()
    app.geometry("%dx%d%+d%+d" % (800, 600, 0, 0))
    SliceViewer(app).mainloop()
