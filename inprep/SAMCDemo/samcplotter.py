from Tkinter import *
import tkFont

LINE_WIDTH = 3
class Plotter(Canvas):
    def __init__(self, parent, func_to_plot):
        self.background_color = 'blue'
        self.axisfont_color = 'white'
        self.axis_color = 'cyan'
        self.slice_color = 'yellow'
        self.func_color = 'white'
        self.func = func_to_plot
        Canvas.__init__(self, master=parent, bg=self.background_color)

        self.font = tkFont.Font(family='Courier', size=12)
        fm = self.font.metrics()
        self.font_height = fm['ascent'] + fm['descent']

        self.points = []

        self.setXRange(0.0, 1.0, 0.2)
        self.setYRange(0.0, 2.0, 0.5)
        
        self.bind("<Configure>", self.resizeEvent)
        
        self.demo = parent
        
        self.analtype_pos = None   # will be set when plot is rescaled (see resize function)
        self.analtype_font = tkFont.Font(family='Arial', size=24)
        self.analtypemodifier_font = tkFont.Font(family='Arial', size=18)
        fm = self.font.metrics()
        self.analtype_font_height = fm['ascent'] + fm['descent']
        
        self.distname_pos = None   # will be set when plot is rescaled (see resize function)
        self.distname_font = tkFont.Font(family='Arial', size=12)

    def resizeEvent(self, event):
        self.resize(event.width, event.height, 40)
        self.repaint()

    def resize(self, new_width, new_height, new_margin):
        self.plotw = new_width
        self.ploth = new_height
        self.plotm = new_margin
        assert self.font_height < new_margin, 'height of font too large for specified plot margin'
        self.offset = new_margin/2
        self.left = self.plotm
        self.right = self.plotw - self.plotm
        self.bottom = self.ploth - self.plotm
        self.top = self.plotm
        self.xaxislen = self.plotw - 2*self.plotm
        self.yaxislen = self.ploth - 2*self.plotm
        self.analtype_pos = (self.right - 4*self.offset, 0.5*self.top)
        self.distname_pos = ((self.right + self.left)/2.0, 0.5*self.top)
        #Canvas.config(self, width=new_width, height=new_height)        
        
    def repaint(self):
        Canvas.create_rectangle(self, 0, 0, self.plotw, self.ploth, fill=self.background_color);
        self.drawAxes()
        self.drawLevels()
        self.plotFunc(self.func, 100)
        self.plotAllPoints()

    def reset(self):
        self.points = []
        import gc
        gc.collect()
        self.repaint()
        
    def setXRange(self, lowerbound, upperbound, incr):
        self.xlow = lowerbound
        self.xhigh = upperbound
        self.xincr = incr
        self.xticks = 1 + int((upperbound - lowerbound)/incr)
        
    def setYRange(self, lowerbound, upperbound, incr):
        self.ylow = lowerbound
        self.yhigh = upperbound
        self.yincr = incr
        self.yticks = 1 + int((upperbound - lowerbound)/incr)
        
    def drawAxes(self):
        Canvas.create_line(self, self.left, self.bottom, self.right, self.bottom, width=LINE_WIDTH, fill=self.axis_color)
        Canvas.create_line(self, self.left, self.bottom, self.left, self.top, width=LINE_WIDTH, fill=self.axis_color)
        for i in range(self.xticks):
            v = self.xlow + float(i)*self.xincr
            x = self.left + self.xaxislen*((v - self.xlow)/(self.xhigh - self.xlow))
            y = self.bottom + self.offset
            Canvas.create_text(self, x, y, text=str('%.1f' % v), font=self.font, fill=self.axisfont_color)
        for i in range(self.yticks):
            v = self.ylow + float(i)*self.yincr
            x = self.left - self.offset
            y = self.bottom - self.yaxislen*((v - self.ylow)/(self.yhigh - self.ylow))
            Canvas.create_text(self, x, y, text=str('%.1f' % v), font=self.font, fill=self.axisfont_color)

    def drawLevels(self):
        if self.analtype_pos is not None:
            if self.demo.using_samc:
                anal_type = 'SAMC'
                if self.demo.biased_pi:
                    anal_modifier = '(biased)'
                else:
                    anal_modifier = '(unbiased)'
            else:
                anal_type = 'MCMC'
                anal_modifier = ''
            Canvas.create_text(self, self.analtype_pos[0], self.analtype_pos[1], text=str(anal_type), justify=CENTER, font=self.analtype_font, fill=self.axisfont_color)
            Canvas.create_text(self, self.analtype_pos[0], self.analtype_pos[1] + 1.5*self.analtype_font_height, text=str(anal_modifier), justify=CENTER, font=self.analtypemodifier_font, fill=self.axisfont_color)
            dist_name = self.demo.d.getName()
            Canvas.create_text(self, self.distname_pos[0], self.distname_pos[1], text=str(dist_name), anchor=CENTER, font=self.distname_font, fill=self.axisfont_color)
        if self.demo.using_samc and self.demo.m is not None:
            prev_y_lo = self.bottom
            m = self.demo.m
            sz = self.demo.sample_size
            if sz == 0:
                sz = 1
            # For distributions labeled as top_heavy, the level lines are scrunched together closely at the bottom
            # and leave no room to show thetas and counts. For these distributions, put all theta and count reports
            # in the top part of the plot area
            # m = 4
            # i = 3  top - 1*font_height = top - (m - 3)*font_height
            # i = 2  top - 2*font_height = top - (m - 2)*font_height
            # i = 1  top - 3*font_height = top - (m - 1)*font_height
            # i = 0  top - 4*font_height = top - (m - 0)*font_height 
            for i in range(m - 1):
                q = self.demo.theta[i]
                p = self.demo.est_pi[i]
                f = self.demo.levels[i]
                y_hi = self.bottom - self.yaxislen*((f - self.ylow)/(self.yhigh - self.ylow))
                if self.demo.d.top_heavy:
                    y_avg = self.top + self.font_height + 2.0*self.font_height*float(m - i)
                else:
                    y_lo = prev_y_lo
                    y_avg = (y_lo + y_hi)/2.0
                Canvas.create_line(self, self.left, y_hi, self.right, y_hi, width=LINE_WIDTH, fill='white')
                if self.demo.d.levels_on_right:
                    Canvas.create_text(self, self.right - 2*self.offset, y_avg, text=str('level %d, theta = %.3f, count = %d' % (i,q,p)), anchor=E, font=self.font, fill=self.axisfont_color)
                else:
                    Canvas.create_text(self, self.left + self.offset, y_avg, text=str('level %d, theta = %.3f' % (i,q)), anchor=W, font=self.font, fill=self.axisfont_color)
                    Canvas.create_text(self, self.right - 2*self.offset, y_avg, text=str('count = %d' % p), anchor=E, font=self.font, fill=self.axisfont_color)
                prev_y_lo = y_hi
            i = m - 1
            q = self.demo.theta[i]
            p = self.demo.est_pi[i]
            y_hi = self.top
            if self.demo.d.top_heavy:
                y_avg = self.top + self.font_height + 2.0*self.font_height
            else:
                y_lo = prev_y_lo
                y_avg = (y_lo + y_hi)/2.0
            if self.demo.d.levels_on_right:
                Canvas.create_text(self, self.right - 2*self.offset, y_avg, text=str('level %d, theta = %.3f, count = %d' % (i,q,p)), anchor=E, font=self.font, fill=self.axisfont_color)
            else:
                Canvas.create_text(self, self.left + self.offset, y_avg, text=str('level %d, theta = %.3f' % (i,q)), anchor=W, font=self.font, fill=self.axisfont_color)
                Canvas.create_text(self, self.right - 2*self.offset, y_avg, text=str('count = %d' % p), anchor=E, font=self.font, fill=self.axisfont_color)

    def xtranslate(self, x):
        return int(self.left + self.xaxislen*((x - self.xlow)/(self.xhigh - self.xlow)))

    def ytranslate(self, y):
        return int(self.bottom - self.yaxislen*((y - self.ylow)/(self.yhigh - self.ylow)))

    def plotCurrentSliceWidth(self, xleft, xright, xvert, y):
        x0 = self.xtranslate(xleft)
        x1 = self.xtranslate(xright)
        xv = self.xtranslate(xvert)
        y = self.ytranslate(y)
        Canvas.create_line(self, xleft, y, xright, y, fill=self.slice_color, width=LINE_WIDTH)

    def getExpandedSliceDim(self, xleft, xright, xincr):
        n = int(0.5 + (xright - xleft)/xincr) # int((xright - xleft)//xincr)
        ticks = [xleft + xincr*i for i in range(n + 1)]
        tf = ticks[0]
        tl = ticks[-1]
        return self.xtranslate(tf), self.xtranslate(tl)

    def plotOneSliceUnit(self, xleft, xright, xincr, xvert, y, which = None, showTicks=True):
        """
        Plots only one slice unit. If which is None, the unit plotted is the
        one that spans the vertical slice (at x-coordinate xvert). If which
        is an integer, plots the unit whose index is that integer. The return 
        value is the tuple (n,k), where n is the number of slice units
        and k the index of the unit plotted. 
        """
        x0 = self.xtranslate(xleft)
        x1 = self.xtranslate(xright)
        xv = self.xtranslate(xvert)
        y = self.ytranslate(y)
        n = int(0.5 + (xright - xleft)/xincr) # int((xright - xleft)//xincr)
        ticks = [xleft + xincr*i for i in range(n + 1)]
        last_incr = len(ticks) - 1
        for k, tick in enumerate(ticks):
            xleft = self.xtranslate(tick)
            xright = self.xtranslate(tick + xincr)
            
            # If no index was supplied, check to see if the current unit spans
            # the vertical slice. If so, arrange for this unit to be plotted
            if (xleft < xv) and (xright > xv) and (which is None):
                which = k
                
            if which == k:
                # Draw the horizontal line
                Canvas.create_line(self, xleft, y, xright, y, fill=self.slice_color, width=LINE_WIDTH)
                if showTicks or k == 0:
                    # Draw tick mark at left end
                    Canvas.create_line(self, xleft, y-2, xleft, y+2, fill=self.slice_color, width=LINE_WIDTH)
                if showTicks or k == last_incr:
                    # Draw tick mark at right end
                    Canvas.create_line(self, xright, y-2, xright, y+2, fill=self.slice_color, width=LINE_WIDTH)
        assert which is not None, " ".join([str(i) for i in [x0, xv, x1, xleft, xright, xincr, xvert, y]])
        return n, k, which
        
    def plotSlice(self, xleft, xright, xincr, y):
        x0 = self.xtranslate(xleft)
        x1 = self.xtranslate(xright)
        y = self.ytranslate(y)
        Canvas.create_line(self, x0, y, x1, y, fill=self.slice_color, width=LINE_WIDTH)
        n = int((xright - xleft)//xincr)
        ticks = [xleft + xincr*i for i in range(n + 1)]
        for tick in ticks:
            x = self.xtranslate(tick)
            Canvas.create_line(self, x, y-2, x, y+2, fill=self.slice_color, width=LINE_WIDTH)
        
    def plotLine(self, x0, y0, x1, y1, line_color=None):
        ix0 = self.xtranslate(x0)
        ix1 = self.xtranslate(x1)
        iy0 = self.ytranslate(y0)
        iy1 = self.ytranslate(y1)
        if line_color==None:
            Canvas.create_line(self, ix0, iy0, ix1, iy1, fill=self.slice_color, width=LINE_WIDTH)
        else:
            Canvas.create_line(self, ix0, iy0, ix1, iy1, fill=line_color, width=LINE_WIDTH)
        
    def plotLeftPointingArrow(self, xval, yval, color='white', shaft_length=10, arrowhead_length=5):
        # Draw shaft of arrow
        ix0 = self.xtranslate(xval)
        iy0 = self.ytranslate(yval)
        ix1 = ix0 + shaft_length
        iy1 = iy0
        Canvas.create_line(self, ix0, iy0, ix1, iy1, fill=color, width=LINE_WIDTH)
        
        # Draw upper half of arrowhead
        ix1 = ix0 + arrowhead_length
        iy1 = iy0 + arrowhead_length
        Canvas.create_line(self, ix0, iy0, ix1, iy1, fill=color, width=LINE_WIDTH)
        
        # Draw lower half of arrowhead
        ix1 = ix0 + arrowhead_length
        iy1 = iy0 - arrowhead_length
        Canvas.create_line(self, ix0, iy0, ix1, iy1, fill=color, width=LINE_WIDTH)

    def plotPoint(self, xval, yval, color='yellow', radius=3):
        ok = True
        if xval < self.xlow or xval > self.xhigh:
            ok = False
        if yval < self.ylow or yval > self.yhigh:
            ok = False
        if ok:
            x = self.xtranslate(xval)
            y = self.ytranslate(yval)
            Canvas.create_oval(self, x-radius, y-radius, x+radius, y+radius, fill=color, outline=color)

    def plotAllPoints(self):
        if len(self.points) > 0:
            for p in self.points[:-1]:
                self.plotPoint(p[0], p[1])
            p = self.points[-1]
            self.plotPoint(p[0], p[1], 'red')
        
    def plotPointOnFunc(self, f, xval):
        yval = f(xval)
        x = self.xtranlslate(xval)
        y = self.ytranslate(yval)
        Canvas.create_oval(self, x-2, y-2, x+2, y+2, fill='red', outline='white')
        
    def plotFunc(self, f, npoints):
        xval = self.xlow
        yval = f(xval)
        x0 = self.xtranslate(xval)
        y0 = self.ytranslate(yval)
        #print 'xlow=%f, xhigh=%f' % (self.xlow, self.xhigh)
        for i in range(npoints + 1):
            xval = self.xlow + (self.xhigh - self.xlow)*float(i)/float(npoints)
            #print 'xval=%f' % xval
            yval = f(xval)
            x = self.xtranslate(xval)
            y = self.ytranslate(yval)
            Canvas.create_line(self, x0, y0, x, y, fill=self.func_color, width=LINE_WIDTH)
            x0 = x
            y0 = y

