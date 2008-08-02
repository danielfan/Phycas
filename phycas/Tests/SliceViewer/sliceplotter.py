from Tkinter import *
import tkFont

class Plotter(Canvas):
    def __init__(self, parent, func_to_plot):
        self.background_color = 'blue'
        self.axisfont_color = 'white'
        self.axis_color = 'cyan'
        self.slice_color = 'yellow'
        self.func_color = 'white'
        self.func = func_to_plot
        Canvas.__init__(self, master=parent, bg=self.background_color)

        self.font = tkFont.Font(family='Courier', size=8)
        fm = self.font.metrics()
        self.font_height = fm['ascent'] + fm['descent']

        self.points = []

        self.setXRange(0.0, 1.0, 0.2)
        self.setYRange(0.0, 2.0, 0.5)
        
        self.bind("<Configure>", self.resizeEvent)

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
        #Canvas.config(self, width=new_width, height=new_height)        
        
    def repaint(self):
        Canvas.create_rectangle(self, 0, 0, self.plotw, self.ploth, fill=self.background_color);
        self.drawAxes()
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
        Canvas.create_line(self, self.left, self.bottom, self.right, self.bottom, width=2, fill=self.axis_color)
        Canvas.create_line(self, self.left, self.bottom, self.left, self.top, width=2, fill=self.axis_color)
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

    def xtranslate(self, x):
        return int(self.left + self.xaxislen*((x - self.xlow)/(self.xhigh - self.xlow)))

    def ytranslate(self, y):
        return int(self.bottom - self.yaxislen*((y - self.ylow)/(self.yhigh - self.ylow)))

    def plotSlice(self, xleft, xright, xincr, y):
        x0 = self.xtranslate(xleft)
        x1 = self.xtranslate(xright)
        y = self.ytranslate(y)
        Canvas.create_line(self, x0, y, x1, y, fill=self.slice_color)
        n = int((xright - xleft)//xincr)
        ticks = [xleft + xincr*i for i in range(n + 1)]
        for tick in ticks:
            x = self.xtranslate(tick)
            Canvas.create_line(self, x, y-2, x, y+2, fill=self.slice_color)
        
    def plotLine(self, x0, y0, x1, y1, line_color=None):
        ix0 = self.xtranslate(x0)
        ix1 = self.xtranslate(x1)
        iy0 = self.ytranslate(y0)
        iy1 = self.ytranslate(y1)
        if line_color==None:
            Canvas.create_line(self, ix0, iy0, ix1, iy1, fill=self.slice_color, width=2)
        else:
            Canvas.create_line(self, ix0, iy0, ix1, iy1, fill=line_color, width=2)
        
    def plotPoint(self, xval, yval, color='yellow'):
        ok = True
        if xval < self.xlow or xval > self.xhigh:
            ok = False
        if yval < self.ylow or yval > self.yhigh:
            ok = False
        if ok:
            x = self.xtranslate(xval)
            y = self.ytranslate(yval)
            Canvas.create_oval(self, x-2, y-2, x+2, y+2, fill=color, outline=color)

    def plotAllPoints(self):
        for p in self.points:
            self.plotPoint(p[0], p[1])
        
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
            Canvas.create_line(self, x0, y0, x, y, fill=self.func_color)
            x0 = x
            y0 = y

