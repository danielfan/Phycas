# TODO:
# o need to abort if not reading a flex rates parameter file

import sys,math
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch

class ColumnPlot():
    def __init__(self, fn, max_rate, max_prob):
        self.max_rate = max_rate
        self.max_prob = max_prob
        self.page_width, self.page_height = letter

        # create a pdf file (a "canvas")
        self.canvas = canvas.Canvas(fn, pagesize=letter)

        #fonts = self.canvas.getAvailableFonts()
        #for f in fonts:
        #    print f

        self.fontSize = 12
        self.fontName = 'Helvetica'
        self.canvas.setFont(self.fontName, self.fontSize)        
        self.M_width  = self.canvas.stringWidth('M', self.fontName, self.fontSize)
        self.M_height = self.fontSize*inch/72.0

        self.bar_width = 0.25*inch
        self.tick_len = 0.1*inch

        # plot is 6 in wide and 4 in high
        self.x0 = 1.25*inch   # 8.5 - 6 = 2.5, left/right margin both 1.25
        self.y0 = 3.5*inch    # 11 - 4 = 7, top/bottom margin both 3.5
        self.prop_plot_width_used = 0.9
        self.prop_plot_height_used = 0.9
        self.plot_width = 6.0*inch
        self.prop_plot_width_used = 0.9
        self.plot_height = 4.0*inch
        self.prop_plot_height_used = 0.9
        self.x1 = self.page_width - 1.25*inch
        self.y1 = self.page_height - 3.5*inch

        # draw vertical axis at x = x0
        self.canvas.setStrokeColorRGB(0, 0, 0)    # black
        self.canvas.setFillColorRGB(0, 0, 0)    # black
        self.canvas.line(self.x0, self.y0, self.x0, self.y1)

        # draw horizontal axis at y = y0
        self.canvas.setStrokeColorRGB(0, 0, 0)    # black
        self.canvas.setFillColorRGB(0, 0, 0)    # black
        self.canvas.line(self.x0, self.y0, self.x1, self.y0)

        # scale horizontal axis so that it extends a little beyond max_rate
        # want xscale*max_rate = plot_width*prop_plot_width_used
        self.xscale = self.plot_width*self.prop_plot_width_used/self.max_rate

        # scale vertical axis so that it extends a little beyond max_prob    
        # want yscale*max_prob = plot_height*prop_plot_height_used
        self.yscale = self.plot_height*self.prop_plot_height_used/self.max_prob

        # label x axis
        s = '%.5f' % (0.5*self.max_rate)
        x = self.x0 + 0.5*self.max_rate*self.xscale
        y = self.y0 - 1.5*self.M_height
        self.canvas.drawCentredString(x, y, s)
        y = self.y0
        self.canvas.line(x, y, x, y - self.tick_len)

        s = '%.5f' % (self.max_rate)
        x = self.x0 + self.max_rate*self.xscale
        y = self.y0 - 1.5*self.M_height
        self.canvas.drawCentredString(x, y, s)
        y = self.y0
        self.canvas.line(x, y, x, y - self.tick_len)
        
        # label y axis
        s = '%.5f' % (0.5*self.max_prob)
        x = self.x0 - self.M_width
        y = self.y0 + 0.5*self.max_prob*self.yscale
        self.canvas.drawRightString(x, y - self.M_height/2.0, s)
        x = self.x0
        self.canvas.line(x, y, x - self.tick_len, y)
        
        s = '%.5f' % (self.max_prob)
        x = self.x0 - self.M_width
        y = self.y0 + self.max_prob*self.yscale
        self.canvas.drawRightString(x, y - self.M_height/2.0, s)
        x = self.x0
        self.canvas.line(x, y, x - self.tick_len, y)

    def plotOneSample(self, ncat, rates, probs):
        for i in range(ncat):
            r = rates[i]
            rx = self.xscale*r
            p = probs[i]
            py = self.yscale*p
            x = self.x0 + rx
            y = self.y0
            self.canvas.rect(x, y, self.bar_width, py, stroke=1, fill=0)

    def finalize(self):
        # finalize the page
        self.canvas.showPage()

        # close the pdf file    
        self.canvas.save()

def equalWidthColumnsMethod(fn):    
    lines = file(fn, 'r').readlines()
    headers = lines[1].split()
    ncat_pos = headers.index('ncat')
    first_rate_pos = headers.index('rates_probs')
    max_rate = 0.0
    max_prob = 0.0
    samples = []
    for line in lines[2:]:
        parts = line.split()
        ncat = int(parts[ncat_pos])
        gen = int(parts[0])
        rates = []
        probs = []
        for i in range(ncat):
            rate = float(parts[first_rate_pos + i])
            if rate > max_rate:
                max_rate = rate
            rates.append(rate)
            prob = float(parts[first_rate_pos + i + ncat])
            if prob > max_prob:
                max_prob = prob
            probs.append(prob)
        samples.append((ncat, rates, probs))

    plotter = ColumnPlot('flexplot.pdf', max_rate, max_prob)
    for ncat, rates, probs in samples:
        plotter.plotOneSample(ncat, rates, probs)
    plotter.finalize()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Specify a (parameter) file name when invoking flexplot'
        print 'e.g. python flexplot.py myfile.nex.p'
        sys.exit()
    fn = sys.argv[1]
    print 'Analyzing file',fn

    equalWidthColumnsMethod(fn)    
