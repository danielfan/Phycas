import wx, os.path, re
import colorsys
from math import cos, sin, radians

class PlotData(object):
    def __init__(self):
        self.data                 = None
        self.title                = None
        self.xinfo                = None
        self.yinfo                = None
        self.lines                = None
        self.points               = None
        self.title_font_name      = None
        self.title_font_height    = None
        self.label_font_name      = None
        self.label_font_height    = None
        self.line_width           = None
        self.line_cap_style       = None
        self.paper_width_inches   = None
        self.paper_height_inches  = None
        self.right_margin_inches  = None
        self.left_margin_inches   = None
        self.top_margin_inches    = None
        self.bottom_margin_inches = None

class PlotPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1)
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.parent = parent
        
        self.title_font  = None       
        self.label_font  = None       
        
    def OnPaint(self, evt):
        if not the_app:
            return
        pdc = wx.PaintDC(self)
        try:
            dc = wx.GCDC(pdc)
        except:
            dc = pdc
            
        # Get dimensions of the client area of this window in pixels
        client_width, client_height = self.GetClientSizeTuple()
            
        # font families: DEFAULT, DECORATIVE, ROMAN, SCRIPT, SWISS, MODERN, TELETYPE
        # font styles: NORMAL, SLANT, ITALICS
        # font weights: NORMAL, LIGHT, BOLD
        if not self.title_font:
            self.title_font = wx.Font(the_app.pd.title_font_height, wx.SWISS, wx.NORMAL, wx.NORMAL, faceName=the_app.pd.title_font_name)
        if not self.label_font:
            self.label_font = wx.Font(the_app.pd.label_font_height, wx.SWISS, wx.NORMAL, wx.NORMAL, faceName=the_app.pd.label_font_name)

        self.redrawPlotFromExistingData(dc)
        
    def addText(self, dc, x, y, font, text):
        dc.SetFont(font)
        dc.DrawText(text, x, y)

    def addLine(self, dc, x0, y0, x1, y1, line_width, line_style = 'solid', cap_style = 'rounded'):
        dc.DrawLine(x0, y0, x1, y1)
        
    def addRectangle(self, dc, x0, y0, w, h, line_width = 1, line_style = 'solid'):
        dc.DrawRectangle(x0, y0, w, h)
        
    def getWidthAndHeight(self, dc, font, string):
        w, h, descent, external_leading = dc.GetFullTextExtent(string, font=font)
        return (w,h)
        
    def redrawPlotFromExistingData(self, dc):
        xmin    = the_app.pd.xinfo[0]
        xmax    = the_app.pd.xinfo[1]
        xdiff   = xmax - xmin
        xdiv = the_app.pd.xinfo[2]
        xfmt = '%%.%df' % the_app.pd.xinfo[3]
        xlabels = [xfmt % (xmin + xdiff*float(i)/float(xdiv)) for i in range(xdiv + 1)]
        #print 'xlabels =',xlabels
        
        ymin    = the_app.pd.yinfo[0]
        ymax    = the_app.pd.yinfo[1]
        ydiff   = ymax - ymin
        ydiv = the_app.pd.yinfo[2]
        yfmt = '%%.%df' % the_app.pd.yinfo[3]
        ylabels = [yfmt % (ymin + ydiff*float(i)/float(ydiv)) for i in range(ydiv + 1)]
        
        inch          = 72.0
        spacer        = 5.0
        half_tick     = 2.0

        title_height = 0.0
        title_width = 0.0
        if len(the_app.pd.title) > 0:
            title_width, title_height = self.getWidthAndHeight(dc, self.title_font, the_app.pd.title)

        xlabel_width = 0.0
        for label in xlabels:
            width, xlabel_height = self.getWidthAndHeight(dc, self.label_font, label) 
            if width > xlabel_width:
                xlabel_width = width

        ylabel_width = 0.0
        for label in ylabels:
            width, height = self.getWidthAndHeight(dc, self.label_font, label)
            if width > ylabel_width:
                ylabel_width = width
        
        right_margin  = inch*the_app.pd.right_margin_inches
        left_margin   = inch*the_app.pd.left_margin_inches
        top_margin    = inch*the_app.pd.top_margin_inches
        bottom_margin = inch*the_app.pd.bottom_margin_inches

        plot_bottom   = bottom_margin + xlabel_height + 2.0*spacer
        plot_left     = left_margin + ylabel_width + spacer
        plot_width    = inch*the_app.pd.paper_width_inches - plot_left - right_margin - xlabel_width/2.0
        plot_height   = inch*the_app.pd.paper_height_inches - plot_bottom - top_margin - xlabel_height/2.0 - 3.0*title_height

        # draw the title
        self.addText(dc, plot_left + (plot_width - title_width)/2.0, plot_bottom + plot_height + 2.0*title_height, self.title_font, the_app.pd.title)

        #test_str = '\xc4\xc5\xa7\xa8\xa9\xaa'
        #self.addText(dc, plot_left + (plot_width - title_width)/2.0, plot_bottom + plot_height + title_height, self.symbol_font, test_str)
        
        # draw the x-axis
        self.addLine(dc, plot_left, plot_bottom, plot_left + plot_width, plot_bottom, line_width = 1, cap_style = the_app.pd.line_cap_style)

        # draw the x-axis tick marks and labels
        for i,label in enumerate(xlabels):
            x = plot_width*float(i)/float(len(xlabels) - 1)
            width, height = self.getWidthAndHeight(dc, self.label_font, label)
            half_width = width/2.0
            self.addText(dc, plot_left + x - half_width, plot_bottom - xlabel_height - 2.0*spacer, self.label_font, label)
            self.addLine(dc, plot_left + x, plot_bottom - half_tick, plot_left + x, plot_bottom + half_tick, line_width = 1, cap_style = the_app.pd.line_cap_style)
        
        # draw the y-axis
        self.addLine(dc, plot_left, plot_bottom, plot_left, plot_bottom + plot_height, line_width = 1, cap_style = the_app.pd.line_cap_style)

        # draw the y-axis tick marks and labels
        for i,label in enumerate(ylabels):
            y = plot_height*float(i)/float(len(ylabels) - 1)
            ylabel = y - xlabel_height/2.0
            self.addText(dc, left_margin, plot_bottom + ylabel, self.label_font, label)
            self.addLine(dc, plot_left - half_tick, plot_bottom + y, plot_left + half_tick, plot_bottom + y, line_width = 1, cap_style = the_app.pd.line_cap_style)
            
        # draw points and/or lines
        for line_data in the_app.pd.data:
            first_point = True
            only_point = (len(line_data) == 1)
            for point in line_data:
                if first_point:
                    x0 = plot_width*(point[0] - xmin)/xdiff
                    y0 = plot_height*(point[1] - ymin)/ydiff
                    if the_app.pd.points:
                        self.addRectangle(dc, plot_left + x0 - spacer/2.0,
                                          plot_bottom + y0 - spacer/2.0,
                                          spacer,
                                          spacer,
                                          1,
                                          'solid')
                    if the_app.pd.lines and only_point and not the_app.pd.points:
                        # show something, if only just a line so short that it is effectively a point
                        self.addLine(dc, plot_left + x0, plot_bottom + y0, plot_left + x0 + the_app.pd.line_width/2.0, plot_bottom + y0, line_width = the_app.pd.line_width, cap_style = the_app.pd.line_cap_style)
                    #print 'x0 =',point[0],', y0 =',point[1]
                    first_point = False
                else:
                    x = plot_width*(point[0] - xmin)/xdiff
                    y = plot_height*(point[1] - ymin)/ydiff
                    if the_app.pd.lines:
                        self.addLine(dc, plot_left + x0, plot_bottom + y0, plot_left + x, plot_bottom + y, line_width = the_app.pd.line_width, cap_style = the_app.pd.line_cap_style)
                    if the_app.pd.points:
                        self.addRectangle(dc, plot_left + x - spacer/2.0,
                                          plot_bottom + y - spacer/2.0,
                                          spacer,
                                          spacer,
                                          1,
                                          'solid')
                    #print 'x  =',point[0],', y  =',point[1]
                    x0 = x
                    y0 = y
            #raw_input('check')

            
class wxAWTY(wx.Frame):
    def __init__(self, parent, ID, title, 
                pos   = (10, 10),   # wx.DefaultPosition, 
                size  = (800, 600),   # wx.DefaultSize,
                style = wx.DEFAULT_FRAME_STYLE):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes data members.

        """
        wx.Frame.__init__(self, parent, ID, title, pos, size, style)
        
        # Create status bar
        self.status_bar = self.CreateStatusBar()
        
        # Create menu bar 
        main_menu = wx.MenuBar() 
        
        # Create File menu with Exit item
        file_menu = wx.Menu()
        item = file_menu.Append(wx.ID_EXIT, 'E&xit', 'Quit program')
        self.Bind(wx.EVT_MENU, self.OnMenuExit, item)
        main_menu.Append(file_menu, '&File')
        
        # Create Help menu with About item
        help_menu = wx.Menu()
        item = help_menu.Append(wx.ID_ABOUT, '&About', 'About this program')
        self.Bind(wx.EVT_MENU, self.OnMenuAbout, item)
        main_menu.Append(help_menu, '&Help')
        
        # Make it so
        self.SetMenuBar(main_menu)

        # Create client window that will do all the heavy lifting
        self.panel = PlotPanel(self)
        
    def OnMenuExit(self, event):
        self.Destroy()
        
    def OnMenuAbout(self, event):
        help_text = 'AWTY within Phycas\nWritten by David L. Swofford, Paul O. Lewis and Mark T. Holder\nJuly, 2008'
        wx.MessageDialog(self, help_text, 'About "Are We There Yet?"', wx.OK).ShowModal()

class wxAWTYApp(wx.App):
    def __init__(self):
        # specify redirect=_defRedirect to get the default behavior, which is True for Windows/MacOS and False 
        # otherwise. I've set it to False here because the default causes error messages to die with the frame 
        # window, which only appears for an instant (see line 7732 of <site-packages>/wx-2.8-msw-unicode/wx/_core.py 
        # for the definition of class App)
        wx.App.__init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True)
        
    def OnInit(self):
        frame = wxAWTY(None, -1, "Are We There Yet?")
        self.SetTopWindow(frame)
        frame.Show(True)
        return True

    def scatterPlot(self, data,
                    title                = '',
                    xinfo                = (0.0, 1.0, 10, 1),  # (min, max, divisions, precision)
                    yinfo                = (0.0, 1.0, 10, 1),  # (min, max, divisions, precision)
                    lines                = True,
                    points               = True,
                    title_font           = 'Helvetica',
                    title_font_height    = 14,
                    label_font           = 'Helvetica',
                    font_height          = 12,
                    line_width           = 1,
                    line_cap_style       = 'rounded',
                    paper_width_inches   = 11.0,
                    paper_height_inches  = 8.5,
                    right_margin_inches  = 1.0,
                    left_margin_inches   = 1.0,
                    top_margin_inches    = 1.0,
                    bottom_margin_inches = 1.0):
        self.pd                      = PlotData()
        self.pd.data                 = data
        self.pd.title                = title
        self.pd.xinfo                = xinfo
        self.pd.yinfo                = yinfo
        self.pd.lines                = lines
        self.pd.points               = points
        self.pd.title_font_name      = title_font
        self.pd.title_font_height    = title_font_height
        self.pd.label_font_name      = label_font
        self.pd.label_font_height    = font_height
        self.pd.line_width           = line_width
        self.pd.line_cap_style       = line_cap_style
        self.pd.paper_width_inches   = paper_width_inches
        self.pd.paper_height_inches  = paper_height_inches
        self.pd.right_margin_inches  = right_margin_inches
        self.pd.left_margin_inches   = left_margin_inches
        self.pd.top_margin_inches    = top_margin_inches
        self.pd.bottom_margin_inches = bottom_margin_inches

# See <phycasdev>/phycas/Phycas/SumTImpl.py and search for CreateAWTYApp to find where this is used
# After the call to CreateAWTYApp, a call is made to wxAWTYApp::scatterPlot so that by the time 
# a window is created the data to be plotted are already in the structure the_app.pd
# There must be a better way to do this (i.e. without creating a global)
the_app = None
def CreateAWTYApp():
    global the_app
    the_app = wxAWTYApp()
    return the_app
    
if __name__ == '__main__':
    the_app = wxAWTYApp()
    the_app.MainLoop()

            
