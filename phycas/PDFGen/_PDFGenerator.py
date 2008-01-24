import os.path, re

class PDFObject(object):
    def __init__(self, pdf_generator, object_type):
        self.pdf = pdf_generator
        self.type = object_type
        self.number = self.pdf.getNextObjectNumber()
        self.pdf.objects.append(self)
    
class PDFCatalogObject(PDFObject):
    def __init__(self, pdf_generator):
        PDFObject.__init__(self, pdf_generator, 'Catalog')
        self.outlines_obj = 2
        self.pages_obj = 3

    def write(self, outf, terminator):
        outstr = '%d 0 obj%c' % (self.number, terminator)
        outstr += '    <<%c' % terminator
        outstr += '        /Type /%s%c' % (self.type, terminator)
        outstr += '        /Outlines %d 0 R%c' % (self.outlines_obj, terminator)
        outstr += '        /Pages %d 0 R%c' % (self.pages_obj, terminator)
        outstr += '    >>%c' % terminator
        outstr += '    endobj%c' % terminator
        outf.write(outstr)
        return len(outstr)

class PDFOutlinesObject(PDFObject):
    def __init__(self, pdf_generator):
        PDFObject.__init__(self, pdf_generator, 'Outlines')

    def write(self, outf, terminator):
        outstr = '%d 0 obj%c' % (self.number, terminator)
        outstr += '    <<%c' % terminator
        outstr += '        /Type /%s%c' % (self.type, terminator)
        outstr += '        /Count 0%c' % terminator
        outstr += '    >>%c' % terminator
        outstr += '    endobj%c' % terminator
        outf.write(outstr)
        return len(outstr)

class PDFProcSetObject(PDFObject):
    def __init__(self, pdf_generator):
        PDFObject.__init__(self, pdf_generator, 'ProcSet')

    def write(self, outf, terminator):
        outstr = '%d 0 obj%c' % (self.number, terminator)
        outstr += '    [/PDF /Text]%c' % terminator
        outstr += '    endobj%c' % terminator
        outf.write(outstr)
        return len(outstr)

class PDFPagesObject(PDFObject):
    def __init__(self, pdf_generator):
        PDFObject.__init__(self, pdf_generator, 'Pages')
        self.kids = []

    def addPageObj(self, page_obj):
        self.kids.append('%d 0 R' % page_obj)

    def write(self, outf, terminator):
        outstr = '%d 0 obj%c' % (self.number, terminator)
        outstr += '    <<%c' % terminator
        outstr += '        /Type /%s%c' % (self.type, terminator)
        outstr += '        /Kids [ %s ]%c' % (' '.join(self.kids), terminator)
        outstr += '        /Count %d%c' % (len(self.kids), terminator)
        outstr += '    >>%c' % terminator
        outstr += '    endobj%c' % terminator
        outf.write(outstr)
        return len(outstr)

class PDFPageObject(PDFObject):
    def __init__(self, pdf_generator, parent_page):
        PDFObject.__init__(self, pdf_generator, 'Page')
        self.parent_num = parent_page.number
        self.contents = []
        self.fonts_used = {}    # fonts_used.keys() can be used to avoid duplicates

    def addContent(self, obj):
        self.contents.append('%d 0 R' % obj.number)
        if obj.type == 'Text':
            self.fonts_used[(obj.font_num, obj.font_obj)] = 1

    def write(self, outf, terminator):
        outstr = '%d 0 obj%c' % (self.number, terminator)
        outstr += '    <<%c' % terminator
        outstr += '        /Type /%s%c' % (self.type, terminator)
        outstr += '        /Parent %d 0 R%c' % (self.parent_num, terminator)
        outstr += '        /MediaBox [ 0 0 %d %d ]%c' % (self.pdf.page_width, self.pdf.page_height, terminator)
        if len(self.contents) == 1:
            outstr += '        /Contents %s%c' % (self.contents[0], terminator)
        elif len(self.contents) > 1:
            outstr += '        /Contents [ %s ]%c' % (' '.join(self.contents), terminator)
        outstr += '        /Resources << /ProcSet %d 0 R%c' % (self.pdf.getProcSetObjNumber(), terminator)
        if len(self.fonts_used) > 0:
            # Note: do not use separate /Font entries for each font. Instead, put them all in one /Font entry
            # Using separate /Font entries will result in a message like this in Acrobat Reader:
            # "Can not find a font in the Resource dictionary - using Helvetica instead"
            outstr += '                      /Font << %c' % terminator
            for fnum,onum in self.fonts_used.keys():
                outstr += '                               /F%d %d 0 R%c' % (fnum, onum, terminator)
            outstr += '                            >>%c' % terminator
        outstr += '                   >>%c' % terminator
        outstr += '    >>%c' % terminator
        outstr += '    endobj%c' % terminator
        outf.write(outstr)
        return len(outstr)

class PDFFontObject(PDFObject):
    def __init__(self, pdf_generator, font_family):
        PDFObject.__init__(self, pdf_generator, 'Font')
        self.font_family = font_family
        self.font_num = self.pdf.getNextFontNumber()
        self.pdf.fonts.append(self)

    def write(self, outf, terminator):
        outstr = '%d 0 obj%c' % (self.number, terminator)
        outstr += '    <<%c' % terminator
        outstr += '        /Type /%s%c' % (self.type, terminator)
        outstr += '        /Subtype /Type1%c' % terminator
        outstr += '        /Name /F%d%c' % (self.font_num, terminator)
        outstr += '        /BaseFont /%s%c' % (self.font_family, terminator)
        #outstr += '        /Encoding /MacRomanEncoding%c' % terminator
        outstr += '        /Encoding /WinANSIEncoding%c' % terminator
        outstr += '    >>%c' % terminator
        outstr += '    endobj%c' % terminator
        outf.write(outstr)
        return len(outstr)

class PDFTextObject(PDFObject):
    def __init__(self, pdf_generator, x, y, font_family, font_size, text):
        PDFObject.__init__(self, pdf_generator, 'Text')
        self.font_num, self.font_obj = self.pdf.findFont(font_family)
        self.font_size = font_size
        self.x = int(x)
        self.y = int(y)
        self.text = text

    def write(self, outf, terminator):
        txtstr  = 'BT%c' % terminator
        txtstr += '    /F%d %d Tf%c' % (self.font_num, self.font_size, terminator)
        txtstr += '    %d %d Td%c' % (self.x, self.y, terminator)
        txtstr += '    (%s) Tj%c' % (self.text, terminator)
        txtstr += 'ET'
        
        outstr = '%d 0 obj%c' % (self.number, terminator)
        outstr += '    <<%c' % terminator
        outstr += '        /Length %d%c' % (len(txtstr), terminator)
        outstr += '    >>%c' % terminator
        outstr += 'stream%c' % terminator
        outstr += '%s%c' % (txtstr, terminator)
        outstr += 'endstream%c' % terminator
        outstr += '    endobj%c' % terminator
        outf.write(outstr)
        return len(outstr)

class PDFLineObject(PDFObject):
    def __init__(self, pdf_generator, x0, y0, x1, y1, line_width, line_style):
        PDFObject.__init__(self, pdf_generator, 'Line')
        self.x0 = float(x0)
        self.y0 = float(y0)
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.width = float(line_width)
        self.dash_pattern = line_style == 'dotted' and '[3] 0 d' or '[] 0 d'
        self.line_cap_style = 1 # 0 = square ends; 1 = rounded ends; 2 = projecting square ends 

    def write(self, outf, terminator):
        objstr  = '    %.1f w%c' % (self.width, terminator)
        objstr += '    %s%c' % (self.dash_pattern, terminator)
        objstr += '    %d J%c' % (self.line_cap_style, terminator)
        objstr += '    %.1f %.1f m%c' % (self.x0, self.y0, terminator)
        objstr += '    %.1f %.1f l%c' % (self.x1, self.y1, terminator)
        objstr += '    S'
        
        outstr = '%d 0 obj%c' % (self.number, terminator)
        outstr += '    <<%c' % terminator
        outstr += '        /Length %d%c' % (len(objstr), terminator)
        outstr += '    >>%c' % terminator
        outstr += 'stream%c' % terminator
        outstr += '%s%c' % (objstr, terminator)
        outstr += 'endstream%c' % terminator
        outstr += '    endobj%c' % terminator
        outf.write(outstr)
        return len(outstr)

class PDFRectangleObject(PDFObject):
    def __init__(self, pdf_generator, x0, y0, width, height, line_width, line_style):
        PDFObject.__init__(self, pdf_generator, 'Rectangle')
        self.x0 = float(x0)
        self.y0 = float(y0)
        self.w = float(width)
        self.h = float(height)
        self.lw = float(line_width)
        self.dash_pattern = line_style == 'dotted' and '[3] 0 d' or '[] 0 d'
        self.line_cap_style = 1 # 0 = square ends; 1 = rounded ends; 2 = projecting square ends 

    def write(self, outf, terminator):
        objstr  = '    %.1f w%c' % (self.lw, terminator)
        objstr += '    %d J%c' % (self.line_cap_style, terminator)
        objstr += '    %s%c' % (self.dash_pattern, terminator)
        objstr += '    %.1f %.1f %.1f %.1f re%c' % (self.x0, self.y0, self.w, self.h, terminator)
        objstr += '    S'
        
        outstr = '%d 0 obj%c' % (self.number, terminator)
        outstr += '    <<%c' % terminator
        outstr += '        /Length %d%c' % (len(objstr), terminator)
        outstr += '    >>%c' % terminator
        outstr += 'stream%c' % terminator
        outstr += '%s%c' % (objstr, terminator)
        outstr += 'endstream%c' % terminator
        outstr += '    endobj%c' % terminator
        outf.write(outstr)
        return len(outstr)

class PDFGenerator(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Class that is responsible for creating PDF files. This PDF generator
    is very lightweight, having only the capabilities needed to draw
    trees. It is also not very efficient (separate objects are created
    for each line, text string, etc., rather than combining them into one
    stream object. It would not be too difficult, however, to join
    objects if efficiency becomes a problem.

    """
    def __init__(self, page_width_inches, page_height_inches):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes data members.

        >>> from phycas import *
        >>> inch = 72.0
        >>> pdf = PDFGen.PDFGenerator(8.5, 11.0)
        >>> pdf.overwrite = True
        >>> pdf.newPage()
        >>> pdf.addText(3.0*inch, 6.0*inch, 'Times-Italic', 36, 'Hello, World!')
        >>> pdf.addLine(3.0*inch, 5.9*inch, 5.8*inch, 5.9*inch, 4)
        >>> pdf.saveDocument('test.pdf')

        """
        self.objects     = []
        self.fonts       = []
        self.char_width  = {}
        self.xheight     = {}
        self.filename    = None
        self.overwrite   = False
        self.terminator  = '\x0A' # line feed
        self.inch        = 72.0
        self.page_width  = int(self.inch*page_width_inches)   # 612 equals (72 points/inch)*(8.5 inches)
        self.page_height = int(self.inch*page_height_inches)  # 792 equals (72 points/inch)*(11 inches)

        # Create a catalog object        
        self.catalog    = PDFCatalogObject(self)
        
        # Create an outlines (a.k.a. bookmarks) root object
        self.outlines   = PDFOutlinesObject(self)
        self.catalog.outlines_obj = self.outlines.number

        # Create a procset object (obsolete, included for backward compatability)
        self.procset    = PDFProcSetObject(self)

        # Create a pages object
        self.pages      = PDFPagesObject(self)
        self.catalog.pages_obj = self.pages.number

        self.curr_page = None

    def newPage(self):                
        self.curr_page = PDFPageObject(self, self.pages)
        self.pages.addPageObj(self.curr_page.number)

    def addText(self, x, y, font_family, font_size, text):
        assert self.curr_page, 'call newPage() function before adding first content'
        text = PDFTextObject(self, x, y, font_family, font_size, text)
        self.curr_page.addContent(text)

    def addLine(self, x0, y0, x1, y1, line_width, line_style = 'solid'):
        assert self.curr_page, 'call newPage() function before adding first content'
        line = PDFLineObject(self, x0, y0, x1, y1, line_width, line_style)
        self.curr_page.addContent(line)

    def addRectangle(self, x0, y0, w, h, line_width = 1, line_style = 'solid'):
        assert self.curr_page, 'call newPage() function before adding first content'
        rect = PDFRectangleObject(self, x0, y0, w, h, line_width, line_style)
        self.curr_page.addContent(rect)

    def getNextObjectNumber(self):
        return len(self.objects) + 1

    def getNextFontNumber(self):
        return len(self.fonts) + 1

    def getProcSetObjNumber(self):
        return self.procset.number

    def findFont(self, font_family):
        for f in self.fonts:
            if font_family == f.font_family:
                return (f.font_num, f.number)
        # font apparently does not yet exist, so add it now 
        f = PDFFontObject(self, font_family)
        return (f.font_num, f.number)

    def getXHeight(self, font_face):
        if not font_face in self.xheight:
            self.loadAFMFile(font_face)
        return self.xheight[font_face]
        
    def calcStringWidth(self, font_face, s):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes the width of string s using the specified font_face, which
        should be one of the standard 14 fonts: Times-Roman, Helvetica,
        Courier, Symbol, Times-Bold, Helvetica-Bold, Courier-Bold,
        ZapfDingbats, Times-Italic, Helvetica-Oblique, Courier-Oblique,
        Times-BoldItalic, Helvetica-BoldOblique, or Courier-BoldOblique.

        """
        if not font_face in self.char_width:
            self.loadAFMFile(font_face)
        wdict = self.char_width[font_face]
        wx = 0.0
        for c in s:
            #print c,'=',wdict[ord(c)],'-->',25.4*36.0*wdict[ord(c)]/72.0,'mm'
            wx += wdict[ord(c)]
        return wx

    def loadAFMFile(self, font_face):
        if font_face in self.char_width:
            assert font_face in self.xheight, 'character widths for font %s loaded but not xheight' % font_face
            return  # already been loaded
        mydir = os.path.dirname(os.path.abspath(__file__))
        filename = '%s/AFM/%s.afm' % (mydir,font_face)
        assert os.path.exists(filename), 'Font metrics file %s not found' % filename
        afmdata = open(filename,'r').read()
        # Find the xheight and store
        r = re.compile('XHeight (\d+)', re.MULTILINE)
        m = r.search(afmdata)
        if m:
            xheight = float(m.group(1))/1000.0
        else:
            xheight = 0.450 # use value for Times-Roman as default
        self.xheight[font_face] = xheight
        # Get all lines representing individual character metrics
        r = re.compile('StartCharMetrics\s+\d+(.+)EndCharMetrics', re.MULTILINE | re.DOTALL)
        m = r.search(afmdata)
        metrics = m.group(1) #.splitlines()
        r = re.compile('C (\d+) ; WX (\d+)', re.MULTILINE)
        m = r.findall(metrics)
        self.char_width[font_face] = {}
        for x in m:
            i = int(x[0])
            w = float(x[1])/1000.0
            if i > 0:
                self.char_width[font_face][i] = w

    def saveDocument(self, filename):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates the PDF document.

        """
        if os.path.exists(filename) and not self.overwrite:
            print "PDF file '%s' could not be saved because it already exists" % filename
            print 'Either set overwrite to True or delete/rename the existing file and try again'
            return
        outf = open(filename, 'wb')
        outstr = '%%PDF-1.4%c' % self.terminator
        outf.write(outstr)
        cum_bytes = len(outstr)
        offsets = []
        for o in self.objects:
            offsets.append(cum_bytes)
            cum_bytes += o.write(outf, self.terminator)
        pdf_size = len(offsets) + 1
        xrefstr = 'xref%c' % self.terminator
        xrefstr += '0 %d%c' % (pdf_size, self.terminator)
        xrefstr += '0000000000 65535 f %c' % self.terminator
        for o in offsets:
            xrefstr += '%010d 00000 n %c' % (o, self.terminator)
        outf.write(xrefstr)
        outf.write('trailer%c' % self.terminator)
        outf.write('    << /Size %d%c' % (pdf_size, self.terminator))
        outf.write('       /Root 1 0 R%c' % self.terminator)
        outf.write('    >>%c' % self.terminator)
        outf.write('startxref%c' % self.terminator)
        outf.write('%d%c' % (cum_bytes, self.terminator))
        outf.write('%%EOF')
        outf.close()

    def scatterPlot(self, data,
                    xinfo = (0.0, 1.0, 10, 1),  # (min, max, divisions, precision)
                    yinfo = (0.0, 1.0, 10, 1),  # (min, max, divisions, precision)
                    lines = True,
                    label_font = 'Helvetica',
                    font_height = 12,
                    line_width = 1,
                    paper_width_inches = 11.0,
                    paper_height_inches = 8.5,
                    right_margin_inches = 1.0,
                    left_margin_inches = 1.0,
                    top_margin_inches = 1.0,
                    bottom_margin_inches = 1.0):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a scatter plot using the supplied data. The xinfo and yinfo
        variables should both be 4-tuples of the form (minimum value, maximum
        value, number of divisions, number of decimal places). This
        information will be used to label the corresponding axis. For example,
        if xinfo = (0.0, 1.0, 10, 1) these labels would be evenly spaced along
        the x-axis: '0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7',
        '0.8', '0.9', '1.0'. The data variable should be in the form of a list
        of lists of tuples. Each outer list defines a separate line in the
        plot; e.g. [line1,line2,line3]. The next level defines the points in a
        line; e.g. line1 = [point1, point2, point3]. The final level defines
        the x- and y-coordinates of a point; e.g. point1 = (x,y). The values
        x and y should be in the same units used to label the axes. That is,
        if xmin = 0.0 and xmax = 1.0, then the x values should be between 0.0
        and 1.0.

        """
        xmin    = xinfo[0]
        xmax    = xinfo[1]
        xdiff   = xmax - xmin
        xdiv = xinfo[2]
        xfmt = '%%.%df' % xinfo[3]
        xlabels = [xfmt % (xmin + xdiff*float(i)/float(xdiv)) for i in range(xdiv + 1)]
        #print 'xlabels =',xlabels
        
        ymin    = yinfo[0]
        ymax    = yinfo[1]
        ydiff   = ymax - ymin
        ydiv = yinfo[2]
        yfmt = '%%.%df' % yinfo[3]
        ylabels = [yfmt % (ymin + ydiff*float(i)/float(ydiv)) for i in range(ydiv + 1)]
        #print 'ylabels =',ylabels
        
        self.newPage()
        
        inch          = 72.0
        spacer        = 5.0
        half_tick     = 2.0

        xlabel_height = font_height*self.getXHeight(label_font)
        xlabel_width = 0.0
        for label in xlabels:
            width = font_height*self.calcStringWidth(label_font, label)
            if width > xlabel_width:
                xlabel_width = width

        ylabel_width = 0.0
        for label in ylabels:
            width = font_height*self.calcStringWidth(label_font, label)
            if width > ylabel_width:
                ylabel_width = width
        
        right_margin  = inch*right_margin_inches
        left_margin   = inch*left_margin_inches
        top_margin    = inch*top_margin_inches
        bottom_margin = inch*bottom_margin_inches

        plot_bottom   = bottom_margin + xlabel_height + 2.0*spacer
        plot_left     = left_margin + ylabel_width + spacer
        plot_width    = inch*paper_width_inches - plot_left - right_margin - xlabel_width/2.0
        plot_height   = inch*paper_height_inches - plot_bottom - top_margin - xlabel_height/2.0

        #self.addRectangle(left_margin, bottom_margin, inch*paper_width_inches - left_margin - right_margin, inch*paper_height_inches - top_margin - bottom_margin, 1, 'dotted')

        # draw the x-axis
        self.addLine(plot_left, plot_bottom, plot_left + plot_width, plot_bottom, line_width)

        # draw the x-axis tick marks and labels for 0, ..., ntrees
        for i,label in enumerate(xlabels):
            x = plot_width*float(i)/float(len(xlabels) - 1)
            half_width = font_height*self.calcStringWidth(label_font, label)/2.0
            self.addText(plot_left + x - half_width, plot_bottom - xlabel_height - 2.0*spacer, label_font, font_height, label)
            self.addLine(plot_left + x, plot_bottom - half_tick, plot_left + x, plot_bottom + half_tick, line_width)
        
        # draw the y-axis
        self.addLine(plot_left, plot_bottom, plot_left, plot_bottom + plot_height, line_width)

        # draw the y-axis tick marks and labels for 0.0, 0.1, ..., 1.0
        for i,label in enumerate(ylabels):
            y = plot_height*float(i)/float(len(ylabels) - 1)
            ylabel = y - xlabel_height/2.0
            self.addText(left_margin, plot_bottom + ylabel, label_font, font_height, label)
            self.addLine(plot_left - half_tick, plot_bottom + y, plot_left + half_tick, plot_bottom + y, line_width)
            
        # draw points or lines
        for line_data in data:
            first_point = True
            for point in line_data:
                if first_point:
                    x0 = plot_width*(point[0] - xmin)/xdiff
                    y0 = plot_height*(point[1] - ymin)/ydiff
                    #print 'x0 =',point[0],', y0 =',point[1]
                    first_point = False
                else:
                    x = plot_width*(point[0] - xmin)/xdiff
                    y = plot_height*(point[1] - ymin)/ydiff
                    self.addLine(plot_left + x0, plot_bottom + y0, plot_left + x, plot_bottom + y, line_width)
                    #print 'x  =',point[0],', y  =',point[1]
                    x0 = x
                    y0 = y
            #raw_input('check')

if __name__ == '__main__':
    # The 14 standard fonts guaranteed to be available in all PDF consumer applications:
    #   Times-Roman      Helvetica             Courier             Symbol
    #   Times-Bold       Helvetica-Bold        Courier-Bold        ZapfDingbats
    #   Times-Italic     Helvetica-Oblique     Courier-Oblique
    #   Times-BoldItalic Helvetica-BoldOblique Courier-BoldOblique
    inch = 72.0
    page_width = 8.5*inch
    page_height = 11.0*inch
    pdf = PDFGenerator(page_width, page_height)
    #pdf.overwrite = True
    #pdf.newPage()
    #pdf.addText(3*inch, 6*inch, 'Times-Italic', 36, 'Hello, World!')
    #pdf.addLine(3*inch, 5.9*inch, 5.8*inch, 5.9*inch, 4)
    #pdf.saveDocument('test.pdf')
    print 'width of "Hello there!" in 12 pt Times-Italic',12.0*pdf.calcStringWidth('Times-Italic','Hello there!')
    print 'width of "Hello there!" in 12 pt Helvetica',12.0*pdf.calcStringWidth('Helvetica','Hello there!')
