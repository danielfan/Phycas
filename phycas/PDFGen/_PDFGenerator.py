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
            for fnum,onum in self.fonts_used.keys():
                outstr += '                      /Font << /F%d %d 0 R >>%c' % (fnum, onum, terminator)
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
        outstr += '        /Encoding /MacRomanEncoding%c' % terminator
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
    def __init__(self, pdf_generator, x0, y0, x1, y1, line_width):
        PDFObject.__init__(self, pdf_generator, 'Line')
        self.x0 = float(x0)
        self.y0 = float(y0)
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.width = float(line_width)
        self.line_cap_style = 1 # 0 = square ends; 1 = rounded ends; 2 = projecting square ends 

    def write(self, outf, terminator):
        objstr  = '    %.1f w%c' % (self.width, terminator)
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

    def addLine(self, x0, y0, x1, y1, line_width):
        assert self.curr_page, 'call newPage() function before adding first content'
        line = PDFLineObject(self, x0, y0, x1, y1, line_width)
        self.curr_page.addContent(line)

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
        afmdata = open(filename,'r').read()
        # Find the xheight and store
        r = re.compile('XHeight (\d+)', re.MULTILINE)
        m = r.search(afmdata)
        self.xheight[font_face] = float(m.group(1))/1000.0
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