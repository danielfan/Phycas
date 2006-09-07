from Tkinter import *
from tkMessageBox import askokcancel
from tkFileDialog import askopenfilename
import ImageTk,Image,JpegImagePlugin
import ReadNexus
import os
	
useImage = True
class NexusReaderTester(Frame):
    def __init__(self, parent=None):
        # call the Frame constructor
        Frame.__init__(self, parent)
        self.pack()

        # create a Canvas on which to show the phycas logo
        fileString= os.path.join(os.environ.get('PHYCAS_ROOT'), 'images', 'phycaslogo.jpg')
        #print fileString
        if useImage:
        	self.imgobj = ImageTk.PhotoImage(file = fileString)
        	self.imgcanvas = Canvas(self, width=167, height=100, bg='red')
        	self.imgcanvas.pack()
        	self.imgcanvas.create_image(0, 0, image=self.imgobj, anchor=NW)

        # add a Button so that user can browse for a nexus file to open
        buttonString = useImage and 'Browse for nexus file...' or str(fileString)
        Button(self, text = buttonString, command=self.findnexus, width=25).pack(side=TOP)

    def findnexus(self):
        # get a file name from the user
        fn = askopenfilename(initialdir='.')

        # call the embedded NCL code to read and store the file
        r = ReadNexus.NexusReader(4)
        #r.ReadFile(r'C:\Synchronized\Phycas\pythonextend\readnxs\nyldna4.nex')
        try:
            r.ReadFile(fn)
        except RuntimeError, e:
            askokcancel('Error:', '%s' % e.__str__())
        else:
            askokcancel('Success!', 'Number of characters read = %d' % r.GetNChar())

if __name__ == '__main__':
    root = Tk()
    app = NexusReaderTester(root)
    root.mainloop()
