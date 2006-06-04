# A tree viewer to use for debugging purposes, especially for debugging likelihood
# calculations and MCMC moves involving trees. Shows a graphical representation of 
# the tree as it is laid out in memory. Right now (June 4, 2006), only a hard-coded
# default tree is shown, but eventually one will be able to pass in a tree complete
# with conditional likelihood arrays (CLAs) and see, graphically, which CLAs are
# currently valid, which need to be recalculated, which ones are flagged as needing
# recalcuation, and which node is the current likelihood root node. In addition,
# nodes (and their associated edges) that are selected will be shown in a different
# color, so it will be possible to see which edges are affected by moves such as the
# Larget-Simon move.

from Phycas import *
from Tkinter import *
from tkMessageBox import askokcancel, askquestion
from tkSimpleDialog import askstring
from tkFileDialog import askopenfilename
import tkFont
import math

class TreeCanvas(Canvas):
    def xtranslate(self, x):
        return int(self.left + self.xscaler*x)

    def ytranslate(self, y):
        return int(self.bottom - self.yscaler*y)

    def plotNode(self, xval, yval, name, color='yellow', radius=2):
        x = self.xtranslate(xval)
        y = self.ytranslate(yval)
        #w = self.font.measure(name)
        #Canvas.create_oval(self, x-w, y-w, x+w, y+w, fill='black', outline='white')
        Canvas.create_text(self, x, y, text=name, font=self.font, fill=self.nodename_text_color)

    def plotEdge(self, fromx, fromy, tox, toy, radius=10, color='yellow'):
        x0 = self.xtranslate(fromx)
        y0 = self.ytranslate(fromy)
        x = self.xtranslate(tox)
        y = self.ytranslate(toy)

        # Leave space at ends for node identifier
        if x == x0:
            if y > y0:
                theta = math.pi/2.0
            else:
                theta = -math.pi/2.0
        else:
            theta = math.atan(float(y - y0)/float(x - x0))
            if x < x0:
                theta += math.pi
        dx = float(radius)*math.cos(theta)
        dy = float(radius)*math.sin(theta)
        
        Canvas.create_line(self, x0+dx, y0+dy, x-dx, y-dy, fill=color)

    def drawTree(self):
        # Only recalculate tree height if tree changes (don't want to do this 
        # extra work if user simply resizes the window, for example)
        if self.rescale:
            self.tree_height = self.tree.calcTotalHeight()
            self.yscaler = self.usableh/self.tree_height
            self.tree_width = float(self.tree.getNTips() - 2)
            self.xscaler = self.usablew/self.tree_width

        # Do a postorder traversal to draw the tree
        x = self.tree_width
        nd = self.tree.getFirstPostorder()
        while not nd.isRoot():
            if nd.isTip():
                nd.setY(0.0)
                nd.setX(x)
                x -= 1.0
            else:
                lchild_x = nd.getLeftChild().getX()
                rchild_x = nd.getX()
                nd.setX((lchild_x + rchild_x)/2.0)

            # If nd is rightmost child of its parent, initialize parent's x and y values
            parent = nd.getParent()
            if not nd.getRightSib():
                parent.setX(nd.getX())
                parent.setY(0.0)

            # Make sure that parent's y equals total distance from parent to
            # furthest tip in the lineage that includes nd
            nd_edgelen = (self.tree.hasEdgeLens() and nd.getEdgeLen() or 1.0)
            height_thru_nd = nd_edgelen
            if not nd.isTip():
                height_thru_nd += nd.getY()
            if height_thru_nd > parent.getY():
                parent.setY(height_thru_nd)

            nd = nd.getNextPostorder()

        # Proceed up tree in preorder fashion to draw
        root = self.tree.getFirstPreorder()
        x = root.getLeftChild().getX()
        y = 0.0
        self.plotNode(x, y, root.getNodeName(), self.node_unselected_fill_color)
        nd = nd.getNextPreorder()

        while nd:
            parent = nd.getParent()
            x0 = parent.getX()
            y0 = self.tree_height - parent.getY()
            nd_edgelen = (self.tree.hasEdgeLens() and nd.getEdgeLen() or 1.0)
            nd.setY(parent.getY() - nd_edgelen)
            x = nd.getX()
            y = self.tree_height - nd.getY()
            self.plotEdge(x0, y0, x, y, 3*self.font_Mwidth//2, self.edge_unselected_color)
            self.plotNode(x, y, nd.getNodeName(), self.node_unselected_fill_color)
            nd = nd.getNextPreorder()
        
    def __init__(self, parent):
        self.background_color = 'blue'

        self.rescale = True
        self.tree_height = 0.0
        self.tree_width = 0.0
        self.xscaler = 1.0
        self.yscaler = 1.0
        self.plot_margin = 20        

        # variables associated with the tree being displayed
        self.default_tree_topology = '(1,(((6,7)5,(9,(11,12)10)8)4,13)3,(15,(17,(19,20)18,21)16)14)2'
        self.tree = Phylogeny.Tree()
        self.tree.buildFromString(self.default_tree_topology)

        # variables associated with showing edges
        self.edge_selected_color = 'cyan'
        self.edge_unselected_color = 'cyan'

        # variables associated with showing nodes
        self.node_selected_fill_color = 'red'
        self.node_unselected_fill_color = 'cyan'
        self.node_stroke_color = 'cyan'

        # variables associated with showing status of conditional likelihood arrays
        self.CLA_radius = 2             # radius of circle plotted for each CLA
        self.CLA_valid_color = 'green'  # color indicating a CLA that is currently valid
        self.CLA_invalid_color = 'red'  # color indicating a CLA that needs to be recalculated
        
        Canvas.__init__(self, master=parent, bg=self.background_color)

        # font-related
        self.font = tkFont.Font(family='Courier', size=8)
        fm = self.font.metrics()
        self.font_height = fm['ascent'] + fm['descent']
        self.font_Mwidth = self.font.measure('M')
        self.nodename_text_color = 'cyan'

        self.bind("<Configure>", self.resizeEvent)

    def resizeEvent(self, event):
        self.resize(event.width, event.height, self.plot_margin)
        self.repaint()

    def resize(self, new_width, new_height, new_margin):
        self.plotw = new_width
        self.ploth = new_height
        self.plotm = new_margin
        assert self.font_height < new_margin, 'height of font too large for specified plot margin'
        #self.offset = new_margin/2
        self.usablew = self.plotw - 2*self.plotm
        self.usableh = self.ploth - 2*self.plotm
        self.left = self.plotm
        self.top = self.plotm
        self.right = self.plotw - self.plotm
        self.bottom = self.ploth - self.plotm
        self.hcenter = self.plotw/2
        self.vcenter = self.ploth/2
        #Canvas.config(self, width=new_width, height=new_height)        
        
    def repaint(self):
        Canvas.create_rectangle(self, 0, 0, self.plotw, self.ploth, fill=self.background_color);
        self.drawTree()

    def reset(self):
        import gc       # garbage collector
        gc.collect()    # should capture return value, which is number of unreachable objects found
        self.repaint()

class TreeViewer(Frame):
    def __init__(self, parent=None):
        Frame.__init__(self, parent, bg='yellow')
        self.pack(expand=YES, fill=BOTH)

        # variables related to how much is displayed
        self.show_node_numbers = False
        self.show_CLAs = False

        # create a frame to hold the menu buttons
        menuf = Frame(self, bg='magenta')
        menuf.pack(expand=NO, fill=X)
        
        # create the File menu button
        self.filemb = Menubutton(menuf, text='File', relief=RAISED, anchor=W, borderwidth=0)
        self.filemb.pack(expand=NO, fill=X, side=LEFT)
        self.filemb.menu = Menu(self.filemb, tearoff=0)
        self.filemb['menu'] = self.filemb.menu
        self.filemb.menu.add_command(label='Quit', command=self.quit)
        
        # create the Options menu button
        self.samplemb = Menubutton(menuf, text='Options', relief=RAISED, anchor=W, borderwidth=0)
        self.samplemb.pack(expand=YES, fill=X, side=LEFT)
        self.samplemb.menu = Menu(self.samplemb, tearoff=0)
        self.samplemb['menu'] = self.samplemb.menu
        self.samplemb.menu.add_command(label='Toggle node numbers', command=self.toggleNodeNumbers)
        self.samplemb.menu.add_command(label='Toggle CLAs', command=self.toggleCLAs)

        # create the canvas
        self.plotter = TreeCanvas(parent=self)
        self.plotter.pack(side=TOP, expand=YES, fill=BOTH)
        
        # create the status label
        self.status_label = Label(self, justify=LEFT, relief=SUNKEN, height=1, anchor=W, text='Ready')
        self.status_label.pack(side=TOP, expand=NO, fill=X)
        
        # bind some keys to the application (doesn't matter which widget has focus when you use bind_all
        #self.bind_all("<KeyPress-o>", self.keybdOverrelaxedOneStep)
        #self.bind_all("<KeyPress-n>", self.keybdOneStep)
        #self.bind_all("<Shift-KeyPress-N>", self.keybdManySteps)
        #self.bind_all("<KeyPress-a>", self.keybdAdaptSampler)
        #self.bind_all("<KeyPress-r>", self.keybdReset)
        #self.bind_all("<KeyPress-q>", self.keybdQuit)

        # configure event is bound only to the main frame
        #self.bind("<Configure>", self.resizing)

    def toggleNodeNumbers(self):
        if self.show_node_numbers:
            self.show_node_numbers = False
            self.status_label.config(text='node numbers hidden')
        else:
            self.show_node_numbers = True
            self.status_label.config(text='node numbers visible')
        
    def toggleCLAs(self):
        if self.show_CLAs:
            self.show_CLAs = False
            self.status_label.config(text='CLAs hidden')
        else:
            self.show_CLAs = True
            self.status_label.config(text='CLAs visible')

if __name__ == '__main__':
    TreeViewer().mainloop()
