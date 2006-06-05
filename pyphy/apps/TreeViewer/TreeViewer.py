# A tree viewer to use for debugging purposes, especially for debugging likelihood
# calculations and MCMC moves involving trees. Shows a graphical representation of 
# the tree as it is laid out in memory.
#
# Features:
# o Initially, node numbers are used to identify nodes, but pressing 'n' toggles
#   the display of node names rather than numbers
# o Nodes that are selected (along with their associated edges) are shown in red
#   (this is useful for showing, for example, which edges were modified by a
#   Larget-Simon move)
# o The label of the node currently serving at the likelihood root is shown in
#   magenta, whereas other nodes are gray (todo: allow nodes other than node 0 to
#   serve as the likelihood root)
# o A default tree is shown initially (todo: allow tree to be replaced)
# o Parental and filial conditional likelihood array status is indicated by colored
#   circles at the two ends of each edge (todo: always gray now, need to detect status)
#   Key to colors:
#     gray:  CLA status could not be determined (e.g. no TipData or InternalData structures)
#     green: CLA is valid
#     red:   CLA is invalid and will be recalculated upon next use

# Useful colors (for others, see http://www.mindspring.com/~squicker/colors.html):
white       = '#ffffff' 
black       = '#000000'
red         = '#ff0000' 
magenta     = '#ff00ff'
maroon      = '#800000'
green       = '#00ff00' 
dkgreen     = '#008000'
teal        = '#008080'
cyan        = '#00ffff'
blue        = '#0000ff'
purple      = '#800080'
navy        = '#000080'
midnight    = '#00009C'
gray        = '#808080' 
silver      = '#c0c0c0' 
brown       = '#5C3317'
olive       = '#808000'
yellow      = '#ffff00'

# The values of these variables determine the color scheme
color_plot_background   = midnight
color_undefined_cla     = silver
color_valid_cla         = green
color_invalid_cla       = red
color_selected_edge     = red
color_unselected_edge   = silver
color_selected_node     = red
color_unselected_node   = silver
color_likelihood_root   = magenta

from Phycas import *
from threading import *
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

    def plotNode(self, xval, yval, number, label, node_color=color_unselected_node, radius=2):
        x = self.xtranslate(xval)
        y = self.ytranslate(yval)
        color = (number == self.likelihood_root_nodenum and color_likelihood_root or node_color)
        Canvas.create_text(self, x, y, text=str(label), font=self.font, fill=color)

    def plotEdge(self, parent_x, parent_y, child_x, child_y, radius=10, edge_color=color_unselected_edge, parental_color=color_undefined_cla, filial_color=color_undefined_cla):
        x0 = self.xtranslate(parent_x)
        y0 = self.ytranslate(parent_y)
        x = self.xtranslate(child_x)
        y = self.ytranslate(child_y)

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
        
        # Draw the edge itself
        Canvas.create_line(self, x0+dx, y0+dy, x-dx, y-dy, fill=edge_color)

        # Draw the parental CLA marker
        Canvas.create_oval(self, x0+dx-2, y0+dy-2, x0+dx+2, y0+dy+2, fill=parental_color, outline=parental_color)

        # Draw the filial CLA marker
        Canvas.create_oval(self, x-dx-2, y-dy-2, x-dx+2, y-dy+2, fill=filial_color, outline=filial_color)

    def getEdgeLen(self, nd):
        return (self.tree.hasEdgeLens() and nd.getEdgeLen() or 1.0)

    def drawTree(self):
        # acquire will block other threads while tree is drawn (this prevents trying
        # to drawing a tree that is in the process of being modified - probably not a
        # good idea)
        self.tree_mutex.acquire()
        
        # Do a postorder traversal to gather information
        self.tree_width = float(self.tree.getNTips() - 2)
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
            height_thru_nd = self.getEdgeLen(nd)
            if not nd.isTip():
                height_thru_nd += nd.getY()
            if height_thru_nd > parent.getY():
                parent.setY(height_thru_nd)

            nd = nd.getNextPostorder()

        # Do a preorder traversal to draw the tree
        root = self.tree.getFirstPreorder()
        nd = nd.getNextPreorder()

        # Compute scaling factors for x and y axis based on tree height and
        # the width and height of the plot area
        self.tree_height = nd.getY() + self.getEdgeLen(nd)
        self.yscaler = self.usableh/self.tree_height
        self.xscaler = self.usablew/self.tree_width

        # Draw the root node
        x = root.getLeftChild().getX()
        y = 0.0
        id = self.use_node_names and root.getNodeName() or root.getNodeNumber()
        color = root.isSelected() and color_selected_node or color_unselected_node
        self.plotNode(x, y, root.getNodeNumber(), id, color)

        while nd:
            parent = nd.getParent()
            x0 = parent.getX()
            y0 = self.tree_height - parent.getY()
            nd.setY(parent.getY() - self.getEdgeLen(nd))
            x = nd.getX()
            y = self.tree_height - nd.getY()
            color = nd.isSelected() and color_selected_edge or color_unselected_edge
            par_color, fil_color = self.checkCLAstatus(nd)
            self.plotEdge(parent_x=x0, parent_y=y0, child_x=x, child_y=y, radius=2*self.font_Mwidth,
                          edge_color=color, parental_color=par_color, filial_color=fil_color)
            id = self.use_node_names and nd.getNodeName() or nd.getNodeNumber()
            color = nd.isSelected() and color_selected_node or color_unselected_node
            self.plotNode(x, y, nd.getNodeNumber(), id, color)
            nd = nd.getNextPreorder()

        # Release the lock so other threads can play with the tree            
        self.tree_mutex.release()

    def checkCLAstatus(self, nd):
        parental_color = color_undefined_cla
        filial_color = color_undefined_cla
        if nd.isTip():
            td = nd.getTipData()
            if td:
                parental_color = td.parentalCLAValid() and color_valid_cla or color_invalid_cla
        else:
            id = nd.getInternalData()
            if id:
                parental_color = id.parentalCLAValid() and color_valid_cla or color_invalid_cla
                filial_color = id.filialCLAValid() and color_valid_cla or color_invalid_cla
        return parental_color, filial_color
        
    def __init__(self, parent, tree, tree_lock, width, height):
        self.tree_mutex = tree_lock
        Canvas.__init__(self, master=parent, bg=color_plot_background, width=width, height=height)

        self.tree = tree
        self.tree_modified = True
        self.tree_height = 0.0
        self.tree_width = 0.0
        self.xscaler = 1.0
        self.yscaler = 1.0
        self.plot_margin = 20        

        # variables associated with the tree being displayed
        #self.default_tree_topology = '(A,(((B,C)U,(D,(E,F)W)V)T,G)S,(H,(I,(J,K)Z,L)Y)X)R'
        #self.tree = Phylogeny.Tree()
        #self.tree.buildFromString(self.default_tree_topology)
        self.likelihood_root_nodenum = 0

        # variables associated with showing nodes
        self.use_node_names = False

        # variables associated with showing status of conditional likelihood arrays
        self.CLA_radius = 2             # radius of circle plotted for each CLA
        
        # font-related
        self.font = tkFont.Font(family='Courier', size=8)
        fm = self.font.metrics()
        self.font_height = fm['ascent'] + fm['descent']
        self.font_Mwidth = self.font.measure('M')
        self.nodename_text_color = color_unselected_node

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
        Canvas.create_rectangle(self, 0, 0, self.plotw, self.ploth, fill=color_plot_background);
        self.drawTree()

    def reset(self):
        import gc       # garbage collector
        gc.collect()    # should capture return value, which is number of unreachable objects found
        self.repaint()

class TreeViewer(Frame,Thread):
    def __init__(self, tree, mutex, parent=None):
        Thread.__init__(self)
        Frame.__init__(self, parent)
        self.pack(expand=YES, fill=BOTH)

        # variables related to how much is displayed
        self.show_CLAs = False

        # create a frame to hold the menu buttons
        menuf = Frame(self)
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
        self.samplemb.menu.add_command(label='Toggle node numbers', command=self.toggleNodeNamesNumbers)
        self.samplemb.menu.add_command(label='Toggle CLAs', command=self.toggleCLAs)

        # create the canvas
        canvasw = int(0.67*self.winfo_screenwidth())
        canvash = int(0.67*self.winfo_screenheight())
        self.plotter = TreeCanvas(parent=self, tree=tree, tree_lock=mutex, width=canvasw, height=canvash)
        self.plotter.pack(side=TOP, expand=YES, fill=BOTH)
        
        # create the status label
        self.status_label = Label(self, justify=LEFT, relief=SUNKEN, height=1, anchor=W, text='Ready')
        self.status_label.pack(side=TOP, expand=NO, fill=X)
        
        # bind some keys to the application (doesn't matter which widget has focus when you use bind_all
        self.bind_all("<KeyPress-n>", self.keybdToggleNodeNamesNumbers)
        #self.bind_all("<Shift-KeyPress-N>", self.keybdManySteps)

        # configure event is bound only to the main frame
        #self.bind("<Configure>", self.resizing)
        
    def keybdToggleNodeNamesNumbers(self, event):
        self.toggleNodeNamesNumbers()

    def toggleNodeNamesNumbers(self):
        if self.plotter.use_node_names:
            self.plotter.use_node_names = False
            self.status_label.config(text='now showing node numbers')
        else:
            self.plotter.use_node_names = True
            self.status_label.config(text='now showing node names')
        self.plotter.repaint()
        
    def toggleCLAs(self):
        if self.show_CLAs:
            self.show_CLAs = False
            self.status_label.config(text='CLAs hidden')
        else:
            self.show_CLAs = True
            self.status_label.config(text='CLAs visible')

    def close(self):
        self.quit()
        
    def refresh(self, message):
        self.status_label.config(text=message)
        self.plotter.repaint()
        
    def run(self):
        mainloop()

if __name__ == '__main__':
    newick = '(A,(((B,C)U,(D,(E,F)W)V)T,G)S,(H,(I,(J,K)Z,L)Y)X)R'
    t = Phylogeny.Tree()
    t.buildFromString(newick)

    m = Lock()
    tv = TreeViewer(tree=t, mutex=m)

    # Call start() method of base class Thread
    # This invokes TreeViewer.run() in its own thread
    tv.start()
    
