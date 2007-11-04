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

from phycas import *
#from threading import *
from Tkinter import *
from tkMessageBox import askokcancel, askquestion, showinfo
from tkSimpleDialog import askstring
from tkFileDialog import askopenfilename
import tkFont
import math

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
color_plot_background       = midnight
color_undefined_cla         = silver
color_valid_cla             = green
color_valid_cached_dot      = dkgreen
color_invalid_cla           = red
color_invalid_cached_dot    = maroon
color_selected_edge         = white
color_unselected_edge       = silver
color_selected_node         = white
color_unselected_node       = silver
color_likelihood_root       = magenta

# The text displayed in the Help | About dialog box
helptext = """
Knobs on ends of edges represent conditional likelihood arrays (CLAs)
Valid CLAs are GREEN
Invalid CLAs are RED
Undefined CLAs are GRAY
Dotted CLAs are cached
Likelihood root node is PINK

Keyboard shortcuts:
  h - opens this help dialog box
  n - toggles between node numbers and names
  q - quits application

Currently, edge are not shown proportional to their lengths.
"""

class TreeCanvas(Canvas):
    # Wrapping Canvas calls in these display* functions in order to later make
    # it easier to draw to a PDF file rather than the screen
    def displayText(self, x, y, text, font, color):
        Canvas.create_text(self, x, y, text=text, font=font, fill=color)
        
    def displayLine(self, x0, y0, x, y, color, thickness):
        Canvas.create_line(self, x0, y0, x, y, fill=color, width=thickness)

    def displayFilledOval(self, x0, y0, x, y, color):
        Canvas.create_oval(self, x0, y0, x, y, fill=color, outline=color)

    def displayFilledRectangle(self, x, y, width, height, color):
        Canvas.create_rectangle(self, x, y, width, height, fill=color)
        
    def xtranslate(self, x):
        return int(self.left + self.xscaler*x)

    def ytranslate(self, y):
        return int(self.bottom - self.yscaler*y)

    def plotNode(self, xval, yval, number, label, node_color=color_unselected_node):
        x = self.xtranslate(xval)
        y = self.ytranslate(yval)
        color = (number == self.likelihood_root_nodenum and color_likelihood_root or node_color)
        self.displayText(x, y, text=str(label), font=self.font, color=color)

    def plotEdge(self, parent_x, parent_y, child_x, child_y,
                 thickness=1,
                 nodenum_radius=10, cla_radius=3, cached_radius=1,
                 edge_color=color_unselected_edge,
                 parental_color=color_undefined_cla,
                 parental_cached_color=color_undefined_cla,
                 filial_color=color_undefined_cla,
                 filial_cached_color=color_undefined_cla):
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
        dx = float(nodenum_radius)*math.cos(theta)
        dy = float(nodenum_radius)*math.sin(theta)
        
        # Draw the edge itself
        self.displayLine(x0+dx, y0+dy, x-dx, y-dy, edge_color, thickness)

        # Draw the parental CLA marker
        self.displayFilledOval(x0+dx-cla_radius, y0+dy-cla_radius, x0+dx+cla_radius, y0+dy+cla_radius, color=parental_color)
        self.displayFilledOval(x0+dx-cached_radius, y0+dy-cached_radius, x0+dx+cached_radius, y0+dy+cached_radius, color=parental_cached_color)

        # Draw the filial CLA marker
        self.displayFilledOval(x-dx-cla_radius, y-dy-cla_radius, x-dx+cla_radius, y-dy+cla_radius, color=filial_color)
        self.displayFilledOval(x-dx-cached_radius, y-dy-cached_radius, x-dx+cached_radius, y-dy+cached_radius, color=filial_cached_color)

    def getEdgeLen(self, nd):
        return (self.use_edgelens and nd.getEdgeLen() or 1.0)

    def drawTree(self):
        # The angle at which edges are drawn is theta. If v is the length of an
        # edge, then dx = v cos(theta) and dy = v sin(theta).
        #
        # X         Y
        # \        /      |    This is a drawing of one pair of tips (X, Y)
        #  \      /       |    and their parent (Z). The length of the
        #   \   v/       dy    right child's edge is v. The distance between
        #    \  / theta   |    sibling distance (the intersibling distance, 
        #     \/____      |    or isd) is the horizontal distance between
        #     Z  dx            X and Y. The isd determines theta.
        #
        # The plot area is w pixels wide and h pixels high, so we need to find
        # the angle that would allow us to draw a tree that has the same ratio
        # of w to h. The width of the tree when plotted using a given theta is
        # equal to (left_sum + right_sum) * cos(theta), where left_sum is the
        # sum of all edge lengths from the subroot (only child of root) to the
        # subroot's leftmost descendant, and right_sum is the sum of all edge
        # lengths from the subroot node to its rightmost descendant (which will
        # always be the tree's first postorder node). Likewise, the height of
        # the tree when plotted will be the subroot's edge length plus
        # longest_path * sin(theta), where longest_path is the sum of the edge
        # lengths in the path from the subroot to the furthest leaf.

        # acquire will block other threads while tree is drawn (this prevents trying
        # to drawing a tree that is in the process of being modified - probably not a
        # good idea)
        #if self.tree_mutex:
        #    self.tree_mutex.acquire()
        
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
            par_color, par_dot_color, fil_color, fil_dot_color = self.checkCLAstatus(nd)
            self.plotEdge(parent_x=x0, parent_y=y0, child_x=x, child_y=y,
                          thickness=nd.isSelected() and 3 or 1,
                          nodenum_radius=2*self.font_Mwidth,
                          cla_radius=4, cached_radius=2, 
                          edge_color=color,
                          parental_color=par_color, parental_cached_color=par_dot_color,
                          filial_color=fil_color, filial_cached_color=fil_dot_color)
            id = self.use_node_names and nd.getNodeName() or nd.getNodeNumber()
            color = nd.isSelected() and color_selected_node or color_unselected_node
            self.plotNode(x, y, nd.getNodeNumber(), id, color)
            nd = nd.getNextPreorder()

        # Release the lock so other threads can play with the tree            
        #if self.tree_mutex:
        #    self.tree_mutex.release()

    def checkCLAstatus(self, nd):
        parental_color = color_undefined_cla
        filial_color = color_undefined_cla
        parental_cached_color = color_undefined_cla
        filial_cached_color = color_undefined_cla
        if nd.isTip():
            td = nd.getTipData()
            if td:
                if td.parentalCLAValid():
                    parental_color = color_valid_cla
                    if td.parentalCLACached():
                        parental_cached_color = color_valid_cached_dot
                    else:
                        parental_cached_color = color_valid_cla
                else:
                    parental_color = color_invalid_cla
                    if td.parentalCLACached():
                        parental_cached_color = color_invalid_cached_dot
                    else:
                        parental_cached_color = color_invalid_cla
        else:
            id = nd.getInternalData()
            if id:
                if id.parentalCLAValid():
                    parental_color = color_valid_cla
                    if id.parentalCLACached():
                        parental_cached_color = color_valid_cached_dot
                    else:
                        parental_cached_color = color_valid_cla
                else:
                    parental_color = color_invalid_cla
                    if id.parentalCLACached():
                        parental_cached_color = color_invalid_cached_dot
                    else:
                        parental_cached_color = color_invalid_cla
                if id.filialCLAValid():
                    filial_color = color_valid_cla
                    if id.filialCLACached():
                        filial_cached_color = color_valid_cached_dot
                    else:
                        filial_cached_color = color_valid_cla
                else:
                    filial_color = color_invalid_cla
                    if id.filialCLACached():
                        filial_cached_color = color_invalid_cached_dot
                    else:
                        filial_cached_color = color_invalid_cla
        return parental_color, parental_cached_color, filial_color, filial_cached_color
        
    #def __init__(self, parent, tree, tree_lock, width, height):
    def __init__(self, parent, tree, width, height):
        #self.tree_mutex = tree_lock
        Canvas.__init__(self, master=parent, bg=color_plot_background, width=width, height=height)

        self.tree = tree
        self.use_edgelens = False   # by default, all edges are drawn as equal in length (call useEdgelens() to change this)
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
        self.displayFilledRectangle(0, 0, self.plotw, self.ploth, color=color_plot_background)
        self.drawTree()

    def reset(self):
        import gc       # garbage collector
        gc.collect()    # should capture return value, which is number of unreachable objects found
        self.repaint()

    def useEdgelens(self):
        # Can only honor this request if tree has edge lengths
        self.use_edgelens = self.tree.hasEdgeLens()
        if self.use_edgelens:
            self.repaint()

#class TreeViewer(Frame,Thread):
class TreeViewer(Frame):
    #def __init__(self, tree, mutex=None, parent=None):
    def __init__(self, tree, msg='Tree Viewer for Debugging', parent=None):
        #Thread.__init__(self)
        Frame.__init__(self, parent)
        self.pack(expand=YES, fill=BOTH)

        # set the window title
        self.winfo_toplevel().title(msg)

        # always position the main window 50 pixels from top and 50 pixels from left
        self.winfo_toplevel().geometry("+%d+%d" % (50, 50))

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
        self.samplemb.pack(expand=NO, fill=X, side=LEFT)
        self.samplemb.menu = Menu(self.samplemb, tearoff=0)
        self.samplemb['menu'] = self.samplemb.menu
        self.samplemb.menu.add_command(label='Toggle node numbers', command=self.toggleNodeNamesNumbers)

        # create the Help menu button
        self.helpmb = Menubutton(menuf, text='Help', relief=RAISED, anchor=W, borderwidth=0)
        self.helpmb.pack(expand=YES, fill=X, side=LEFT)
        self.helpmb.menu = Menu(self.helpmb, tearoff=0)
        self.helpmb['menu'] = self.helpmb.menu
        self.helpmb.menu.add_command(label='About', command=self.helpAbout)
        
        # create the canvas
        canvasw = int(0.67*self.winfo_screenwidth())
        canvash = int(0.67*self.winfo_screenheight())
        #self.plotter = TreeCanvas(parent=self, tree=tree, tree_lock=mutex, width=canvasw, height=canvash)
        self.plotter = TreeCanvas(parent=self, tree=tree, width=canvasw, height=canvash)
        self.plotter.pack(side=TOP, expand=YES, fill=BOTH)
        
        # create the status label
        self.status_label = Label(self, justify=LEFT, relief=SUNKEN, height=1, anchor=W, text='Ready')
        self.status_label.pack(side=TOP, expand=NO, fill=X)
        
        # bind some keys to the application (doesn't matter which widget has focus when you use bind_all
        self.bind_all("<KeyPress-n>", self.keybdToggleNodeNamesNumbers)
        self.bind_all("<KeyPress-h>", self.keybdHelpAbout)
        self.bind_all("<KeyPress-q>", self.keybdQuit)
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
        
    def keybdHelpAbout(self, event):
        self.helpAbout()

    def helpAbout(self):
        showinfo('About the Phycas TreeViewer', helptext)

    def keybdQuit(self, event):
        self.close()

    def close(self):
        self.quit()
        
    def refresh(self, message):
        self.status_label.config(text=message)
        self.plotter.repaint()
        
    def setLikelihoodRoot(self, like_root_node_num):
        if like_root_node_num < 0:
            self.plotter.likelihood_root_nodenum = None
        else:
            self.plotter.likelihood_root_nodenum = like_root_node_num
        
    def run(self):
        mainloop()

if __name__ == '__main__':
    newick = '(A,(((B,C)U,(D,(E,F)W)V)T,G)S,(H,(I,(J,K)Z,L)Y)X)R'
    t = Phylogeny.Tree()
    t.buildFromString(newick)

    #tv = TreeViewer(tree=t)
    TreeViewer(tree=t, msg='Showing default tree').run()

    # Call start() method of base class Thread
    # This invokes TreeViewer.run() in its own thread
    #tv.start()
    
