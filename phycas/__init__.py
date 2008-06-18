'''Root of the phycas package.'''
#__all__ = [
#    'Conversions',
#    'DataMatrix', 
#    'Likelihood', 
#    'PDFGen',
#    'Phylogeny',
#    'ProbDist', 
#    'ReadNexus', 
#    'PhycasA',
#    'Examples',
#    ]

import Conversions
import DataMatrix
import Likelihood
import PDFGen
import Phylogeny
import ProbDist
import ReadNexus
from Phycas import Phycas
_user_ini_checked = False
if not _user_ini_checked:
    import os
    _user_ini_checked = True
    p = os.path.expanduser("~/.phycas/startup.py")
    if os.path.exists(p):
        execfile(p)
phycas = Phycas()


import textwrap


###############################################################################
#   Taken from
#   http://mail.python.org/pipermail/python-list/2006-February/365594.html

def ttysize():
    try:
        import fcntl
        import struct
        import termios
        buf = 'abcdefgh'
        buf = fcntl.ioctl(0, termios.TIOCGWINSZ, buf)
        row, col, rpx, cpx = struct.unpack('hhhh', buf)
        return row, col
    except:
        return None
###############################################################################


class PhycasCmdOpts(object):
    _help_str_wrapper = textwrap.TextWrapper()
    _help_fmt_str = "%-19s %-19s %s"
    
    def __init__(self, command, args):
        self.help_info = {}
        self.optionsInOrder = args
        for opt in args:
            name, default, help_str = opt
            command.__dict__[name] = default
            self.help_info[name] =  [default, help_str]
    def __str__(self):
        PhycasCmdOpts._reset_term_width()
        opts_help = []
        for i in self.optionsInOrder:
            hs = "\n".join(PhycasCmdOpts._help_str_wrapper.wrap(i[2]))
            s = PhycasCmdOpts._help_fmt_str % (i[0], i[1], hs[40:])
            opts_help.append(s)
        return "\n".join(opts_help)

    def _reset_term_width():
        tty_sz = ttysize()
        if tty_sz:
            PhycasCmdOpts._set_terminal_width(tty_sz[1])
    _reset_term_width = staticmethod(_reset_term_width)
    def _set_terminal_width(w):
        new_width = max(60, w)
        PhycasCmdOpts._help_str_wrapper.width = new_width
        PhycasCmdOpts._help_str_wrapper.subsequent_indent = " "*40
        PhycasCmdOpts._help_str_wrapper.initial_indent = " "*40
    _set_terminal_width = staticmethod(_set_terminal_width)
PhycasCmdOpts._set_terminal_width(80)

class Sumt:
    def __init__(self, p):
        args = (   ("outgroup_taxon",      None,           "Set to the taxon name of the tip serving as the outgroup for display rooting purposes (note: at this time outgroup can consist of just one taxon)"),
                   ("input_tree_file",     None,           "Set to the name of the input tree file. This setting should not be None at the time the sumt method is called."),
                   ("trees_prefix",        'sumt_trees',   "The output tree file in which all distinct tree topologies are saved along with the majority-rule consensus tree will be named <sumt_trees_prefix>.tre and the corresponding pdf file containing graphical representations of these trees will be named <sumt_trees_prefix>.pdf. This setting cannot be None when the sumt method is called."),
                   ("splits_prefix",       'sumt_splits',  "The pdf file showing plots depicting split posteriors through time and split sojourns will be named <sumt_splits_prefix>.pdf. If None, this analysis will be skipped."),
                   ("output_replace",      False,          "If True, output files will be replaced automatically if they exist; if False, a random integer will be added to the name so that the name no longer matches an existing file"),
                   ("burnin",              1,              "Number of trees to skip in sumt_input_tree_file"),
                   ("equal_brlens",        False,          "If True, trees in pdf file will be drawn with branch lengths equal, making support values easier to see; if set to True, consider setting pdf_scalebar_position = None (scalebar is irrelevant in this case)"),
                   ("tree_credible_prob",  0.95,           "Include just enough trees in the <sumt_trees_prefix>.tre and <sumt_trees_prefix>.pdf files such that the cumulative posterior probability is greater than this value"),
                   ("rooted",              False,          "Set to True if trees in sumt_input_tree_file are rooted; otherwise, leave set to default value of False to assume trees are unrooted"),
                )
        self.options = PhycasCmdOpts(self, args)
        self.phycas = p

    def __call__(self):
        self.phycas.sumt_outgroup_taxon     = self.outgroup_taxon
        self.phycas.sumt_input_tree_file    = self.input_tree_file   
        self.phycas.sumt_trees_prefix       = self.trees_prefix      
        self.phycas.sumt_splits_prefix      = self.splits_prefix     
        self.phycas.sumt_output_replace     = self.output_replace    
        self.phycas.sumt_burnin             = self.burnin            
        self.phycas.sumt_equal_brlens       = self.equal_brlens      
        self.phycas.sumt_tree_credible_prob = self.tree_credible_prob
        self.phycas.sumt_rooted             = self.rooted            
        self.phycas.sumt()
        
    def help(self, opt=None):
       print self.options
       
sumt = Sumt(phycas)

#print 'importing phycas...'
