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
import sys
from Phycas import Phycas
from Phycas.PhycasCommand import PhycasCmdOpts, PhycasCommand, BoolArgValidate, IntArgValidate
_user_ini_checked = False

if not _user_ini_checked:
    import os
    _user_ini_checked = True
    p = os.path.expanduser("~/.phycas/startup.py")
    if os.path.exists(p):
        execfile(p)

phycas = Phycas()
def error_msg(msg):
    sys.stderr.write("Error: %s\n" % msg)

def phycas_except_hook(t, v, tb):
    error_msg(v)

sys.excepthook = phycas_except_hook

class Sumt(PhycasCommand):
    def __init__(self, p):
        args = (   ("outgroup_taxon",      None,           "Set to the taxon name of the tip serving as the outgroup for display rooting purposes (note: at this time outgroup can consist of just one taxon)"),
                   ("input_tree_file",     None,           "Set to the name of the input tree file. This setting should not be None at the time the sumt method is called."),
                   ("trees_prefix",        'sumt_trees',   "The output tree file in which all distinct tree topologies are saved along with the majority-rule consensus tree will be named <sumt_trees_prefix>.tre and the corresponding pdf file containing graphical representations of these trees will be named <sumt_trees_prefix>.pdf. This setting cannot be None when the sumt method is called."),
                   ("splits_prefix",       'sumt_splits',  "The pdf file showing plots depicting split posteriors through time and split sojourns will be named <sumt_splits_prefix>.pdf. If None, this analysis will be skipped."),
                   ("output_replace",      False,          "If True, output files will be replaced automatically if they exist; if False, a random integer will be added to the name so that the name no longer matches an existing file", BoolArgValidate),
                   ("burnin",              1,              "Number of trees to skip in sumt_input_tree_file", IntArgValidate),
                   ("equal_brlens",        False,          "If True, trees in pdf file will be drawn with branch lengths equal, making support values easier to see; if set to True, consider setting pdf_scalebar_position = None (scalebar is irrelevant in this case)", BoolArgValidate),
                   ("tree_credible_prob",  0.95,           "Include just enough trees in the <sumt_trees_prefix>.tre and <sumt_trees_prefix>.pdf files such that the cumulative posterior probability is greater than this value"),
                   ("rooted",              False,          "Set to True if trees in sumt_input_tree_file are rooted; otherwise, leave set to default value of False to assume trees are unrooted", BoolArgValidate),
                )
        PhycasCommand.__init__(self, p, args)

    def __call__(self, **kwargs):
        for key, value in kwargs.iteritems():
            setattr(self, key, value)
        import SumTImpl
        tree_summarizer = SumTImpl.TreeSummarizer(self)
        tree_summarizer.consensus()
       
sumt = Sumt(phycas)

#print 'importing phycas...'
