from phycas.Phycas.PhycasCommand import *
from phycas.Phycas.SumTImpl import TreeSummarizer
class SumT(PhycasCommand):
    def __init__(self, p):
        args = (   ("outgroup_taxon",      None,           "Set to the taxon name of the tip serving as the outgroup for display rooting purposes (note: at this time outgroup can consist of just one taxon)"),
                   ("trees",               None,           "A source of trees (list of trees or to the name of the input tree file) to be summarized. This setting should not be None at the time the sumt method is called.", TreeSourceValidate),
                   ("trees_prefix",        'sumt_trees',   "The output tree file in which all distinct tree topologies are saved along with the majority-rule consensus tree will be named <sumt_trees_prefix>.tre and the corresponding pdf file containing graphical representations of these trees will be named <sumt_trees_prefix>.pdf. This setting cannot be None when the sumt method is called."),
                   ("splits_prefix",       'sumt_splits',  "The pdf file showing plots depicting split posteriors through time and split sojourns will be named <sumt_splits_prefix>.pdf. If None, this analysis will be skipped."),
                   ("output_replace",      False,          "If True, output files will be replaced automatically if they exist; if False, a random integer will be added to the name so that the name no longer matches an existing file", BoolArgValidate),
                   ("burnin",              1,              "Number of trees from the input list of trees to skip", IntArgValidate(min=0)),
                   ("equal_brlens",        False,          "If True, trees in pdf file will be drawn with branch lengths equal, making support values easier to see; if set to True, consider setting pdf_scalebar_position = None (scalebar is irrelevant in this case)", BoolArgValidate),
                   ("tree_credible_prob",  0.95,           "Include just enough trees in the <sumt_trees_prefix>.tre and <sumt_trees_prefix>.pdf files such that the cumulative posterior probability is greater than this value"),
                   ("rooted",              False,          "Set to True if trees are rooted; otherwise, leave set to default value of False to assume trees are unrooted", BoolArgValidate),
                )
        PhycasCommand.__init__(self, p, args)

    def __call__(self, **kwargs):
        self.set(**kwargs)
        tree_summarizer = TreeSummarizer(self.phycas, self)
        tree_summarizer.consensus()
        
