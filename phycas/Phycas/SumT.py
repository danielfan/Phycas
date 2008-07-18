from phycas.Utilities.PhycasCommand import *
from phycas.Phycas.SumTImpl import TreeSummarizer
class SumT(PhycasCommand):
    def __init__(self, p):
        args = (   ("outgroup_taxon",      None,           "Set to the taxon name of the tip serving as the outgroup for display rooting purposes (note: at this time outgroup can consist of just one taxon)"),
                   ("trees",               None,           "A source of trees (list of trees or to the name of the input tree file) to be summarized. This setting should not be None at the time the sumt method is called.", TreeSourceValidate),
                   ("burnin",              1,              "Number of trees from the input list of trees to skip", IntArgValidate(min=0)),
                   ("tree_credible_prob",  0.95,           "Include just enough trees in the <sumt_trees_prefix>.tre and <sumt_trees_prefix>.pdf files such that the cumulative posterior probability is greater than this value", ProbArgValidate()),
                   ("useGUI",              True,           "If True, and if wxPython is installed, a graphical user interface (GUI) will be used to display, and allow manipulation of, AWTY plots", BoolArgValidate),
                   ("rooted",              False,          "Set to True if trees are rooted; otherwise, leave set to default value of False to assume trees are unrooted", BoolArgValidate),
                )
        o = PhycasCommandOutputOptions()
        o.__dict__["_help_order"] = ["trees", "splits"]
        t = TreeOutputSpec('sumt_trees', "The output tree file in which all distinct tree topologies are saved along with the majority-rule consensus tree will be named <sumt_trees_prefix>.tre and the corresponding pdf file containing graphical representations of these trees will be named <sumt_trees_prefix>.pdf. This setting cannot be None when the sumt method is called.")
        t._options = PhycasCmdOpts(t, (("equal_brlens", False, "If True, trees in pdf file will be drawn with branch lengths equal, making support values easier to see; if set to True, consider setting pdf_scalebar_position = None (scalebar is irrelevant in this case)", BoolArgValidate),))
        o.__dict__["trees"] = t
        s = PDFOutputSpec('sumt_splits', "The pdf file showing plots depicting split posteriors through time and split sojourns will be named <sumt_splits_prefix>.pdf. If None, this analysis will be skipped.")
        o.__dict__["splits"] = s
        PhycasCommand.__init__(self, p, args, "sumt", "The sumt command is used to summarize a collection of trees (usually trees that have been produced by an MCMC simulation).", o)

    def __call__(self, **kwargs):
        self.set(**kwargs)
        tree_summarizer = TreeSummarizer(self.phycas, self)
        tree_summarizer.consensus()
        
