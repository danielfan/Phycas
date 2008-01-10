import os,sys,math
from phycas.TreeViewer import *
from phycas import *

class TreeSummarizer(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Saves consensus tree and QQ plots in pdf files.
    
    """
    def __init__(self, phycas):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes TreeSummarizer object.
        
        """
        self.phycas = phycas

    def readTreesFromFile(self):
        self.phycas.reader.readFile(self.phycas.sumt_tfile_name)
        self.taxon_labels = self.phycas.reader.getTaxLabels()
        self.stored_treenames = []
        self.stored_newicks = []
        for t in self.phycas.reader.getTrees():
            self.stored_treenames.append(t.name) # should use hash to store both names and newicks
            self.stored_newicks.append(t.newick)
        self.phycas.phycassert(len(self.stored_newicks) > 0, 'expecting a trees block defining at least one tree in the file %s' % self.phycas.sumt_tfile_name)

    def consensus(self):
        # Open pdf_treefile and read trees therein
        self.readTreesFromFile()

        # Build each tree and determine its height
        t = Phylogeny.Tree()
        max_height = 0.0
        num_trees = 0
        num_trees_considered = 0
        num_new_splits = 0
        num_old_splits = 0
        split_map = {}
        first_ss = None
        for newick in self.stored_newicks:
            if num_trees < self.phycas.sumt_burnin:
                num_trees += 1
                continue
            num_trees += 1
            num_trees_considered += 1
            t.buildFromString(newick, True)
            t.rectifyNames(self.taxon_labels)
            if self.phycas.sumt_outgroup_taxon:
                num = t.findTipByName(self.phycas.sumt_outgroup_taxon)
                self.phycas.phycassert(num, 'could not root tree using specified outgroup: no tip having name "%s" could be found' % self.phycas.sumt_outgroup_taxon)
                t.rerootAtTip(num)
            t.recalcAllSplits(t.getNTips())
            nd = t.getFirstPreorder()
            assert nd.isRoot(), 'first preorder node should be the root'
            split_vec = []
            while True:
                nd = nd.getNextPreorder()
                if not nd:
                    break
                elif nd.isInternal():
                    s = nd.getSplit()
                    if s.isBitSet(0):
                        print 'Inverting split...'
                        print '  before =',s.createPatternRepresentation()
                        s.invertSplit()
                        print '  after =',s.createPatternRepresentation()
                    ss = s.createPatternRepresentation()
                    print 'nd number =',nd.getNodeNumber()
                    print 'ss =',ss
                    t.unselectAllNodes()
                    nd.selectNode()
                    TreeViewer.TreeViewer(tree=t, msg='Current tree').run()
                    raw_input('press enter to continue...')
                    if ss in split_vec:
                        raw_input('bad')
                    else:
                        split_vec.append(ss);
                    if ss in split_map.keys():
                        split_map[ss] += 1
                        num_old_splits += 1
                    else:
                        split_map[ss] = 1
                        num_new_splits += 1
                        if num_new_splits == 1:
                            first_ss = ss

        print 'Number of trees considered =',num_trees_considered        
        print 'Number of splits in map =',len(split_map.keys())                
        print 'Number of new splits =',num_new_splits            
        print 'Number of old splits =',num_old_splits

        for k in split_map.keys():
            p = float(split_map[k])/float(num_trees_considered)
            if p > 0.5:
                print '%s => %d' % (k,split_map[k])

        # Build each tree again and save in PDF file            
        #pdf = PDFGenerator(self.pdf_page_width, self.pdf_page_height)
        #pdf.overwrite = True
        #for newick in self.stored_newicks:
        #    tree.buildFromString(newick, True)
        #    tree.rectifyNames(self.taxon_labels)
        #    if self.pdf_outgroup_taxon:
        #        num = tree.findTipByName(self.pdf_outgroup_taxon)
        #        tree.rerootAtTip(num)
        #    if self.pdf_ladderize:
        #        if self.pdf_ladderize == 'right':
        #            tree.ladderizeRight()
        #        else:
        #            tree.ladderizeLeft()
        #    tree.rectifyNames(self.taxon_labels)
        #    pdf.newPage()
        #    self.tree2pdf(pdf, tree, max_height)
        #pdf.saveDocument(self.phycas.sumt_contree_pdf_file)
