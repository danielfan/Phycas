import os,sys,math,random
#from phycas.PDFGen import *
from phycas.TreeViewer import *
from phycas import *

class Brownian(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Used to compute variance-covariance matrix for rooted trees assuming
    a Brownian-motion model in which variances are proportional to branch
    length.
    
    """
    def __init__(self, phycas):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes Brownian object by assigning supplied phycas object to a
        data member variable.
        
        """
        self.phycas = phycas

    def readTreesFromFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reads all trees from the file named in self.phycas.brownian_input_tree_file.
        Stores tree descriptions in self.stored_tree_defs and taxon labels in
        self.taxon_labels.
        
        """
        self.phycas.reader.readFile(self.phycas.brownian_input_tree_file)
        self.taxon_labels = self.phycas.reader.getTaxLabels()
        self.self.stored_tree_defs = self.phycas.reader.getTrees()
        self.phycas.phycassert(len(self.stored_tree_defs) > 0, 'expecting a trees block defining at least one tree in the file %s' % self.phycas.sumt_input_tree_file)

    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes and displays the variance-covariance matrix for each tree in
        brownian_input_tree_file. Variances on the main diagonal are distances
        from the root to each tip, and covariances off the diagonal are
        distances from the root to the most recent common ancestor of the pair
        of taxa represented by the row and column.
        
        """
        # Check to make sure user specified an input tree file
        self.phycas.phycassert(self.phycas.brownian_input_tree_file, 'brownian_input_tree_file must be specified before varcov method is called')
            
        num_trees_considered = 0

        # Open brownian_input_tree_file and read trees therein
        self.phycas.output('\nReading trees from file %s...' % self.phycas.brownian_input_tree_file)
        self.readTreesFromFile()

        self.phycas.output('Key to taxa:')
        for i,taxlabel in enumerate(self.taxon_labels):
            self.phycas.output('%4d %s' % (i+1,taxlabel))

        ntax = len(self.taxon_labels)            

        # Build each tree and add the splits and tree topolgies found there to the
        # dictionary of splits (split_map) and the dictionary of tree topologies
        # (tree_map), respectively
        #self.phycas.output('\nEvaluating trees...')
        t = Phylogeny.Tree()
        t.setRooted()

        # values used for display purposes
        num_stored_trees = len(self.stored_tree_defs)
        self.phycas.phycassert(num_stored_trees > 0, 'Specified tree file (%s) contained no stored trees' % self.phycas.sumt_input_tree_file)
        
        for tree_def in self.stored_tree_defs:
            num_trees_considered += 1

            # Create variance-covariance matrix
            matrix = []
            for i in range(ntax):
                matrix.append([0.0]*ntax)
            
            # Build the tree
            #raw_input('stopped')
            tree_def.buildTree(t)
            t.rectifyNames(self.taxon_labels)
            ntips = t.getNTips()
            t.recalcAllSplits(ntips)

            #tv = TreeViewer.TreeViewer(t)
            #tv.run()

            nd = t.getFirstPreorder()
            assert nd.isRoot(), 'the first preorder node should be the root'
            while True:
                nd = nd.getNextPreorder()
                if not nd:
                    break
                else:
                    # Grab the edge length
                    edge_len = nd.getEdgeLen()
                    
                    # Grab the split and invert it if necessary to attain a standard polarity
                    s = nd.getSplit()
                        
                    # Create a string representation of the split
                    ss = s.createPatternRepresentation()

                    on = s.getOnList()
                    num_on = len(on)
                    onstr = ','.join([('%d' % i) for i in on])

                    #self.phycas.output('%12.5f %s %d %2s %s' % (edge_len, ss, nd.getNodeNumber(), nd.getNodeName(), onstr))
                    self.phycas.phycassert(num_on > 0, 'error: split found with no bits set in Brownian.run function')
                    if num_on == 1:
                        v = on[0]
                        matrix[v][v] += edge_len
                    else:
                        for v in on:
                            matrix[v][v] += edge_len
                        for i in range(num_on - 1):
                            v = on[i]
                            for j in range(i + 1, num_on):
                                w = on[j]
                                matrix[v][w] += edge_len
                                matrix[w][v] += edge_len

            self.phycas.output('\nVar-Cov matrix for tree %d\n' % num_trees_considered)
            tmpstr = '%12s' % ' '
            for j in range(ntax):
                tmpstr += ' %12d' % (j+1)
            self.phycas.output('%s' % tmpstr)
            for i in range(ntax):
                tmpstr = '%12d' % (i+1)
                for j in range(ntax):
                    tmpstr += ' %12.5f' % matrix[i][j]
                self.phycas.output('%s' % tmpstr)
                    
        self.phycas.output('\nBrownian finished.')
