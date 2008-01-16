import os,sys,math,random
from phycas.PDFGen import *
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
        Initializes TreeSummarizer object by assigning supplied phycas object
        to a data member variable.
        
        """
        self.phycas = phycas
        if not self.phycas.sumt_input_tree_file:
            self.phycas.sumt_input_tree_file = self.phycas.getPrefix()+'.t'
        if not self.phycas.sumt_output_tree_file:
            self.phycas.sumt_output_tree_file = self.phycas.getPrefix()+'.distinct.t'

    def readTreesFromFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reads all trees from the file named in self.phycas.sumt_input_tree_file.
        Stores tree descriptions in self.stored_newicks and taxon labels in
        self.taxon_labels.
        
        """
        self.phycas.reader.readFile(self.phycas.sumt_input_tree_file)
        self.taxon_labels = self.phycas.reader.getTaxLabels()
        self.stored_treenames = []
        self.stored_newicks = []
        for t in self.phycas.reader.getTrees():
            self.stored_treenames.append(t.name) # should use hash to store both names and newicks
            self.stored_newicks.append(t.newick)
        self.phycas.phycassert(len(self.stored_newicks) > 0, 'expecting a trees block defining at least one tree in the file %s' % self.phycas.sumt_input_tree_file)

    def assignEdgeLensAndSupportValues(self, tree, split_map, total_samples):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Assigns edge lengths in supplied tree using posterior mean split
        weights stored in split_vect. The variable split_vect associates a
        list with each split. The first element of this list is the number of
        trees considered in which the split was found, and the second element
        is the sum of the edge lengths of the split over all trees in which it
        was found. The estimated edge length applied to the tree is the
        second element of the split divided by the first element.
        
        """
        tree.recalcAllSplits(tree.getNTips())
        nd = tree.getFirstPreorder()
        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            s = nd.getSplit()
            ss = s.createPatternRepresentation()
            v = split_map[ss]
            assert v, 'could not find edge length information for the split %s' % ss
            support = float(v[0])/float(total_samples)
            nd.setSupport(support)
            if self.phycas.sumt_equal_brlens:
                edge_len = 1.0
            else:
                edge_len = float(v[1])/float(v[0])
            nd.setEdgeLen(edge_len)

    def save_trees(self, tree_file, pdf_file, short_names, full_names, trees):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Saves the supplied tree names (tree_names) and newick descriptions
        (trees) to a nexus tree file named tree_file. If a file by that name
        exists, a new name will be invented by adding a random integer
        between 1 and 1 million to the end of the supplied tree file name.
        
        """
        # Open the output pdf file
        pdf = PDFGenerator(self.phycas.pdf_page_width, self.phycas.pdf_page_height)
        alt_pdf_file = pdf_file
        if not alt_pdf_file:
            alt_pdf_file = 'sumt.pdf'
        if self.phycas.sumt_output_replace:
            pdf.overwrite = True
        else:
            while os.path.exists(alt_pdf_file):
                alt_pdf_file = '%s_%d.pdf' % (pdf_file, random.randint(1,1000000))
            if alt_pdf_file != pdf_file:
                self.phycas.output('Supplied pdf file name %s exists, using %s instead' % (pdf_file, alt_pdf_file))
        
        # Open the output tree file        
        alt_tree_file = tree_file
        if not alt_tree_file:
            alt_tree_file = 'sumt.tre'
        if not self.phycas.sumt_output_replace:
            while os.path.exists(alt_tree_file):
                alt_tree_file = '%s_%d.tre' % (tree_file, random.randint(1,1000000))
            if alt_tree_file != tree_file:
                self.phycas.output('Supplied tree file name %s exists, using %s instead' % (tree_file, alt_tree_file))
        treef = open(alt_tree_file, 'w')
        
        treef.write('#nexus\n\n')
        treef.write('begin trees;\n')
        treef.write('  translate\n')
        ntax = len(self.taxon_labels)
        for i,label in enumerate(self.taxon_labels):
            if ' ' in label:
                tmp = "'%s'" % label
            else:
                tmp = label
            treef.write("  %d %s%s\n" % (i+1, tmp, (i < ntax - 1) and ',' or ''))
        treef.write('  ;\n')
        for short,full,tree in zip(short_names, full_names, trees):
            # Reroot (if requested) tree and output newick description to tree file
            if self.phycas.sumt_outgroup_taxon:
                num = list(self.taxon_labels).index(self.phycas.sumt_outgroup_taxon)
                if num == len(self.taxon_labels):
                    self.phycas.output('Warning: could not root tree using specified outgroup: no tip having name "%s" could be found' % self.phycas.sumt_outgroup_taxon)
                    self.phycas.output('Here is the list of all taxon labels:')
                    for label in self.taxon_labels:
                        self.phycas.output('  %s' % label)
                tree.rerootAtTip(num)
            treef.write('  tree %s = [&U] %s;\n' % (short, tree.makeNewick()))

            # Add taxon names to tip nodes and output to pdf file
            tree.rectifyNames(self.taxon_labels)
            if self.phycas.pdf_ladderize:
                if self.phycas.pdf_ladderize == 'right':
                    tree.ladderizeRight()
                else:
                    tree.ladderizeLeft()
            pdf.newPage()
            self.phycas.tree2pdf(pdf, tree, full, show_support = True)
            
        pdf.saveDocument(pdf_file)
        treef.write('end;\n')
        treef.close()
        
    def consensus(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This is the main member function of the class. It computes and outputs
        a summary of splits found in the sampled trees, as well as the
        distinct tree topologies found. It also computes a majority rule
        consensus tree and saves this to the file sumt_output_tree_file if
        this variable is not set to None.
        
        """
        num_trees = 0
        num_trees_considered = 0
        split_map = {}
        tree_map = {}

        # Open sumt_tfile_name and read trees therein
        self.phycas.output('\nReading trees from file %s...' % self.phycas.sumt_input_tree_file)
        self.readTreesFromFile()

        # Build each tree and add the splits and tree topolgies found there to the
        # dictionary of splits (split_map) and the dictionary of tree topologies
        # (tree_map), respectively
        self.phycas.output('Compiling lists of tree topologies and splits...')
        t = Phylogeny.Tree()
        for newick in self.stored_newicks:
            if num_trees < self.phycas.sumt_burnin:
                num_trees += 1
                continue
            num_trees += 1
            num_trees_considered += 1

            # Create a tree key, which is a list of internal node splits that can be
            # used to uniquely identify a tree topology
            tree_key = []

            # Build the tree
            t.buildFromString(newick, True)
            t.rectifyNames(self.taxon_labels)
            t.recalcAllSplits(t.getNTips())

            # Each split found in any tree is associated with a list by the dictionary split_map
            # The list is organized as follows:
            # - element 0 is the number of times the split was seen over all sampled trees
            #     (this provides the posterior probability of this split when divided by the
            #     number of trees sampled)
            # - element 1 holds the sum of edge lengths for this split (this provides the posterior
            #     mean edge length corresponding to this split when divided by the value in element 0)
            # - elements 2... are, for internal nodes, the indices of trees in which the split was
            #     found (this part is not absent for terminal nodes, which must appear in every tree)
            # From the tree index list we can compute the following quantities:
            # 1. the number of sojourns for each split (a sojourn is a series of consecutive
            #    samples in which a split was found preceded and followed by at least one sample
            #    not containing the split)
            # 2. sliding window and cumulative plots, such as those produced by AWTY
            nd = t.getFirstPreorder()
            assert nd.isRoot(), 'the first preorder node should be the root'
            split_vec = []
            while True:
                nd = nd.getNextPreorder()
                if not nd:
                    break
                else:
                    # Determine whether this split represents an internal or tip node
                    tip_node = nd.isTip() or nd.getParent().isRoot()
                    internal_node = not tip_node
                    
                    # Grab the edge length
                    edge_len = nd.getEdgeLen()
                    
                    # Grab the split and invert it if necessary to attain a standard polarity
                    s = nd.getSplit()
                    if s.isBitSet(0):
                        s.invertSplit()
                        
                    # Create a string representation of the split
                    ss = s.createPatternRepresentation()
                    
                    # Add string represention of the split to the tree_key list, which
                    # will be used to uniquely identify the tree topology
                    if internal_node:                        
                        tree_key.append(ss)

                    # Update the dictionary entry corresponding to this split
                    if ss in split_map.keys():
                        # split has been seen before
                        split_map[ss][0] += 1
                        split_map[ss][1] += edge_len
                        if internal_node:
                            split_map[ss].append(num_trees_considered)
                    else:
                        # split not yet seen
                        if tip_node:
                            split_map[ss] = [1, edge_len]
                        else:
                            split_map[ss] = [1, edge_len, num_trees_considered]
                            
            # Update tree_map, which is a map with keys equal to lists of internal node splits
            # and values equal to 2-element lists containing the frequency and newick tree
            # description
            tree_key.sort()
            k = tuple(tree_key)
            if k in tree_map.keys():
                tree_map[k][0] += 1
            else:
                tree_map[k] = [1,newick]

        self.phycas.output('\nSummary of sampled trees:')
        self.phycas.output('-------------------------')
        self.phycas.output('Tree file = %s' % self.phycas.sumt_input_tree_file)
        self.phycas.output('Total number of trees in file = %d' % num_trees)
        self.phycas.output('Number of trees considered = %d' % num_trees_considered)
        self.phycas.output('Number of distinct tree topologies found = %d' % len(tree_map.keys()))
        self.phycas.output('Number of distinct splits found = %d' % (len(split_map.keys()) - t.getNTips()))

        # Sort the splits from highest posterior probabilty to lowest        
        split_vect = split_map.items()
        split_vect.sort(cmp = lambda x,y: cmp(y[1][0], x[1][0]))

        # Output summary of splits
        self.phycas.output('\nSplit, frequency, mean edge len, first sojourn start, last sojourn end, and number of sojourns:')
        first_below_50 = None
        for i,(k,v) in enumerate(split_vect):
            trivial_split = len(v) == 2 and True or False
            
            # Split frequency is simply the first element of the list
            split_freq = v[0]

            # Split percent is split_freq divided by the number of trees considered
            split_pct = 100.0*float(split_freq)/float(num_trees_considered)

            # split_weight is the sum of edge lengths v[1] divided by the split_freq
            split_weight = float(v[1])/float(split_freq)

            # Identify the index of the first split having a posterior probability
            # lower than 0.5. The list up to this point will be used to generate the
            # majority-rule consensus tree
            if not first_below_50 and split_pct < 50.0:
                first_below_50 = i

            # Determine first sojourn (the third element of the list)
            first_sojourn_start = trivial_split and 1 or v[2]

            # Determine last sojourn (the final element of the list)
            last_sojourn_end = trivial_split and num_trees_considered or v[-1]

            # Determine the number of sojourns (this must be counted)
            num_sojourns = 1
            if not trivial_split:
                in_sojourn = True
                prev = v[2]
                for curr in v[3:]:
                    if curr - prev > 1:
                        num_sojourns += 1
                    prev = curr
                
            self.phycas.output('%s\t%d\t%.1f\t%.5f\t%d\t%d\t%d' % (k, split_freq, split_pct, split_weight, first_sojourn_start, last_sojourn_end, num_sojourns))

        # Build 50% majority rule tree if requested
        self.phycas.output('\nSaving majority-rule consensus tree to file %s' % self.phycas.sumt_output_tree_file)
        majrule = Phylogeny.Tree()
        tm = Phylogeny.TreeManip(majrule)
        majrule_splits = []
        for k,v in split_vect[:first_below_50]:
            if len(v) > 2:
                majrule_splits.append(k)
        tm.buildTreeFromSplitVector(majrule_splits, ProbDist.ExponentialDist(10))
        self.assignEdgeLensAndSupportValues(majrule, split_map, num_trees_considered)
        summary_short_name_list = ['majrule']
        summary_full_name_list = ['Majority-rule Consensus']
        summary_tree_list = [majrule]

        # Output summary of tree topologies
        self.phycas.output('\nTree topology, frequency, posterior probability, cumulative probability:')
        cum_prob = 0.0
        tree_vect = tree_map.items()
        tree_vect.sort(cmp = lambda x,y: cmp(y[1][0], x[1][0]))
        for i,(k,v) in enumerate(tree_vect):
            post_prob = float(v[0])/float(num_trees_considered)
            cum_prob += post_prob
            self.phycas.output('%d\t%d\t%.3f\t%.3f' % (i+1, v[0], post_prob, cum_prob))
            if self.phycas.sumt_output_tree_file:
                t = Phylogeny.Tree()
                t.buildFromString(v[1], True)   # these trees come from a nexus file, so taxon numbers will start at 0, hence the True second argument
                self.assignEdgeLensAndSupportValues(t, split_map, num_trees_considered)
                t.stripNodeNames()
                summary_short_name_list.append('%d_%d_of_%d' % (i+1,v[0],num_trees_considered))
                summary_full_name_list.append('Frequency %d of %d' % (v[0],num_trees_considered))
                summary_tree_list.append(t)
        self.phycas.phycassert(self.phycas.sumt_output_tree_file, 'sumt_output_tree_file should not be None')
        self.phycas.output('\nSaving tree topologies in file %s' % self.phycas.sumt_output_tree_file)
        self.save_trees(self.phycas.sumt_output_tree_file, self.phycas.sumt_output_pdf_file, summary_short_name_list, summary_full_name_list, summary_tree_list)
