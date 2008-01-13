import os,sys,math,random
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

    def save_trees(self, file_name, tree_names, trees):
        alt_file_name = file_name
        while os.path.exists(alt_file_name):
            alt_file_name = '%s_%d.tre' % (file_name, random.randint(1,1000000))

        if alt_file_name != file_name:
            self.phycas.output('Supplied tree file name %s exists, using %s instead' % (file_name, alt_file_name))
        treef = open(alt_file_name, 'w')
        treef.write('#nexus\n\n')
        treef.write('begin trees;\n')
        for name,tree in zip(tree_names, trees):
            treef.write('  tree %s = [&U] %s;\n' % (name, tree.makeNewick()))
        treef.write('end;\n')
        treef.close()        
        
    def consensus(self):
        num_trees = 0
        num_trees_considered = 0
        split_map = {}
        tree_map = {}

        # Open pdf_treefile and read trees therein
        self.readTreesFromFile()

        # Build each tree and build up split_map and tree_map
        t = Phylogeny.Tree()
        for newick in self.stored_newicks:
            if num_trees < self.phycas.sumt_burnin:
                num_trees += 1
                continue
            tree_key = []
            num_trees += 1
            num_trees_considered += 1
            t.buildFromString(newick, True)
            t.rectifyNames(self.taxon_labels)
            #if self.phycas.sumt_outgroup_taxon:
            #    num = t.findTipByName(self.phycas.sumt_outgroup_taxon)
            #    self.phycas.phycassert(num, 'could not root tree using specified outgroup: no tip having name "%s" could be found' % self.phycas.sumt_outgroup_taxon)
            #    t.rerootAtTip(num)
            t.recalcAllSplits(t.getNTips())

            # Each split found in any tree is associated with a list by the dictionary split_map
            # The list holds a series of integers, each of which is the index of a tree in which
            # that split was found. From such lists, we can compute the following quantities:
            # 1. the number of sojourns for each split (a sojourn is a series of consecutive
            #    samples in which a split was found preceded and followed by at least one sample
            #    not containing the split)
            # 2. the posterior probability of the split (list length divided by number of trees
            #    sampled)
            # 3. sliding window and cumulative plots, such as those produced by AWTY
            nd = t.getFirstPreorder()
            assert nd.isRoot(), 'first preorder node should be the root'
            nd = nd.getNextPreorder()   # subroot split is a trivial split, which can be ignored
            split_vec = []
            while True:
                nd = nd.getNextPreorder()
                if not nd:
                    break
                elif nd.isInternal():
                    s = nd.getSplit()
                    if s.isBitSet(0):
                        s.invertSplit()
                    ss = s.createPatternRepresentation()
                    tree_key.append(ss)
                    if ss in split_map.keys():
                        split_map[ss].append(num_trees_considered)
                        split_map[ss][0] += 1
                    else:
                        # split not yet seen, so create list for this split
                        # Note: first element reserved for storing number of
                        # times split appears in trees sampled (i.e. first
                        # element will equal length of list minus 1).
                        split_map[ss] = [1,num_trees_considered]
            tree_key.sort()
            k = tuple(tree_key)
            if k in tree_map.keys():
                tree_map[k][0] += 1
            else:
                tree_map[k] = [1,newick]

        self.phycas.output('\nSummary of sampled trees:')
        self.phycas.output('-------------------------')
        self.phycas.output('Tree file = %s' % self.phycas.sumt_tfile_name)
        self.phycas.output('Total number of trees in file = %d' % num_trees)
        self.phycas.output('Number of trees considered = %d' % num_trees_considered)
        self.phycas.output('Number of distinct tree topologies found = %d' % len(tree_map.keys()))
        self.phycas.output('Number of distinct splits found = %d' % len(split_map.keys()))

        self.phycas.output('\nSplit, frequency, first sojourn start, last sojourn end, and number of sojourns:')
        split_vect = split_map.items()
        split_vect.sort(cmp = lambda x,y: cmp(y[1][0], x[1][0]))
        first_below_50 = None
        for i,(k,v) in enumerate(split_vect):
            # Split frequency is simply the first element of the list
            split_freq = v[0]

            # Split percent is split_freq divided by the number of trees considered
            split_pct = 100.0*float(split_freq)/float(num_trees_considered)

            if not first_below_50 and split_pct < 50.0:
                first_below_50 = i

            # Determine first sojourn is simply the second element of the list
            first_sojourn_start = v[1]

            # Determine first sojourn is simply the final element of the list
            last_sojourn_end = v[-1]

            # Number of sojourns must be counted
            num_sojourns = 1
            in_sojourn = True
            prev = v[1]
            for curr in v[2:]:
                if curr - prev > 1:
                    num_sojourns += 1
                prev = curr
                
            self.phycas.output('%s\t%d\t%.1fd\t%d\t%d\t%d' % (k, split_freq, split_pct, first_sojourn_start, last_sojourn_end, num_sojourns))

        # Build 50% majority rule tree
        if self.phycas.sumt_majrule_tree_file:
            self.phycas.output('\nSaving majority-rule consensus tree to file %s' % self.phycas.sumt_majrule_tree_file)
            majrule = Phylogeny.Tree()
            tm = Phylogeny.TreeManip(majrule)
            edge_len_dist = ProbDist.ExponentialDist(10)
            majrule_splits = [k for k,v in split_vect[:first_below_50]]
            tm.buildTreeFromSplitVector(majrule_splits, edge_len_dist)
            self.save_trees(self.phycas.sumt_majrule_tree_file, ('majrule',),(majrule,))

        self.phycas.output('\nTree topology, frequency, posterior probability, cumulative probability:')
        cum_prob = 0.0
        tree_vect = tree_map.items()
        tree_vect.sort(cmp = lambda x,y: cmp(y[1][0], x[1][0]))
        for i,(k,v) in enumerate(tree_vect):
            post_prob = float(v[0])/float(num_trees_considered)
            cum_prob += post_prob
            self.phycas.output('%d\t%d\t%.3f\t%.3f' % (i + 1, v[0], post_prob, cum_prob))
        
        self.phycas.output('\nKey to tree topologies:')
        for i,(k,v) in enumerate(tree_vect):
            self.phycas.output('%d\t%s' % (i + 1, v[1]))

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
        #pdf.saveDocument(self.phycas.sumt_majrule_pdf_file)
