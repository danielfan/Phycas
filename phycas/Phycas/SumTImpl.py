import os,sys,math,random
from phycas.PDFGen import *
#from phycas.TreeViewer import *
from phycas import *

class TreeSummarizer(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Saves consensus tree and QQ plots in pdf files.
    
    """
    def __init__(self, phycas, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes TreeSummarizer object by assigning supplied phycas object
        to a data member variable.
        
        """
        self.phycas = phycas
        self.opts = opts
        self.rooted_trees = opts.rooted
        self.outgroup = opts.outgroup_taxon
        if self.rooted_trees and opts.outgroup_taxon:
            self.outgroup = None
            self.phycas.warning('Specifying True for sumt_rooted is incompatible with specifying\nsumt_outgroup_taxon; I will pretend that you set sumt_outgroup_taxon to None')

    def readTreesFromFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reads all trees from the file named in self.opts.input_tree_file.
        Stores tree descriptions in self.stored_tree_defs and taxon labels in
        self.taxon_labels.
        
        """
        self.phycas.reader.readFile(self.opts.input_tree_file)
        self.taxon_labels = self.phycas.reader.getTaxLabels()
        self.stored_tree_defs = self.phycas.reader.getTrees()
        self.phycas.phycassert(len(self.stored_tree_defs) > 0, 'expecting a trees block defining at least one tree in the file %s' % self.opts.input_tree_file)

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
        #tree.recalcAllSplits(tree.getNTips())
        tree.recalcAllSplits(tree.getNObservables())
        nd = tree.getFirstPreorder()
        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            s = nd.getSplit()
            if (not self.rooted_trees) and s.isBitSet(0):
                s.invertSplit()
            ss = s.createPatternRepresentation()
            try:
                v = split_map[ss]
            except KeyError:
                self.phycas.abort('could not find edge length information for the split %s' % ss)
            support = float(v[0])/float(total_samples)
            nd.setSupport(support)
            if self.opts.equal_brlens:
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
        if self.opts.output_replace:
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
        if not self.opts.output_replace:
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
            if self.outgroup:
                try:
                    num = list(self.taxon_labels).index(self.outgroup)
                except ValueError:
                    self.phycas.output('Warning: could not root tree using specified outgroup: no tip having name "%s" could be found' % self.outgroup)
                    self.phycas.output('Here is the list of all taxon labels:')
                    for label in self.taxon_labels:
                        self.phycas.output('  %s' % label)
                else:
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

    def awtyPlot(self, pdf, split_vect, ntrees, ndivisions = 10):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The supplied split_vect is a vector of key,value pairs, where the key
        is a split pattern string and v=value[2:] is a list of samples in
        which the split was found. The minimum possible value in v is 1, and
        the maximum possible value is ntrees. Consider the following
        oversimplified example involving a split with posterior probability
        0.6 (v has length 6 and ntrees = 10): v = [3,4,5,7,9,10]. Suppose we
        wish to plot the cumulative posterior probability of this split at
        five, evenly-spaced points. Dividing ntrees by 5 gives 2, and these
        are thus the values to be plotted:
        x     y
        --------   x = range(0, ntrees + 1, ntrees//ndivisions) 
        0    0.0     = [0, 2, 4, 6, 8, 10]
        2    0.0   
        4    0.5   To calculate y, let k track the number of values in v
        6    0.5   that are less than x
        8    0.5     
        10   0.6    
        --------
        The vector x can be calculated using the range function; however,
        computing y requires a loop:
        i  x[i]  k  v[k]  y = k/x[i] except first
        ----------------------------
        0    0   0   3           0.0 <- this one always 0.0
        1    2   0   3    0/2  = 0.0
        2    4   2   5    2/4  = 0.5
        3    6   3   7    3/6  = 0.5
        4    8   4   9    4/8  = 0.5
        5   10   6   -    6/10 = 0.6
        ----------------------------
        
        """
        # candidates for phycas variable status
        splits_to_plot = 1000
        ignore_uninteresting = True

        splits_plotted = 0
        trivial_ignored = 0
        uninteresting_ignored = 0
        data = []
        xvect = range(0, ntrees + 1, ntrees//ndivisions)
        for k,value in split_vect:
            if len(value) == 2:
                trivial_ignored += 1
                continue
            if ignore_uninteresting and value[0] == ntrees:
                uninteresting_ignored += 1
                continue
            line_data = [(0.0,0.0)]
            v = value[2:]
            k = 0
            for i,x in enumerate(xvect[1:]):
                while k < len(v) and v[k] <= x:
                    k += 1
                y = float(k)/float(x)
                line_data.append((x,y))
            data.append(line_data)
            
            splits_plotted += 1
            if splits_plotted == splits_to_plot:
                break

        if splits_plotted == 0:
            self.phycas.output('No AWTY plot created because no splits with posteriors less than 1.0 were found')
            return False
        
        pdf.scatterPlot(data, title = 'Split Probabilities Through Time', xinfo = (0,ntrees,10,0), yinfo = (0.0,1.0,10,1))

        try:
            import wx
        except:
            pass
        else:
            wxapp = wxPhycas.CreateAWTYApp()
            wxapp.scatterPlot(data, title = 'Split Probabilities Through Time', xinfo = (0,ntrees,10,0), yinfo = (0.0,1.0,10,1))
            wxapp.MainLoop()
        
        if trivial_ignored + uninteresting_ignored > 0:
            self.phycas.output('%d trivial and %d uninteresting splits were ignored.' % (trivial_ignored, uninteresting_ignored))
        return True
        
    def sojournPlot(self, pdf, split_vect, ntrees):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a plot comprising numerous horizontal lines, each representing
        the trajectory of a non-trivial, "interesting" split throughout the
        MCMC simulation. Trivial splits are those separating a terminal taxon
        from the other taxa, and are aways present. "Interesting" splits are
        those that were absent during some of the non-burn-in period. This
        plot provides a visual picture of mixing in the Markov chain with
        respect to splits.
        
        """
        # candidates for phycas variable status
        splits_to_plot = 1000
        ignore_uninteresting = True

        trivial_ignored = 0
        uninteresting_ignored = 0

        # First, save the sojourn vectors of all the splits to be plotted
        v = []
        splits_plotted = 0
        for k,value in split_vect:
            if len(value) == 2:
                trivial_ignored += 1
                continue
            if ignore_uninteresting and value[0] == ntrees:
                uninteresting_ignored += 1
                continue
            v.append(value[2:])
            splits_plotted += 1
            if splits_plotted == splits_to_plot:
                break

        if splits_plotted == 0:
            self.phycas.output('No sojourn plot created because all non-trivial splits had posterior probability 1.0')
            return False

        # Now create the plot data
        data = []
        for i,s in enumerate(v):
            y = float(i + 1)
            prev_x = s[0]
            line_data = [(prev_x,y)]
            for x in s[1:]:
                if x - prev_x > 1:
                    # a gap has appeared, start a new line and close old one
                    data.append(line_data)
                    line_data = [(x,y)]
                else:
                    line_data.append((x,y))
                prev_x = x
            data.append(line_data)

        # Determine y-axis info (ymax and ydivs) such that ydivs is as close as possible to
        # 10 but divides evenly into ymax - ymin so that there are no fractional y-axis labels
        splits_per_tick = splits_plotted//10
        if splits_per_tick < 2:
            splits_per_tick = 2
        min_ticks = math.ceil(float(splits_plotted)/float(splits_per_tick))
        ymax = float(splits_per_tick)*min_ticks
        ydivs = int(ymax)//splits_per_tick
        if splits_plotted <= 10:
            ydivs = int(ymax)

        pdf.scatterPlot(data, points = False, line_width = 3, line_cap_style = 'square', title = 'Split Sojourns', xinfo = (0,ntrees,10,0), yinfo = (0.0,float(ymax),ydivs,0))
        if trivial_ignored + uninteresting_ignored > 0:
            self.phycas.output('%d trivial and %d uninteresting splits were ignored.' % (trivial_ignored, uninteresting_ignored))

        return True            
        
    def consensus(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This is the main member function of the class. It computes and outputs
        a summary of splits found in the sampled trees, a summary of distinct
        tree topologies found, and a majority rule consensus tree, which is
        saved to the file sumt_output_tree_file if this variable is not None.
        
        """
        #raw_input('debug stop')
        # Check to make sure user specified an input tree file
        self.phycas.phycassert(self.opts.input_tree_file, 'sumt_input_tree_file must be specified before sumt method is called')
        self.phycas.phycassert(self.opts.trees_prefix, 'sumt_trees_prefix must be specified before sumt method is called')
        self.phycas.phycassert(self.opts.trees_prefix != self.opts.splits_prefix, 'sumt_trees_prefix must differ from sumt_splits_prefix')
            
        treefname = '%s.tre' % self.opts.trees_prefix
        pdffname = '%s.pdf' % self.opts.trees_prefix

        num_trees = 0
        num_trees_considered = 0
        split_map = {}
        tree_map = {}

        # Open sumt_tfile_name and read trees therein
        self.phycas.output('\nReading trees from file %s...' % self.opts.input_tree_file)
        self.readTreesFromFile()

        # Build each tree and add the splits and tree topolgies found there to the
        # dictionary of splits (split_map) and the dictionary of tree topologies
        # (tree_map), respectively
        self.phycas.output('Compiling lists of tree topologies and splits...')
        t = Phylogeny.Tree()
        if self.rooted_trees:
            t.setRooted()

        # values used for display purposes
        split_field_width = 0
        num_stored_trees = len(self.stored_tree_defs)
        self.phycas.phycassert(num_stored_trees > 0, 'Specified tree file (%s) contained no stored trees' % self.opts.input_tree_file)
        sojourn_field_width = 2 + math.floor(math.log10(float(num_stored_trees)))
        
        for tree_def in self.stored_tree_defs:
            if num_trees < self.opts.burnin:
                num_trees += 1
                continue
            num_trees += 1
            num_trees_considered += 1

            # Create a tree key, which is a list of internal node splits that can be
            # used to uniquely identify a tree topology
            tree_key = []

            # Build the tree
            tree_def.buildTree(t)
            t.rectifyNames(self.taxon_labels)
            #ntips = t.getNTips()
            ntips = t.getNObservables()
            if ntips > split_field_width:
                # this is necessary only if number of taxa varies from tree to tree
                split_field_width = ntips
            t.recalcAllSplits(ntips)
            treelen = t.edgeLenSum()

            # Each split found in any tree is associated with a list by the dictionary split_map
            #   The list is organized as follows:
            #   - element 0 is the number of times the split was seen over all sampled trees
            #       (this provides the posterior probability of this split when divided by the
            #       number of trees sampled)
            #   - element 1 holds the sum of edge lengths for this split (this provides the posterior
            #       mean edge length corresponding to this split when divided by the value in element 0)
            #   - elements 2... are, for internal nodes, the indices of trees in which the split was
            #       found (this part is omitted for terminal nodes, which must appear in every tree)
            #   From the tree index list we can compute the following quantities for internal splits:
            #   1. the number of sojourns for each internal split (a sojourn is a series of consecutive
            #      samples in which a split was found that is preceded and followed by at least one
            #      sample not containing the split)
            #   2. sliding window and cumulative plots, such as those produced by AWTY
            # Each distinct tree topology is associated with a similar list:
            #   - element 0 is again the number of times the tree topology was seen over all samples
            #   - element 1 holds the newick tree description
            #   - element 2 holds the sum of tree lengths over all sampled trees of this topology
            #   - elements 3... are the indices of trees in which the topology was found. This list
            #       allows the number and extent of each sojourn to be computed
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
                    if (not self.rooted_trees) and s.isBitSet(0):
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
                        entry = split_map[ss]
                        entry[0] += 1
                        entry[1] += edge_len
                        if internal_node:
                            entry.append(num_trees_considered)
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
                # tree topology has been seen before
                entry = tree_map[k]
                entry[0] += 1
                entry[2] += treelen
                entry.append(num_trees_considered)
            else:
                # tree topology has not yet been seen
                tree_map[k] = [1, tree_def, treelen, num_trees_considered]

        self.phycas.output('\nSummary of sampled trees:')
        self.phycas.output('-------------------------')
        self.phycas.output('Tree file = %s' % self.opts.input_tree_file)
        self.phycas.output('Total number of trees in file = %d' % num_trees)
        self.phycas.output('Number of trees considered = %d' % num_trees_considered)
        self.phycas.output('Number of distinct tree topologies found = %d' % len(tree_map.keys()))
        #self.phycas.output('Number of distinct splits found = %d' % (len(split_map.keys()) - t.getNTips()))
        self.phycas.output('Number of distinct splits found = %d' % (len(split_map.keys()) - t.getNObservables()))

        # Sort the splits from highest posterior probabilty to lowest        
        split_vect = split_map.items()
        split_vect.sort(cmp = lambda x,y: cmp(y[1][0], x[1][0]))

        # Output summary of splits
        if sojourn_field_width < 5:
            # must be large enough to accommodate the column header 'freq.'
            sojourn_field_width = 5
        split_fmt_str = '%%%ds' % split_field_width
        sojourn_label_fmt_str = '%%%ds' % sojourn_field_width
        sojourn_fmt_str = '%%%dd' % sojourn_field_width
        self.phycas.output('\nSplit (split), split representation (pattern), frequency (freq.), posterior')
        self.phycas.output('  probability (prob.), mean edge length (weight), first sojourn start (s0),')
        self.phycas.output('  last sojourn end (sk), and number of sojourns (k):')
        split_str = split_fmt_str % 'pattern'
        freq_str = sojourn_label_fmt_str % 'freq.'
        s0_str = sojourn_label_fmt_str % 's0'
        sk_str = sojourn_label_fmt_str % 'sk'
        k_str = sojourn_label_fmt_str % 'k'
        self.phycas.output('%6s %s %s %10s %10s %s %s %s' % ('split', split_str, freq_str, 'prob.', 'weight', s0_str, sk_str, k_str))
        first_below_50 = None
        for i,(k,v) in enumerate(split_vect):
            trivial_split = len(v) == 2 and True or False
            
            # Split frequency is simply the first element of the list
            split_freq = v[0]

            # Split posterior is split_freq divided by the number of trees considered
            split_posterior = float(split_freq)/float(num_trees_considered)

            # split_weight is the sum of edge lengths v[1] divided by the split_freq
            split_weight = float(v[1])/float(split_freq)

            # Identify the index of the first split having a posterior probability
            # lower than 0.5. The list up to this point will be used to generate the
            # majority-rule consensus tree
            if not first_below_50 and split_posterior < 0.5:
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

            split_str = split_fmt_str % k            
            freq_str = sojourn_fmt_str % split_freq
            s0_str = sojourn_fmt_str % first_sojourn_start
            sk_str = sojourn_fmt_str % last_sojourn_end
            k_str = sojourn_fmt_str % num_sojourns
            self.phycas.output('%6d %s %s %10.5f %10.5f %s %s %s' % (i + 1, split_str, freq_str, split_posterior, split_weight, s0_str, sk_str, k_str))

        # Build 50% majority rule tree if requested
        self.phycas.output('\nSaving majority-rule consensus tree...')
        majrule = Phylogeny.Tree()
        if self.rooted_trees:
            majrule.setRooted()
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
        self.phycas.output('\nTree topology (topology), frequency (freq.), mean tree length (TL),')
        self.phycas.output('  first sojourn start (s0), last sojourn end (sk), number of sojourns (k),')
        self.phycas.output('  posterior probability (prob.), cumulative probability (cum.):')
        freq_str = sojourn_label_fmt_str % 'freq.'
        s0_str = sojourn_label_fmt_str % 's0'
        sk_str = sojourn_label_fmt_str % 'sk'
        k_str = sojourn_label_fmt_str % 'k'
        self.phycas.output('%8s %s %10s %s %s %s %10s %10s' % ('topology', freq_str, 'TL', s0_str, sk_str, k_str, 'prob.', 'cum.'))
        cum_prob = 0.0
        done = False
        tree_vect = tree_map.items()
        tree_vect.sort(cmp = lambda x,y: cmp(y[1][0], x[1][0]))
        if self.opts.tree_credible_prob > 0.0:
            for i,(k,v) in enumerate(tree_vect):
                if done:
                    break
                
                # Determine the posterior probability and cumulative posterior probability
                post_prob = float(v[0])/float(num_trees_considered)
                cum_prob += post_prob
                if cum_prob > self.opts.tree_credible_prob:
                    # stop after the current tree
                    done = True

                # Determine the average tree length of this tree topology
                avgTL = v[2]/float(v[0])

                # Determine the sampled tree that began the first sojourn (the third element of the list)
                first_sojourn_start = v[3]
                
                # Determine the sampled tree that ended the last sojourn (the final element of the list)
                last_sojourn_end = v[-1]

                # Determine the number of sojourns (this must be counted)
                num_sojourns = 1
                in_sojourn = True
                prev = v[3]
                for curr in v[4:]:
                    if curr - prev > 1:
                        num_sojourns += 1
                    prev = curr

                # Output summary line for this tree topology                
                freq_str = sojourn_fmt_str % v[0]
                s0_str = sojourn_fmt_str % first_sojourn_start
                sk_str = sojourn_fmt_str % last_sojourn_end
                k_str = sojourn_fmt_str % num_sojourns
                self.phycas.output('%8d %s %10.5f %s %s %s %10.5f %10.5f' % (i+1, freq_str, avgTL, s0_str, sk_str, k_str, post_prob, cum_prob))

                # Save the tree topology (decorated with posterior mean edge lengths) to the tree file
                t = Phylogeny.Tree()
                if self.rooted_trees:
                    t.setRooted()
                v[1].buildTree(t)
                self.assignEdgeLensAndSupportValues(t, split_map, num_trees_considered)
                t.stripNodeNames()
                summary_short_name_list.append('%d_%d_of_%d' % (i+1,v[0],num_trees_considered))
                summary_full_name_list.append('Frequency %d of %d' % (v[0],num_trees_considered))
                summary_tree_list.append(t)

        self.phycas.output('\nSaving distinct tree topologies...')
        self.save_trees(treefname, pdffname, summary_short_name_list, summary_full_name_list, summary_tree_list)

        if self.opts.splits_prefix:
            fn = '%s.pdf' % self.opts.splits_prefix
            pdf = PDFGenerator(11.0, 8.5)
            pdf.overwrite = True
            
            self.phycas.output('\nSaving AWTY plot in file %s...' % fn)
            awty_ok = self.awtyPlot(pdf, split_vect, num_trees_considered)

            self.phycas.output('\nSaving Sojourn plot in file %s...' % fn)
            sojourn_ok = self.sojournPlot(pdf, split_vect, num_trees_considered)

            if awty_ok or sojourn_ok:
                pdf.saveDocument(fn)

        self.phycas.output('\nSumT finished.')


