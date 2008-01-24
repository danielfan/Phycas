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

    def obsoletePDFScatterPlot(self, pdf, xinfo, yinfo, data, lines = True, w = 11.0, h = 8.5):
        xmin    = xinfo[0]
        xmax    = xinfo[1]
        xdiff   = xmax - xmin
        xdiv = xinfo[2]
        xfmt = '%%.%df' % xinfo[3]
        xlabels = [xfmt % (xmin + xdiff*float(i)/float(xdiv)) for i in range(xdiv + 1)]
        print 'xlabels =',xlabels
        
        ymin    = yinfo[0]
        ymax    = yinfo[1]
        ydiff   = ymax - ymin
        ydiv = yinfo[2]
        yfmt = '%%.%df' % yinfo[3]
        ylabels = [yfmt % (ymin + ydiff*float(i)/float(ydiv)) for i in range(ydiv + 1)]
        print 'ylabels =',ylabels
        
        pdf.newPage()
        
        inch          = 72.0
        spacer        = 5.0
        half_tick     = 2.0
        font_height   = float(self.phycas.pdf_plot_label_height)
        label_font    = self.phycas.pdf_plot_label_font
        line_width    = self.phycas.pdf_line_width

        xlabel_height = font_height*pdf.getXHeight(label_font)
        xlabel_width = 0.0
        for label in xlabels:
            width = font_height*pdf.calcStringWidth(label_font, label)
            if width > xlabel_width:
                xlabel_width = width

        ylabel_width = 0.0
        for label in ylabels:
            width = font_height*pdf.calcStringWidth(label_font, label)
            if width > ylabel_width:
                ylabel_width = width
        
        right_margin  = inch*self.phycas.pdf_right_margin
        left_margin   = inch*self.phycas.pdf_left_margin
        top_margin    = inch*self.phycas.pdf_top_margin
        bottom_margin = inch*self.phycas.pdf_bottom_margin

        plot_bottom   = bottom_margin + xlabel_height + 2.0*spacer
        plot_left     = left_margin + ylabel_width + spacer
        plot_width    = inch*w - plot_left - right_margin - xlabel_width/2.0
        plot_height   = inch*h - plot_bottom - top_margin - xlabel_height/2.0

        #pdf.addRectangle(left_margin, bottom_margin, inch*w - left_margin - right_margin, inch*h - top_margin - bottom_margin, 1, 'dotted')

        # draw the x-axis
        pdf.addLine(plot_left, plot_bottom, plot_left + plot_width, plot_bottom, line_width)

        # draw the x-axis tick marks and labels for 0, ..., ntrees
        for i,label in enumerate(xlabels):
            x = plot_width*float(i)/float(len(xlabels) - 1)
            half_width = font_height*pdf.calcStringWidth(label_font, label)/2.0
            pdf.addText(plot_left + x - half_width, plot_bottom - xlabel_height - 2.0*spacer, label_font, font_height, label)
            pdf.addLine(plot_left + x, plot_bottom - half_tick, plot_left + x, plot_bottom + half_tick, line_width)
        
        # draw the y-axis
        pdf.addLine(plot_left, plot_bottom, plot_left, plot_bottom + plot_height, line_width)

        # draw the y-axis tick marks and labels for 0.0, 0.1, ..., 1.0
        for i,label in enumerate(ylabels):
            y = plot_height*float(i)/float(len(ylabels) - 1)
            ylabel = y - xlabel_height/2.0
            pdf.addText(left_margin, plot_bottom + ylabel, label_font, font_height, label)
            pdf.addLine(plot_left - half_tick, plot_bottom + y, plot_left + half_tick, plot_bottom + y, line_width)
            
        # draw points or lines
        # data should be in the form of a list of lists of tuples
        # - each outer list defines a separate line in the plot
        #     e.g. [line1,line2,line3]
        # - the next level defines the points in a line
        #     e.g. line1 = [point1, point2, point3]
        # - the final level defines the x- and y-coordinates of a point
        #     e.g. point1 = (x,y)
        #     x and y should be in the same units used to label the axes. That is, if 
        #     xmin = 0.0 and xmax = 1.0, then the x values should be between 0.0 and 1.0
        for line_data in data:
            first_point = True
            for point in line_data:
                if first_point:
                    x0 = plot_width*(point[0] - xmin)/xdiff
                    y0 = plot_height*(point[1] - ymin)/ydiff
                    #print 'x0 =',point[0],', y0 =',point[1]
                    first_point = False
                else:
                    x = plot_width*(point[0] - xmin)/xdiff
                    y = plot_height*(point[1] - ymin)/ydiff
                    pdf.addLine(plot_left + x0, plot_bottom + y0, plot_left + x, plot_bottom + y, line_width)
                    #print 'x  =',point[0],', y  =',point[1]
                    x0 = x
                    y0 = y
            #raw_input('check')

    def awtyPlot(self, split_vect, ntrees):
        splits_to_plot = 1000
        ignore_uninteresting = True

        pdf = PDFGenerator(11.0, 8.5)
        pdf.overwrite = True
        
        splits_plotted = 0
        trivial_ignored = 0
        uninteresting_ignored = 0
        data = []
        for i,(k,v) in enumerate(split_vect):
            if len(v) == 2:
                trivial_ignored += 1
                continue
            if ignore_uninteresting and v[0] == ntrees:
                uninteresting_ignored += 1
                continue
            line_data = [(0.0,0.0)]
            for j,m in enumerate(v[2:]):
                x = float(m)
                y = float(j + 1)/float(m)
                line_data.append((x,y))
            x = ntrees
            y = float(v[0])/float(ntrees)
            line_data.append((x,y))
            data.append(line_data)
            
            splits_plotted += 1
            if splits_plotted == splits_to_plot:
                break

        pdf.scatterPlot(data, xinfo = (0,ntrees,10,0), yinfo = (0.0,1.0,10,1))
        pdf.saveDocument('awty.pdf')
        if trivial_ignored + uninteresting_ignored > 0:
            self.phycas.output('%d trivial and %d uninteresting splits were ignored.' % (trivial_ignored, uninteresting_ignored))
        
    def obsoleteAWTYPlot(self, split_vect, ntrees):
        #pdf = PDFGenerator(self.phycas.pdf_page_width, self.phycas.pdf_page_height)
        w             = 11.0
        h             = 8.5
        pdf = PDFGenerator(w, h)
        pdf.overwrite = True
        pdf.newPage()

        splits_to_plot = 1000
        ignore_uninteresting = True

        inch          = 72.0
        spacer        = 5.0
        half_tick     = 2.0
        font_height   = float(self.phycas.pdf_plot_label_height)
        label_font    = self.phycas.pdf_plot_label_font
        line_width    = self.phycas.pdf_line_width

        xlabel_height = font_height*pdf.getXHeight(label_font)
        ylabel_width  = font_height*pdf.calcStringWidth(label_font, '0.0')
        
        right_margin  = inch*self.phycas.pdf_right_margin
        left_margin   = inch*self.phycas.pdf_left_margin
        top_margin    = inch*self.phycas.pdf_top_margin
        bottom_margin = inch*self.phycas.pdf_bottom_margin

        plot_bottom   = bottom_margin + xlabel_height + 2.0*spacer
        plot_left     = left_margin + ylabel_width + spacer
        plot_width    = inch*w - plot_left - right_margin
        plot_height   = inch*h - plot_bottom - top_margin

        #pdf.addRectangle(left_margin, bottom_margin, inch*w - left_margin - right_margin, inch*h - top_margin - bottom_margin, 1, 'dotted')

        # draw the x-axis
        pdf.addLine(plot_left, plot_bottom, plot_left + plot_width, plot_bottom, line_width)

        # draw the x-axis tick marks and labels for 0, ..., ntrees
        for i in range(11):
            v = math.ceil(float(i*ntrees)/10.0)
            x = plot_width*float(i)/10.0
            half_width = font_height*pdf.calcStringWidth(label_font, '%d' % v)/2.0
            pdf.addText(plot_left + x - half_width, plot_bottom - xlabel_height - 2.0*spacer, label_font, font_height, '%d' % v)
            pdf.addLine(plot_left + x, plot_bottom - half_tick, plot_left + x, plot_bottom + half_tick, line_width)
        
        # draw the y-axis
        pdf.addLine(plot_left, plot_bottom, plot_left, plot_bottom + plot_height, line_width)

        # draw the y-axis tick marks and labels for 0.0, 0.1, ..., 1.0
        for i in range(11):
            v = 0.1*float(i)
            y = v*plot_height
            ylabel = y - xlabel_height/2.0
            pdf.addText(left_margin, plot_bottom + ylabel, label_font, font_height, '%.1f' % v)
            pdf.addLine(plot_left - half_tick, plot_bottom + y, plot_left + half_tick, plot_bottom + y, line_width)
            
        # draw lines for each split
        splits_plotted = 0
        trivial_ignored = 0
        uninteresting_ignored = 0
        for i,(k,v) in enumerate(split_vect):
            if len(v) == 2:
                trivial_ignored += 1
                continue
            if ignore_uninteresting and v[0] == ntrees:
                uninteresting_ignored += 1
                continue
            x0 = 0.0
            y0 = 0.0
            for j,m in enumerate(v[2:]):
                x = plot_width*float(m)/float(ntrees)
                y = plot_height*float(j + 1)/float(m)
                pdf.addLine(plot_left + x0, plot_bottom + y0, plot_left + x, plot_bottom + y, line_width)
                x0 = x
                y0 = y
            split_postprob = float(v[0])/float(ntrees)
            #self.phycas.output('%6d %6d %.3f %6d %s' % (i,v[0],split_postprob,len(v), k))
            x = plot_width
            y = plot_height*split_postprob
            pdf.addLine(plot_left + x0, plot_bottom + y0, plot_left + x, plot_bottom + y, line_width)
            #pdf.addText(plot_left + plot_width, plot_bottom + y, label_font, font_height, '%.3f' % split_postprob)
            
            splits_plotted += 1
            if splits_plotted == splits_to_plot:
                break

        pdf.saveDocument('awty.pdf')
        if trivial_ignored + uninteresting_ignored > 0:
            self.phycas.output('%d trivial and %d uninteresting splits were ignored.' % (trivial_ignored, uninteresting_ignored))
        
    def consensus(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This is the main member function of the class. It computes and outputs
        a summary of splits found in the sampled trees, a summary of distinct
        tree topologies found, and a majority rule consensus tree, which is
        saved to the file sumt_output_tree_file if this variable is not None.
        
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

        # values used for display purposes
        split_field_width = 0
        num_stored_trees = len(self.stored_newicks)
        self.phycas.phycassert(num_stored_trees > 0, 'Specified tree file (%s) contained no stored trees' % self.phycas.sumt_input_tree_file)
        sojourn_field_width = 2 + math.floor(math.log10(float(num_stored_trees)))
        
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
            ntips = t.getNTips()
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
                tree_map[k] = [1, newick, treelen, num_trees_considered]

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
        self.phycas.output('\nTree topology (topology), frequency (freq.), mean tree length (TL),')
        self.phycas.output('  first sojourn start (s0), last sojourn end (sk), number of sojourns (k),')
        self.phycas.output('  posterior probability (prob.), cumulative probability (cum.):')
        freq_str = sojourn_label_fmt_str % 'freq.'
        s0_str = sojourn_label_fmt_str % 's0'
        sk_str = sojourn_label_fmt_str % 'sk'
        k_str = sojourn_label_fmt_str % 'k'
        self.phycas.output('%8s %s %10s %s %s %s %10s %10s' % ('topology', freq_str, 'TL', s0_str, sk_str, k_str, 'prob.', 'cum.'))
        cum_prob = 0.0
        tree_vect = tree_map.items()
        tree_vect.sort(cmp = lambda x,y: cmp(y[1][0], x[1][0]))
        for i,(k,v) in enumerate(tree_vect):
            # Determine the posterior probability and cumulative posterior probability
            post_prob = float(v[0])/float(num_trees_considered)
            cum_prob += post_prob

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
            if self.phycas.sumt_output_tree_file:
                t = Phylogeny.Tree()
                t.buildFromString(v[1], True)   # these trees come from a nexus file, so taxon numbers will start at 0, hence the True second argument
                self.assignEdgeLensAndSupportValues(t, split_map, num_trees_considered)
                t.stripNodeNames()
                summary_short_name_list.append('%d_%d_of_%d' % (i+1,v[0],num_trees_considered))
                summary_full_name_list.append('Frequency %d of %d' % (v[0],num_trees_considered))
                summary_tree_list.append(t)
        self.phycas.phycassert(self.phycas.sumt_output_tree_file, 'sumt_output_tree_file should not be None')

        self.phycas.output('\nSaving tree topologies in file %s...' % self.phycas.sumt_output_tree_file)
        self.save_trees(self.phycas.sumt_output_tree_file, self.phycas.sumt_output_pdf_file, summary_short_name_list, summary_full_name_list, summary_tree_list)

        self.phycas.output('\nSaving AWTY plot in file awty.pdf...')
        self.awtyPlot(split_vect, num_trees_considered)

        self.phycas.output('\nSumT finished.')
