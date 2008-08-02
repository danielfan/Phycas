from phycas.Utilities.PDFTree import PDFTree

# This example plots the trees in compilation.tre all to the same scale, saving
# the result in compilation.pdf (one tree per page)

# The PDFTree class is used for this (note that the word plotter is arbitrary)
plotter = PDFTree()

# Specify the filename of the tree file containing the trees to plot
plotter.pdf_treefile = 'compilation.tre'

# Specify the filename of the PDF file to output
plotter.pdf_filename = 'compilation.pdf'

# Specify the outgroup taxon to use when plotting trees
plotter.pdf_outgroup_taxon = '40 Cyanophora paradoxa'

# Unlink horizontal and vertical dimensions so trees that are not as wide
# as others can nevertheless be just as tall
plotter.keep_xy_proportional = False

# Unlink the size of the tip labels from the size of the tree
# (i.e. trees that are smaller have the same size tip labels 
# as trees that are larger)
plotter.keep_tip_labels_proportional = False

# Create the compilation.pdf file
plotter.pdftree()
