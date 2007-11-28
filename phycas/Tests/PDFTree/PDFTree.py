from phycas import *

phycas = Phycas()
phycas.pdf_newick                = "('P. fimbriata':0.1,('P. articulata':0.09,'P. parksii':0.04)v:0.1,('P. macrophylla':0.14,(('P. gracilis':0.08,('P. ciliata':0.02,'P. basiramia':0.03)x:0.05)y:0.01,'P. polygama':0.09)z:0.06)w:0.07)u"
phycas.pdf_outgroup_taxon        = 'P. macrophylla' # use None to display tree as is
phycas.pdf_tip_label_font        = 'Times-Italic'   # must be one of the 14 standard font names
phycas.pdf_tip_label_height      = 24               # in units of points
phycas.pdf_scalebar_label_font   = 'Times-Roman'    # must be one of the 14 standard font names
phycas.pdf_scalebar_label_height = 10               # in units of points
phycas.pdf_ladderize             = 'right'          # valid choices: None, 'left', or 'right'
phycas.pdf_scalebar_position     = 'bottom'         # value choices: None, 'top', or 'bottom'
phycas.pdf_page_width            = 8.5              # in units of inches
phycas.pdf_page_height           = 11.0             # in units of inches
phycas.pdf_line_width            = 2.0              # in units of points
phycas.pdf_left_margin           = 1.0              # in units of inches
phycas.pdf_right_margin          = 1.0              # in units of inches
phycas.pdf_top_margin            = 1.0              # in units of inches
phycas.pdf_bottom_margin         = 1.0              # in units of inches
phycas.pdf_filename              = 'test.pdf'
phycas.pdftree()

