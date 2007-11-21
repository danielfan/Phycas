from phycas import *

phycas = Phycas()
phycas.starting_tree_source = 'usertree'
phycas.pdf_newick           = "('P. fimbriata':0.1,('P. articulata':0.09,'P. parksii':0.04)v:0.1,('P. macrophylla':0.14,(('P. gracilis':0.08,('P. ciliata':0.02,'P. basiramia':0.03)x:0.05)y:0.01,'P. polygama':0.09)z:0.06)w:0.07)u"
phycas.pdf_tip_label_font   = 'Times-Italic'
phycas.pdf_tip_label_height = 24
phycas.pdf_page_width       = 8.5
phycas.pdf_page_height      = 11.0
phycas.pdf_line_width       = 2.0
phycas.pdf_left_margin      = 1.5
phycas.pdf_right_margin     = 1.5
phycas.pdf_top_margin       = 1.5
phycas.pdf_bottom_margin    = 1.5
phycas.pdf_filename         = 'test.pdf'
phycas.pdftree()

