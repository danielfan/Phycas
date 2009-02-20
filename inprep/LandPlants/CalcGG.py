from phycas import *
from math import exp

phycas = Phycas()
phycas.data_file_name = 'Yang_and_Rannala_Karol_nomissambig.nex'
phycas.gg_pfile = 'Yang_and_Rannala_Karol_nomissambig.nex.p'
phycas.gg_tfile = 'Yang_and_Rannala_Karol_nomissambig.nex.t'
phycas.gg_outfile = 'ggout.txt'
phycas.gg_nreps = 2
phycas.gg_kvect = [1.0]
phycas.gg_burnin = 2991
phycas.gg_save_postpreds = False
phycas.gg()
phycas.shutdown()