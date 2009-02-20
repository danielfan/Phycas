import os,sys
from phycas import *
data_file_name = 'Yang_and_Rannala_Karol.nex'
internal_prior_means = [10**x for x in [-5,-4,-3,-2,-1,0,1,2,3,4,5]]

if os.path.exists(os.path.abspath('output.txt')):
    print 'Aborting because output.txt already exists (remove or rename it and try again)'
    sys.exit(0)
outf = open('output.txt', 'w')
outf.write('prior mean\tPm\tGm\tDm\n')
outf.close()

for mean in internal_prior_means:
    print
    print '*****'
    print '***** Beginning analysis with internal prior mean',mean
    print '*****'
    analyzer = Phycas()
    analyzer.log_file_name = '%.9f_gg.log' % mean
    analyzer.data_file_name = data_file_name
    analyzer.default_model = 'hky'
    analyzer.num_rates = 5
    analyzer.using_hyperprior = False 
    analyzer.nchains = 1
    analyzer.gg_kvect = [1.0]
    analyzer.gg_burnin = 2556
    analyzer.gg_bin_patterns = True
    analyzer.gg_bincount_filename = '%.9f_bincounts.txt' % mean
    analyzer.gg_nreps = 1
    analyzer.gg_pfile = '%.9f_%s.p' % (mean, data_file_name)
    analyzer.gg_tfile = '%.9f_%s.t' % (mean, data_file_name)
    analyzer.gg()
    Pm = analyzer.gg_Pm
    Gm = analyzer.gg_Gm[0]
    Dm = Pm + Gm
    print '*****'
    print '***** Analysis with internal prior mean',mean,'yielded Pm =',Pm,', Gm =',Gm
    print '*****'

outf = open('output.txt', 'a')
outf.write('%.9f\t%f\t%f\t%f\n' % (mean, Pm, Gm, Dm))
outf.close()
	