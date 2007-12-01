from phycas import *

data_file_name = 'Yang_and_Rannala_Karol.nex'
internal_prior_means = [10**-8]
#internal_prior_means = [10**x for x in range(-8,-6)]

outf = open('output.txt', 'w')
outf.write('prior mean\tPm\tGm\tDm\n')
outf.close()

for mean in internal_prior_means:
    analyzer = Phycas()
    analyzer.data_file_name = data_file_name
    analyzer.default_model = 'hky'
    analyzer.num_rates = 5
    analyzer.using_hyperprior = False 
    analyzer.nchains = 1
    analyzer.gg_kvect = [1.0]
    analyzer.gg_outfile = None
    analyzer.gg_burnin = 10
    analyzer.gg_bin_patterns = True
    analyzer.gg_bincount_filename = '%.9f_bincounts.txt' % mean
    analyzer.gg_nreps = 1
    analyzer.gg_pfile = '%.9f_%s.p' % (mean, data_file_name)
    analyzer.gg_tfile = '%.9f_%s.t' % (mean, data_file_name)
    analyzer.gg()
    Pm = analyzer.gg_Pm
    Gm = analyzer.gg_Gm[0]
    Dm = Pm + Gm
    print '***** Analysis with internal prior mean',mean,'yielded Pm =',Pm,', Gm =',Gm

    outf = open('output.txt', 'a')
    outf.write('%.9f\t%f\t%f\t%f\n' % (mean, Pm, Gm, Dm))
    outf.close()
