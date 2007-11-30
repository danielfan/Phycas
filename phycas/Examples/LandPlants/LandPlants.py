from phycas import *
import math
import os.path,shutil

data_file_name = 'Yang_and_Rannala_Karol.nex'
rnseed = 15397
ncycles = 20000
samplefreq = 20
lsweight = 100
nreps = 1
prior_means = [10**x for x in range(-8,4)]

####################################################################
#################### Analyze the simulated data ####################
####################################################################
def analyze(prior_mean):
    analyzer = Phycas()
    analyzer.random_seed = rnseed
    analyzer.default_model = 'jc'
    analyzer.num_rates = 1
    analyzer.edgelen_dist = ProbDist.ExponentialDist(1.0/prior_mean)
    analyzer.using_hyperprior = False 
    analyzer.allow_polytomies = False
    analyzer.data_file_name = data_file_name
    analyzer.starting_tree_source = 'random'
    analyzer.nchains = 1
    analyzer.ncycles = ncycles
    analyzer.sample_every = samplefreq
    analyzer.report_every = ncycles/20
    analyzer.gg_kvect = [1.0]
    analyzer.gg_outfile = None
    analyzer.ls_move_weight = lsweight
    #analyzer.mcmc()
    analyzer.gg_bin_patterns = True
    analyzer.gg_nreps = nreps
    analyzer.gg_pfile = data_file_name + '.p'
    analyzer.gg_tfile = data_file_name + '.t'
    shutil.copyfile(analyzer.gg_pfile, '%.9f_%s' % (prior_mean, analyzer.gg_pfile))
    shutil.copyfile(analyzer.gg_tfile, '%.9f_%s' % (prior_mean, analyzer.gg_tfile))
    analyzer.gg()
    return (analyzer.gg_Pm, analyzer.gg_Gm[0])

if __name__ == '__main__':
    # Open a file that will save the Gelfand-Ghosh values
    outf = file('output.txt', 'w')  # the 'w' means "open the file for writing"
    outf.write('prior mean\tPm\tGm\tDm\n')    # '\t' means tab, '\n' means newline

    # For each value of prior mean, perform an MCMC analysis
    for mean in prior_means:
        Pm, Gm = analyze(mean)
        Dm = Pm + Gm
        print '**'
        print '** Analysis with prior mean',mean,'yielded Pm =',Pm,', Gm =',Gm
        print '**'

        # Write out a line comprising four values separated by tabs
        # Each %f represents one floating point value from the 4-valued tuple that follows
        # Each '\t' is a tab character, and the final '\n' starts a new line in the file
        outf.write('%.9f\t%f\t%f\t%f\n' % (mean, Pm, Gm, Dm))

    # Close the file and we're done!
    outf.close()
