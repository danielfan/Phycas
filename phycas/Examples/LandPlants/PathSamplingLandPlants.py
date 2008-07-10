from phycas import *
import shutil

data_file_name = 'Yang_and_Rannala_Karol.nex'
rnseed = 15397
ncycles = 600
samplefreq = 10
lsweight = 100
external_prior_mean = 0.1
internal_prior_means = [10**x for x in [3]]

def analyze(internal_prior_mean):
    analyzer = Phycas()
    analyzer.random_seed = rnseed
    analyzer.default_model = 'hky'
    analyzer.num_rates = 3
    analyzer.external_edgelen_dist = ProbDist.ExponentialDist(1.0/external_prior_mean)
    analyzer.internal_edgelen_dist = ProbDist.ExponentialDist(1.0/internal_prior_mean)
    analyzer.using_hyperprior = False 
    analyzer.allow_polytomies = False
    analyzer.data_file_name = data_file_name
    analyzer.outfile_prefix = 'ps_%.9f' % internal_prior_mean
    analyzer.starting_tree_source = 'random'
    analyzer.is_standard_heating = False
    analyzer.adapt_first = 10
    analyzer.verbose = True
    analyzer.tree_scaler_weight = 1
    analyzer.slice_weight = 1
    analyzer.pinvar_model = False
    analyzer.use_flex_model = False
    analyzer.ps_toward_posterior = False
    analyzer.ps_burnin = 0
    analyzer.ps_Q = ncycles
    analyzer.ps_nbetaincr = 100
    analyzer.sample_every = samplefreq
    analyzer.report_every = ncycles/20
    analyzer.ls_move_weight = lsweight
    analyzer.pathsampling()
    shutil.copyfile('%s.p' % outfile_prefix, '%.9f_%s.p' % (internal_prior_mean, data_file_name))
    shutil.copyfile('%s.t' % outfile_prefix, '%.9f_%s.t' % (internal_prior_mean, data_file_name))

if __name__ == '__main__':
    for mean in internal_prior_means:
        analyze(mean)
