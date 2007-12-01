from phycas import *
import shutil

data_file_name = 'Yang_and_Rannala_Karol.nex'
rnseed = 15397
ncycles = 21000
samplefreq = 10
lsweight = 100
external_prior_mean = 0.1
internal_prior_means = [10**x for x in range(-4,0)]

def analyze(internal_prior_mean):
    analyzer = Phycas()
    analyzer.random_seed = rnseed
    analyzer.default_model = 'hky'
    analyzer.num_rates = 5
    analyzer.external_edgelen_dist = ProbDist.ExponentialDist(1.0/external_prior_mean)
    analyzer.internal_edgelen_dist = ProbDist.ExponentialDist(1.0/internal_prior_mean)
    analyzer.edgelen_dist = ProbDist.ExponentialDist(1.0/internal_prior_mean)
    analyzer.using_hyperprior = False 
    analyzer.allow_polytomies = False
    analyzer.data_file_name = data_file_name
    analyzer.starting_tree_source = 'random'
    analyzer.nchains = 1
    analyzer.ncycles = ncycles
    analyzer.sample_every = samplefreq
    analyzer.report_every = ncycles/20
    analyzer.ls_move_weight = lsweight
    analyzer.mcmc()
    shutil.copyfile('%s.p' % data_file_name, '%.9f_%s.p' % (internal_prior_mean, data_file_name))
    shutil.copyfile('%s.t' % data_file_name, '%.9f_%s.t' % (internal_prior_mean, data_file_name))

if __name__ == '__main__':
    for mean in internal_prior_means:
        analyze(mean)
