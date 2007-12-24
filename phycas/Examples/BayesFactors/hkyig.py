from phycas import *
import shutil

data_file_name = '../../Tests/Data/green.nex'
rnseed = 75391
ncycles = 11000
samplefreq = 100
lsweight = 100

def analyze():
    analyzer = Phycas()
    analyzer.outfile_prefix = 'ps_hkyig'
    analyzer.default_model = 'hky'
    analyzer.num_rates = 4
    analyzer.estimate_pinvar = True
    analyzer.random_seed = rnseed
    analyzer.edgelen_dist = ProbDist.ExponentialDist(10.0)
    analyzer.using_hyperprior = False 
    analyzer.allow_polytomies = False
    analyzer.data_file_name = data_file_name
    analyzer.starting_tree_source = 'random'
    analyzer.is_standard_heating = False
    analyzer.adapt_first = 10
    analyzer.verbose = True
    analyzer.tree_scaler_weight = 1
    analyzer.slice_weight = 1
    analyzer.use_flex_model = False
    analyzer.ps_toward_posterior = False
    analyzer.ps_burnin = 0
    analyzer.ps_Q = ncycles
    analyzer.ps_nbetaincr = 100
    analyzer.sample_every = samplefreq
    analyzer.report_every = ncycles/20
    analyzer.ls_move_weight = lsweight
    analyzer.pathsampling()

if __name__ == '__main__':
    analyze()
