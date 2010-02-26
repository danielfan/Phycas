from phycas import *

compute_lnL_for_point = False
data_file_name = 'data186.nex'
master_seed = 18625
partitioned_analysis = True
plus_gamma = True
scubed = True
num_beta_values = 11
estimate_subset_relative_rates = False
reportfreq = 100
if num_beta_values > 0:
    # doing steppingtone sampling
    num_cycles = 2000
    num_cycles_per_sample = 10
else:
    # just running plain mcmc
    num_cycles = 22000
    num_cycles_per_sample = 10

if compute_lnL_for_point:
	v1          = 0.0071183
	v2          = 0.00553953
	v3          = 0.02318026
	v4          = 0.20710501
	v5          = 0.09082643
	v6          = 0.1469528
	v7          = 0.00614295
	v8          = 0.00506076
	v9          = 0.00498758
	v10         = 0.00915324
	v11         = 0.00930544
	v12         = 0.01257619
	v13         = 0.01491869
	v14         = 0.00831953
	v15         = 0.00629958
	v16         = 0.00360241
	v17         = 0.10792657
	v18         = 0.00320366
	v19         = 0.08690338
	v20         = 0.17658233
	v21         = 0.00222961
	v22         = 0.11643803
	v23         = 0.00317141
	v24         = 0.00650089
	v25         = 0.00251004
	v26         = 0.18842601
	v27         = 0.00170542
	v28         = 0.03969147
	v29         = 0.60168904
	v30         = 0.46110498
	v31         = 0.00112575
	m_1         = 0.94288505
	m_2         = 1.05711495
	rAC_P1      = 0.14977019
	rAG_P1      = 0.15797095
	rAT_P1      = 0.17350201
	rCG_P1      = 0.14053361
	rCT_P1      = 0.22226152
	rGT_P1      = 0.15596172
	piA_P1      = 0.29228188
	piC_P1      = 0.24865658
	piG_P1      = 0.24312952
	piT_P1      = 0.21593203
	alpha_P1    = 32.16783073
	rAC_P2      = 0.12683939
	rAG_P2      = 0.15711713
	rAT_P2      = 0.13706274
	rCG_P2      = 0.15962537
	rCT_P2      = 0.22751824
	rGT_P2      = 0.19183712
	piA_P2      = 0.28736845
	piC_P2      = 0.24600072
	piG_P2      = 0.24518242
	piT_P2      = 0.22144841
	alpha_P2    = 23.12501776

setMasterSeed(master_seed)
blob = readFile(data_file_name)
nchar = blob.characters.getMatrix().n_char

tree_descriptions = []
for t in blob.trees:
    tree_descriptions.append(t.newick)
assert(len(tree_descriptions) > 0)

if partitioned_analysis:
    model.type = 'gtr'
    if compute_lnL_for_point:
        model.relrates = [rAC_P1, rAG_P1, rAT_P1, rCG_P1, rCT_P1, rGT_P1]
        model.state_freqs = [piA_P1, piC_P1, piG_P1, piT_P1]
    else:
        model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
        model.state_freqs = [0.25, 0.25, 0.25, 0.25]
    model.update_freqs_separately = False
    model.update_relrates_separately = False
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None
    model.pinvar_model = False
    if plus_gamma:
        model.num_rates = 4
        if compute_lnL_for_point:
            model.gamma_shape = alpha_P1
        else:
            model.gamma_shape = 1.0
        model.gamma_shape_prior = Exponential(1.0/100.0)
    else:
        model.num_rates = 1
    m1 = model()
    
    model.type = 'gtr'
    if compute_lnL_for_point:
        model.relrates = [rAC_P2, rAG_P2, rAT_P2, rCG_P2, rCT_P2, rGT_P2]
        model.state_freqs = [piA_P2, piC_P2, piG_P2, piT_P2]
    else:
        model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
        model.state_freqs = [0.25, 0.25, 0.25, 0.25]
    model.update_freqs_separately = False
    model.update_relrates_separately = False
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None
    model.pinvar_model = False
    if plus_gamma:
        model.num_rates = 4
        if compute_lnL_for_point:
            model.gamma_shape = alpha_P2
        else:
            model.gamma_shape = 1.0
        model.gamma_shape_prior = Exponential(1.0/100.0)
    else:
        model.num_rates = 1
    m2 = model()
    
    half_way = nchar/2
    if compute_lnL_for_point:
        partition.subset_relrates = [m_1,m_2]
    else:       
        partition.subset_relrates = [1.0,1.0]
    if estimate_subset_relative_rates:
        partition.fix_subset_relrates = False
    else:
        partition.fix_subset_relrates = True
    partition.addSubset(subset(1,half_way),         m1, 'first_half')
    partition.addSubset(subset(half_way + 1,nchar), m2, 'last_half')
    partition()
    
else:
    model.type = 'gtr'
    model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
    model.update_freqs_separately = False
    model.update_relrates_separately = False
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None
    model.pinvar_model = False
    if plus_gamma:
        model.num_rates = 4
        model.gamma_shape = 1.0
        model.gamma_shape_prior = Exponential(1.0/100.0)
    else:
        model.num_rates = 1

if num_beta_values == 0:
    pfx = '%s_%s_%s_%s_%d' % (data_file_name,plus_gamma and 'plusG' or 'minusG',partitioned_analysis and 'partitioned' or 'unpartitioned',estimate_subset_relative_rates and 'ssrrvar' or 'ssrrfix',master_seed)
else:
    pfx = '%s_%s_%s_%s_%s_%d_%d' % (data_file_name,plus_gamma and 'plusG' or 'minusG',partitioned_analysis and 'partitioned' or 'unpartitioned',scubed and 'sss' or 'ss',estimate_subset_relative_rates and 'ssrrvar' or 'ssrrfix',num_beta_values,master_seed)

if compute_lnL_for_point:
    like.preorder_edgelens = [v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18, v19, v20, v21, v22, v23, v24, v25, v26, v27, v28, v29, v30, v31]
    like.data_source = blob.characters
    like.tree_source = TreeCollection(newick=tree_descriptions[0])
    lnl = like()
    print 'lnL =',lnl
else:
    mcmc.data_source = blob.characters
    mcmc.out.log.prefix = pfx
    mcmc.out.log.mode = REPLACE
    mcmc.out.trees.prefix = pfx
    mcmc.out.trees.mode = REPLACE
    mcmc.out.params.prefix = pfx
    mcmc.out.params.mode = REPLACE
    mcmc.nchains = 1
    mcmc.adapt_first = 2
    mcmc.ncycles = num_cycles
    mcmc.sample_every = num_cycles_per_sample
    mcmc.report_every = reportfreq;

    mcmc.state_freq_weight = 10
    if scubed:
        mcmc.state_freq_psi = 500.0
        mcmc.state_freq_psi0 = 500.0
    else:
        mcmc.state_freq_psi = 500.0
        mcmc.state_freq_psi0 = 1.0
        
    mcmc.rel_rate_weight = 10
    if scubed:
        mcmc.rel_rate_psi = 500.0
        mcmc.rel_rate_psi0 = 500.0
    else:
        mcmc.rel_rate_psi = 500.0
        mcmc.rel_rate_psi0 = 1.0

    mcmc.subset_relrates_weight = (estimate_subset_relative_rates and 1 or 0)
    if scubed:
        mcmc.subset_relrates_psi = 300.0
        mcmc.subset_relrates_psi0 = 300.0
    else:
        mcmc.subset_relrates_psi = 300.0
        mcmc.subset_relrates_psi0 = 1.0

    mcmc.tree_scaler_weight = 1
    if scubed:
        mcmc.tree_scaler_lambda = 0.5
        mcmc.tree_scaler_lambda0 = 0.5
    else:
        mcmc.tree_scaler_lambda = 0.5
        mcmc.tree_scaler_lambda0 = 1.0

    mcmc.slice_weight = 1

    mcmc.starting_tree_source = TreeCollection(newick=tree_descriptions[0]) 
    mcmc.debugging = False

    #mcmc.ls_move_weight = 100
    mcmc.edge_move_weight = 50
    if scubed:
        mcmc.edge_move_lambda = 1.0
        mcmc.edge_move_lambda0 = 1.0
    else:
        mcmc.edge_move_lambda = 1.0
        mcmc.edge_move_lambda0 = 2.0

    mcmc.fix_topology = True
    if num_beta_values == 0:
        mcmc()
    else:
        ss.nbetavals = num_beta_values
        if scubed:
            ss.shape1 = 1.0
            ss.shape2 = 1.0
            ss.scubed = True
        else:
            ss.shape1 = 0.3
            ss.shape2 = 1.0
            ss.scubed = False
        ss()

    sump.file = pfx + '.p'
    sump.out.log = pfx + '.sump.txt'
    sump.out.log.mode = REPLACE
    if num_beta_values == 0:
        sump.burnin = 200
    else:
        sump.burnin = 1
    sump()

print 'Run %s has finished.' % pfx

ignore = """
    # true values used for simulation
	v1          = 0.0155
	v2          = 0.00948
	v3          = 0.01356
	v4          = 0.21492
	v5          = 0.11729
	v6          = 0.16423
	v7          = 0.00681
	v8          = 0.00833
	v9          = 0.00515
	v10         = 0.0088
	v11         = 0.00445
	v12         = 0.02349
	v13         = 0.00614
	v14         = 0.00789
	v15         = 0.01202
	v16         = 0.00988
	v17         = 0.14025
	v18         = 0.01098
	v19         = 0.06575
	v20         = 0.21438
	v21         = 0.00087
	v22         = 0.16026
	v23         = 0.00928
	v24         = 0.0089
	v25         = 0.00471
	v26         = 0.18732
	v27         = 0.01637
	v28         = 0.03656
	v29         = 0.62031
	v30         = 0.59461
	v31         = 0.00323
	m_1         = 1.0
	m_2         = 1.0
	rAC_P1      = 0.164102
	rAG_P1      = 0.166287
	rAT_P1      = 0.175143
	rCG_P1      = 0.144324
	rCT_P1      = 0.189252
	rGT_P1      = 0.160891
	piA_P1      = 0.276584
	piC_P1      = 0.256870
	piG_P1      = 0.241268
	piT_P1      = 0.225278
	alpha_P1    = 6.990495
	rAC_P2      = 0.164102
	rAG_P2      = 0.166287
	rAT_P2      = 0.175143
	rCG_P2      = 0.144324
	rCT_P2      = 0.189252
	rGT_P2      = 0.160891
	piA_P2      = 0.276584
	piC_P2      = 0.256870
	piG_P2      = 0.241268
	piT_P2      = 0.225278
	alpha_P2    = 6.990495
"""
