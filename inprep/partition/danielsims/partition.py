from phycas import *

compare_to_daniel = False
data_file_name = 'data186.nex'
master_seed = 16223
partitioned_analysis = True
plus_gamma = True
scubed = True
num_beta_values = 0
estimate_subset_relative_rates = True
if num_beta_values > 0:
    # doing steppingtone sampling
    num_cycles = 100
    num_cycles_per_sample = 1
else:
    # just running plain mcmc
    num_cycles = 22000
    num_cycles_per_sample = 10

if compare_to_daniel:
	v1          = 0.0158816
	v2          = 0.019522
	v3          = 0.0135714
	v4          = 0.238507
	v5          = 0.0880031
	v6          = 0.184429
	v7          = 0.00464573
	v8          = 0.0127538
	v9          = 0.00810269
	v10         = 0.0139133
	v11         = 0.00106932
	v12         = 0.0240725
	v13         = 0.0124062
	v14         = 0.0184587
	v15         = 0.015272
	v16         = 0.00713013
	v17         = 0.129027
	v18         = 0.0249049
	v19         = 0.0780493
	v20         = 0.216845
	v21         = 0.00360375
	v22         = 0.15527
	v23         = 0.00783398
	v24         = 0.0187774
	v25         = 0.000985893
	v26         = 0.162698
	v27         = 0.030222
	v28         = 0.0245678
	v29         = 0.570006
	v30         = 0.583448
	v31         = 0.00323834
	m_1         = 1
	m_2         = 1
	rAC_P1      = 0.157064
	rAG_P1      = 0.172612
	rAT_P1      = 0.175898
	rCG_P1      = 0.176805
	rCT_P1      = 0.165537
	rGT_P1      = 0.152083
	piA_P1      = 0.278393
	piC_P1      = 0.247253
	piG_P1      = 0.23245
	piT_P1      = 0.241905
	alpha_P1    = 302.758
	rAC_P2      = 0.157267
	rAG_P2      = 0.181339
	rAT_P2      = 0.168812
	rCG_P2      = 0.172162
	rCT_P2      = 0.160687
	rGT_P2      = 0.159734
	piA_P2      = 0.257363
	piC_P2      = 0.230761
	piG_P2      = 0.255366
	piT_P2      = 0.256509
	alpha_P2    = 17.4487

setMasterSeed(master_seed)
blob = readFile(data_file_name)
nchar = blob.characters.getMatrix().n_char

tree_descriptions = []
for t in blob.trees:
    tree_descriptions.append(t.newick)
assert(len(tree_descriptions) > 0)

if partitioned_analysis:
    model.type = 'gtr'
    if compare_to_daniel:
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
        if compare_to_daniel:
            model.gamma_shape = alpha_P1
        else:
            model.gamma_shape = 1.0
        model.gamma_shape_prior = Exponential(1.0/100.0)
    else:
        model.num_rates = 1
    m1 = model()
    
    model.type = 'gtr'
    if compare_to_daniel:
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
        if compare_to_daniel:
            model.gamma_shape = alpha_P2
        else:
            model.gamma_shape = 1.0
        model.gamma_shape_prior = Exponential(1.0/100.0)
    else:
        model.num_rates = 1
    m2 = model()
    
    half_way = nchar/2
    if compare_to_daniel:
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

if compare_to_daniel:
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
    mcmc.report_every = 100

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
    rAC_P1      = 0.16666667
    rAG_P1      = 0.16666667
    rAT_P1      = 0.16666667
    rCG_P1      = 0.16666667
    rCT_P1      = 0.16666667
    rGT_P1      = 0.16666667
    piA_P1      = 0.25
    piC_P1      = 0.25
    piG_P1      = 0.25
    piT_P1      = 0.25
    alpha_P1    = 1.0
    rAC_P2      = 0.16666667
    rAG_P2      = 0.16666667
    rAT_P2      = 0.16666667
    rCG_P2      = 0.16666667
    rCT_P2      = 0.16666667
    rGT_P2      = 0.16666667
    piA_P2      = 0.25
    piC_P2      = 0.25
    piG_P2      = 0.25
    piT_P2      = 0.25
    alpha_P2    = 1.0
"""
