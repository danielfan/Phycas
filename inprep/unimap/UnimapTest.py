from phycas import *

model.type                  	= 'hky'
model.edgelen_hyperprior		= None

mcmc.random_seed            	= 19375
mcmc.data_source         		= 'simulated.nex'
# mcmc.starting_tree_source   	= TreeCollection(newick='(1:0.096451,2:0.059414,(3:0.069940,(4:0.065929,((5:0.066647,8:0.097704):0.010055,(6:0.063657,((7:0.046537,9:0.052228):0.007358,10:0.059386):0.014523):0.012032):0.008731):0.009339):0.007606)')
#mcmc.starting_tree_source   	= TreeCollection(newick='(1:0.10000,2:0.15000,(3:0.02500,4:0.15000):0.05000)')
mcmc.starting_tree_source   	= TreeCollection(newick='(1:0.20000,2:0.0100,(3:0.200,4:0.2000):0.2000)')
mcmc.ncycles                	= 10
mcmc.sample_every           	= 1
mcmc.report_every           	= 1
mcmc.tree_scaler_weight     	= 0
mcmc.slice_weight           	= 1
mcmc.debugging              	= True
mcmc.fix_topology 				= False
mcmc.edge_move_weight 			= 1
mcmc.edge_move_lambda 			= 0.05
mcmc.out.log.mode				= REPLACE
mcmc.out.params.mode			= REPLACE
mcmc.out.trees.mode				= REPLACE

unimap.mapping_move_weight    	= 1
unimap.unimap_nni_move_weight 	= 1

if False:
	import sys,os
	if os.path.basename(sys.executable) == 'python_d.exe':
		raw_input('debug stop')

unimap()

