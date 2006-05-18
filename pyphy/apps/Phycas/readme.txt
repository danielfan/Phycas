27-Dec-2005
-----------
Phycas is intended to (eventually) be the standard-bearer for the PyPhy 
library. It is a command line Python application that uses the PyPhy
library to perform an MCMC analysis on the sequences in a supplied nexus
data file. It implements the JC and HKY models only, and updates the
topology, edge legths, and model parameters (base frequencies and kappa).
The samples are stored in MrBayes-like *.t and *.p files that can be
input into Tracer for summarization.

24-Jan-2005
-----------
Timings for release 0.1.x using seed 13579:

Description                     Evals      Seconds     Evals/sec     LOCAL     TBR    Other
-------------------------------------------------------------------------------------------
MrBayes 3.1.1                   500001     24.98400    20012.84822  113837  340939    45224
                                499998     25.03100    19975.15081  113263  341453    45284
                                           ------------------------------------------------
                                mean       25.00750    19993.99952  113550  341196    45254
                                                                   |---- 90.9% ---|

Release 0.1.2 (rel_0_1_2, 24-Jan-2006):
---------------------------------------
  Changed seed to 13579, metropolis weight to 300, ncycles to 1500, sample_every to 15,
  and adapt_first to 60. These settings yield runs that are similar to MrBayes 3.1.1 in
  terms of the proportion of topology moves vs. substitution model parameter updates,
  and in terms of the overall number of likelihood evaluations for a 4-taxon dataset.
  
    Evals      Seconds     Evals/sec     LOCAL     TBR    Other
    -----------------------------------------------------------
    503182     23.84065    21106.05095  450000       0    53182
    503182     23.72691    21207.22843  450000       0    53182
               ------------------------------------------------
    mean       23.78378    21156.63969                         
    speedup    5.1%        5.8%        |---- 89.4% ---|        

Release 0.1.3 (rel_0_1_3, 24-Jan-2006):
---------------------------------------
  Added MCMCChainManager::partialEdgeLenPrior function to simplify calculation of partial 
  edge length priors for moves. Not sure why the slight slowdown was observed, but this
  simplified version is much cleaner so I decided to keep it anyway.

	Evals      Seconds     Evals/sec     LOCAL     TBR    Other
	-----------------------------------------------------------
	503182     24.67662    20391.04492  450000       0    53182
	503182     24.80555    20285.06094  450000       0    53182
	-----------------------------------------------------------
	mean       24.74109    20338.05293                      
	speedup    -3.9%       -3.9%

Release 0.5.0 (rel_0_5_0, 20-Feb-2006):
---------------------------------------
Comparison of PAUP* (MLEs), MrBayes 3.1 and Phycas on nyldna4.nex for HKY+Gamma model
Phycas: 
	data_file_name = '../../pyphy/nyldna4.nex'
	starting_tree_source = 'random'
	ncycles = 10000 (long run)
	ncycles = 1500 (short run)
	sample_every = 10 (long run)
	sample_every = 30 (short run)
	adapt_first = 100
	random_seed = '13579'
	using_hky = True
	num_rates = 4
	using_hyperprior = True
	edgelen_prior_mean = 0.1
	verbose = True
	metropolis_weight = 300
	slice_weight = 1
	gg_do = False
MrBayes:
    set autoclose=yes;
    lset nst=2 rates=gamma;
    mcmc nruns=1 nchains=1 ngen=500000 samplefreq=10000 printfreq=10000;
PAUP:   
    set criterion=likelihood;
    lset nst=2 variant=hky basefreq=estimate tratio=estimate rates=gamma shape=estimate;
    alltrees;
    lscores 1;
    describe 1 / brlens;

            Phycas         Phycas     
            long run       short run               MrBayes                 PAUP (MLEs)             
---------------------------------------------------------------------------------------
lnL         -7054.094       -7053.898 (0.369, 46)  -7054.066 (0.324, 46)   -7049.01129         
TL          0.262           0.264 (0.002, 46)      0.261 (0.002, 35)       0.25872  
piA         0.288           0.286 (0.001, 44)      0.287 (0.001, 46)       0.287612 
piC         0.175           0.176 (0.001, 32)      0.177 (0.001, 46)       0.175021 
piG         0.205           0.205 (0.001, 46)      0.204 (0.001, 46)       0.204780 
piT         0.332           0.333 (0.001, 46)      0.333 (0.001, 46)       0.332586 
exp. ratio   ---             ---                    ---                    0.986465 
kappa       2.101           2.102 (0.027, 46)      2.075 (0.003, 46)       2.105420 
shape       0.218           0.205 (0.006, 46)      0.23  (0.008, 46)       0.215132 
hyper       0.22            0.229 (0.015, 35)       ---                      ---    
cycles      10,000          1,500                  500,000                   ---        
samples     1,000           50                     50                        ---    
nevals      3,418,862       515,487                500,000                   ---    
secs        441.65931       66.74547                 ---                

16-May-2006
-----------

Comparison of MrBayes and Phycas running 4-taxon data set (nyldna4.nex)

Description       Evals      Seconds     Evals/sec     
----------------------------------------------------
Phycas (SVN 41)   1009910    46.76983    21593.19246 <- 7.7% faster
MrBayes 3.1.1      999993    49.85900    20056.41910

Phycas breakdown of likelihood evals:
-------------------------------------
  3000*300         = 900000 evals. (89.1%) devoted to Larget-Simon moves
  1009910 - 900000 = 109910 evals. (10.9%) devoted to parameter updates

MrBayes breakdown:
------------------
  0.6817*999993 =  681695.2281 evals. (68.17%) devoted to TBR
  0.2273*999993 =  227298.4089 evals. (22.73%) devoted to LOCAL
                                       90.90%  devoted to topology moves
                                       
  0.0455*999993 =   45499.6815 evals. (4.55%) devoted to updating kappa
  0.0455*999993 =   45499.6815 evals. (4.55%) devoted to updating state freqs
                                       9.10%  devoted to parameter updates

Comparison of MrBayes and Phycas running 10-taxon data set (green.nex)

Description       Evals      Seconds     Evals/sec    moves    params
----------------------------------------------------------------------
Phycas (SVN 41)  1011178    311.67920   3244.29093   89.01%    10.99%
MrBayes 3.1.1     999993    199.125     5021.93597   90.90%     9.10% <- 54.8% faster
                                       