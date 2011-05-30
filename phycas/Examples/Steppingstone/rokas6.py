import sys
from phycas import *
testing_sample = '-t' in sys.argv
larget_simon = '-l' in sys.argv
short = '-e' in sys.argv
skip_regular = '-s' in sys.argv
regular_mcmc = '-m' in sys.argv
debugging = '-d' in sys.argv

if testing_sample:
    skip_regular = True
try:
    jobid = int(sys.argv[1])
    setMasterSeed(13*jobid)
except:
    jobid = 106
    setMasterSeed(13*106)

model.type = 'gtr'
model.state_freqs = [0.25, 0.25, 0.25, 0.25]
model.state_freq_prior = Dirichlet((1.0,1.0,1.0,1.0))
model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
model.relrates = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0]
model.num_rates = 4
model.gamma_shape = 0.5
model.gamma_shape_prior = Exponential(1.0)
model.pinvar_model = False
model.pinvar = 0.5
model.pinvar_prior = Beta(1.0, 1.0)

model.edgelen_hyperprior = None
model.edgelen_prior = Exponential(10.0)

model.update_freqs_separately = False
model.update_relrates_separately = False

mcmc.data_source = 'rokas6.nex'
mcmc.out.log = 'output.%d.txt' % (jobid,)
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = 'trees.%d' % (jobid,)
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = 'params.%d' % (jobid,)
mcmc.out.params.mode = REPLACE
if testing_sample:
    mcmc.ncycles = 100000
    if larget_simon:
        mcmc.draw_directly_from_prior = False
    else:
        mcmc.draw_directly_from_prior = True
elif debugging:
    mcmc.ncycles = 10
elif short:
    mcmc.ncycles = 100
else:
    if regular_mcmc:
        mcmc.ncycles = 10000
    else:
        mcmc.ncycles = 5000

if skip_regular:
    ss.xcycles = -(mcmc.ncycles - 1)
    mcmc.sample_every = 1
else:
    mcmc.sample_every = 10
mcmc.report_every = 1000

# take tiny steps when updating freqs and GTR relative rates
mcmc.state_freq_psi = 3000.0
mcmc.rel_rate_psi   = 3000.0

mcmc.draw_directly_from_prior = True

if not regular_mcmc and jobid < 106:
    mcmc.fix_topology = True
else:
    mcmc.fix_topology = False
    ss.override_fixed_topology_restriction = True

mcmc.allow_polytomies = False

tree_descriptions = [
    '(1:0.063233,(((2:0.050563,3:0.111184):0.003100,4:0.144445):0,5:0.177951):0.011630,6:0.919418)',
    '(1:0.074357,((2:0.050198,3:0.110953):0.002985,4:0.144097):0,(5:0.106388,6:0.867352):0.070752)',
    '(1:0.074815,(((2:0.050493,3:0.111313):0.003187,4:0.144611):0,6:0.928366):0,5:0.178143)',
    '(1:0.074540,((2:0.050306,3:0.111156):0.003051,(4:0.109331,6:0.900682):0.034958):0,5:0.177828)',
    '(1:0.074794,(((2:0.050471,3:0.111299):0.002801,6:0.928186):0.000426,4:0.144596):0,5:0.178127)',
    '(1:0.075022,(((2:0.043865,6:0.924661):0.007145,3:0.111472):0.002432,4:0.144724):0,5:0.178280)',
    '(1:0.074817,((2:0.050762,(3:0.087645,6:0.910727):0.023877):0.002444,4:0.144828):0,5:0.178350)',
    '(1:0.057905,((2:0.043478,3:0.115217):0.000867,(4:0.090512,5:0.125965):0.065361):0.009970,6:0.900067)',
    '(1:0.067400,(2:0.042973,3:0.115552):0.000724,((4:0.089846,5:0.125021):0.022974,6:0.875748):0.044066)',
    '(1:0.067446,(2:0.043040,3:0.115417):0.000676,(4:0.090079,(5:0.082423,6:0.852166):0.043347):0.065513)',
    '(1:0.067517,(2:0.043057,3:0.115653):0.000692,((4:0.080599,6:0.880680):0.010430,5:0.125499):0.065846)',
    '(1:0.067834,((2:0.043517,3:0.115195):0.000988,6:0.907187):0,(4:0.090525,5:0.125981):0.065315)',
    '(1:0.068273,((2:0.038341,6:0.904676):0.005618,3:0.115405):9.585e-05,(4:0.090544,5:0.126029):0.065619)',
    '(1:0.068011,(2:0.043602,(3:0.087372,6:0.885644):0.028060):0.000101,(4:0.090343,5:0.125767):0.065833)',
    '(1:0.063231,(((2:0.050563,3:0.111184):0.003100,5:0.177951):0,4:0.144445):0.011632,6:0.919416)',
    '(1:0.074540,((2:0.050306,3:0.111156):0.003051,5:0.177829):0,(4:0.109330,6:0.900681):0.034959)',
    '(1:0.074815,(((2:0.050493,3:0.111313):0.003187,5:0.178143):0,6:0.928366):0,4:0.144611)',
    '(1:0.074357,((2:0.050198,3:0.110953):0.002985,(5:0.106388,6:0.867352):0.070752):0,4:0.144097)',
    '(1:0.074794,(((2:0.050471,3:0.111299):0.002801,6:0.928186):0.000426,5:0.178127):0,4:0.144596)',
    '(1:0.075022,(((2:0.043865,6:0.924661):0.007145,3:0.111472):0.002432,5:0.178280):0,4:0.144724)',
    '(1:0.074817,((2:0.050762,(3:0.087645,6:0.910727):0.023877):0.002444,5:0.178351):0,4:0.144828)',
    '(1:0.063693,(((2:0.052361,5:0.179130):0,3:0.112339):0,4:0.145515):0.012086,6:0.921096)',
    '(1:0.075470,((2:0.052086,5:0.178960):0,3:0.112283):0,(4:0.109782,6:0.902126):0.035520)',
    '(1:0.075759,(((2:0.052339,5:0.179361):1.042e-08,3:0.112500):0,6:0.930547):0,4:0.145714)',
    '(1:0.075728,(((2:0.052232,5:0.179297):0,6:0.930358):0.000126,3:0.112481):1.084e-08,4:0.145705)',
    '(1:0.075280,((2:0.051947,(5:0.106825,6:0.868609):0.071391):0,3:0.112055):0,4:0.145094)',
    '(1:0.075789,(((2:0.044614,6:0.925480):0.007793,5:0.179184):0,3:0.112372):0,4:0.145549)',
    '(1:0.075586,((2:0.052178,5:0.179248):0,(3:0.088033,6:0.911506):0.024381):0,4:0.145648)',
    '(1:0.060808,((2:0.047574,(3:0.095699,5:0.164378):0.022181):0,4:0.149437):0.010662,6:0.916508)',
    '(1:0.071008,(2:0.047151,(3:0.095241,5:0.163800):0.022645):0,(4:0.110889,6:0.893722):0.038377)',
    '(1:0.071486,((2:0.047620,(3:0.095699,5:0.164402):0.022231):0,6:0.924332):0,4:0.149497)',
    '(1:0.071464,((2:0.041457,6:0.920785):0.006184,(3:0.095724,5:0.164418):0.022179):0,4:0.149458)',
    '(1:0.070825,(2:0.046874,((3:0.094850,5:0.163700):0.005159,6:0.914290):0.018844):0,4:0.150206)',
    '(1:0.070713,(2:0.046853,(3:0.095095,(5:0.101069,6:0.862937):0.062427):0.022784):0,4:0.149483)',
    '(1:0.070917,(2:0.046984,((3:0.079818,6:0.904355):0.015832,5:0.164109):0.022967):0,4:0.150080)',
    '(1:0.063696,(((2:0.052361,4:0.145515):0,3:0.112339):0,5:0.179130):0.012084,6:0.921098)',
    '(1:0.075280,((2:0.051947,4:0.145094):0,3:0.112055):0,(5:0.106825,6:0.868609):0.071391)',
    '(1:0.075595,(((2:0.052305,4:0.145710):1.003e-08,3:0.112482):0.000186,6:0.930303):0,5:0.179236)',
    '(1:0.075759,(((2:0.052339,4:0.145714):0,6:0.930547):0,3:0.112500):1.151e-08,5:0.179361)',
    '(1:0.075470,((2:0.052086,(4:0.109783,6:0.902126):0.035519):0,3:0.112283):0,5:0.178960)',
    '(1:0.075789,(((2:0.044613,6:0.925480):0.007793,4:0.145549):0,3:0.112372):0,5:0.179184)',
    '(1:0.075586,((2:0.052178,4:0.145648):0,(3:0.088033,6:0.911506):0.024381):0,5:0.179248)',
    '(1:0.063363,(((2:0.052024,4:0.145255):0,5:0.178825):0.001494,3:0.111432):0.011347,6:0.920378)',
    '(1:0.074289,(((2:0.051887,4:0.145340):0,5:0.178917):0.000922,6:0.928763):0.001136,3:0.111260)',
    '(1:0.073957,((2:0.051544,4:0.144787):0,(5:0.106696,6:0.867857):0.071163):0.001843,3:0.110951)',
    '(1:0.074360,(((2:0.051924,4:0.145366):0,6:0.929399):0,5:0.178946):0.001967,3:0.111322)',
    '(1:0.074120,((2:0.051677,(4:0.109675,6:0.901281):0.035313):0,5:0.178580):0.001887,3:0.111154)',
    '(1:0.074428,(((2:0.044424,6:0.924542):0.007581,4:0.145212):0,5:0.178785):0.001911,3:0.111226)',
    '(1:0.074525,((2:0.051848,4:0.145396):0,5:0.178952):0.001472,(3:0.087865,6:0.910865):0.023664)',
    '(1:0.057994,((2:0.043661,(4:0.090530,5:0.125984):0.065562):0.000529,3:0.115176):0.010048,6:0.900250)',
    '(1:0.067787,((2:0.043596,(4:0.090549,5:0.125996):0.065497):0.000257,6:0.907337):0.000699,3:0.114987)',
    '(1:0.067309,(2:0.042940,((4:0.089860,5:0.125024):0.022955,6:0.875722):0.044232):0.000840,3:0.115321)',
    '(1:0.067336,(2:0.042994,(4:0.090089,(5:0.082429,6:0.852107):0.043343):0.065649):0.000819,3:0.115181)',
    '(1:0.067403,(2:0.043005,((4:0.080615,6:0.880612):0.010424,5:0.125499):0.065986):0.000845,3:0.115409)',
    '(1:0.067804,((2:0.038170,6:0.904443):0.005482,(4:0.090549,5:0.126016):0.065544):0.000832,3:0.115039)',
    '(1:0.067765,(2:0.043423,(4:0.090350,5:0.125762):0.065787):0.000500,(3:0.087332,6:0.885527):0.027899)',
    '(1:0.063363,(((2:0.052024,5:0.178825):0,4:0.145255):0.001494,3:0.111432):0.011347,6:0.920378)',
    '(1:0.074289,(((2:0.051887,5:0.178917):0,4:0.145340):0.000922,6:0.928763):0.001136,3:0.111260)',
    '(1:0.074120,((2:0.051677,5:0.178580):0,(4:0.109674,6:0.901280):0.035315):0.001887,3:0.111154)',
    '(1:0.074360,(((2:0.051924,5:0.178946):0,6:0.929399):0,4:0.145366):0.001967,3:0.111322)',
    '(1:0.073957,((2:0.051544,(5:0.106696,6:0.867857):0.071163):0,4:0.144787):0.001843,3:0.110951)',
    '(1:0.074428,(((2:0.044424,6:0.924542):0.007581,5:0.178785):0,4:0.145212):0.001911,3:0.111226)',
    '(1:0.074525,((2:0.051848,5:0.178952):0,4:0.145396):0.001472,(3:0.087864,6:0.910865):0.023665)',
    '(1:0.060804,((2:0.047574,4:0.149437):0,(3:0.095698,5:0.164378):0.022181):0.010665,6:0.916505)',
    '(1:0.071486,((2:0.047620,4:0.149497):0,6:0.924332):0,(3:0.095699,5:0.164402):0.022231)',
    '(1:0.071008,(2:0.047151,(4:0.110889,6:0.893722):0.038377):0,(3:0.095241,5:0.163800):0.022645)',
    '(1:0.071464,((2:0.041457,6:0.920785):0.006184,4:0.149458):0,(3:0.095724,5:0.164418):0.022179)',
    '(1:0.070825,(2:0.046874,4:0.150205):0,((3:0.094850,5:0.163700):0.005159,6:0.914290):0.018844)',
    '(1:0.070713,(2:0.046853,4:0.149483):0,(3:0.095095,(5:0.101069,6:0.862937):0.062427):0.022784)',
    '(1:0.070917,(2:0.046984,4:0.150080):0,((3:0.079818,6:0.904355):0.015832,5:0.164109):0.022967)',
    '(1:0.059815,((2:0.045878,(3:0.094333,4:0.128695):0.026269):0,5:0.184525):0.010379,6:0.916429)',
    '(1:0.069493,(2:0.045265,(3:0.093638,4:0.127869):0.026757):0,(5:0.108444,6:0.859143):0.075217)',
    '(1:0.070205,((2:0.045947,(3:0.094335,4:0.128705):0.026289):1.726e-05,6:0.923962):0,5:0.184545)',
    '(1:0.070177,((2:0.040154,6:0.920734):0.005805,(3:0.094360,4:0.128718):0.026260):0,5:0.184558)',
    '(1:0.069710,(2:0.045347,((3:0.093827,4:0.128295):0.011767,6:0.916340):0.015725):0,5:0.185179)',
    '(1:0.069654,(2:0.045333,(3:0.094330,(4:0.099976,6:0.897087):0.028791):0.026431):0,5:0.184881)',
    '(1:0.069762,(2:0.045423,((3:0.077232,6:0.906658):0.017476,4:0.128874):0.026368):0,5:0.185150)',
    '(1:0.059803,((2:0.045878,5:0.184524):0,(3:0.094332,4:0.128694):0.026269):0.010388,6:0.916421)',
    '(1:0.070219,((2:0.045951,5:0.184555):0,6:0.923980):0,(3:0.094337,4:0.128704):0.026290)',
    '(1:0.069493,(2:0.045265,(5:0.108444,6:0.859143):0.075217):0,(3:0.093638,4:0.127869):0.026757)',
    '(1:0.070177,((2:0.040153,6:0.920733):0.005806,5:0.184558):0,(3:0.094360,4:0.128718):0.026260)',
    '(1:0.069710,(2:0.045347,5:0.185179):0,((3:0.093827,4:0.128295):0.011766,6:0.916340):0.015725)',
    '(1:0.069654,(2:0.045333,5:0.184881):0,(3:0.094330,(4:0.099976,6:0.897087):0.028791):0.026431)',
    '(1:0.069762,(2:0.045423,5:0.185150):0,((3:0.077233,6:0.906658):0.017476,4:0.128874):0.026368)',
    '(1:0.054722,(2:0.038249,((3:0.094875,4:0.118170):0.001423,5:0.152567):0.042797):0.010251,6:0.905463)',
    '(1:0.064851,(2:0.033273,6:0.910143):0.005159,((3:0.094907,4:0.118200):0.001420,5:0.152617):0.042764)',
    '(1:0.064670,2:0.037944,(((3:0.094908,4:0.116963):0.001288,5:0.150910):0.010304,6:0.887106):0.034088)',
    '(1:0.064518,2:0.037928,((3:0.094269,4:0.117046):0.001303,(5:0.094188,6:0.844201):0.056744):0.043309)',
    '(1:0.064692,2:0.037985,(((3:0.094898,4:0.116909):0.001411,6:0.891794):0.000159,5:0.150834):0.044123)',
    '(1:0.064635,2:0.037977,((3:0.094883,(4:0.094459,6:0.874990):0.022987):0.000852,5:0.151112):0.043887)',
    '(1:0.064681,2:0.038007,(((3:0.079269,6:0.880449):0.016061,4:0.117171):0.000825,5:0.150945):0.044189)',
    '(1:0.054336,(2:0.037661,(3:0.088978,(4:0.087785,5:0.124253):0.043176):0.034745):0.010151,6:0.891157)',
    '(1:0.064354,(2:0.032906,6:0.895826):0.004951,(3:0.089007,(4:0.087810,5:0.124298):0.043186):0.034714)',
    '(1:0.064219,2:0.037421,((3:0.088324,(4:0.087710,5:0.123810):0.042532):0.008307,6:0.878226):0.027750)',
    '(1:0.064233,2:0.037411,(3:0.088134,((4:0.087210,5:0.123040):0.016078,6:0.862565):0.027939):0.035216)',
    '(1:0.064188,2:0.037378,(3:0.088169,(4:0.087427,(5:0.081579,6:0.836677):0.042230):0.042388):0.035303)',
    '(1:0.064282,2:0.037417,(3:0.088256,((4:0.078512,6:0.864481):0.009807,5:0.123561):0.042599):0.035397)',
    '(1:0.064247,2:0.037462,((3:0.073907,6:0.871394):0.014842,(4:0.087702,5:0.123747):0.042651):0.035450)',
    '(1:0.054708,(2:0.038255,((3:0.094665,5:0.152145):0.001379,4:0.118544):0.043033):0.010256,6:0.905462)',
    '(1:0.064842,(2:0.033268,6:0.910138):0.005170,((3:0.094696,5:0.152194):0.001379,4:0.118573):0.042999)',
    '(1:0.064666,2:0.037954,(((3:0.094904,5:0.150679):0.001005,4:0.117342):0.010515,6:0.887266):0.034146)',
    '(1:0.064628,2:0.037979,((3:0.094574,5:0.150757):0.001080,(4:0.094481,6:0.874944):0.023136):0.043960)',
    '(1:0.064670,2:0.037996,(((3:0.094209,5:0.150033):0,6:0.890497):0.002356,4:0.117133):0.043951)',
    '(1:0.064509,2:0.037936,((3:0.093808,(5:0.094172,6:0.843952):0.056187):0.001687,4:0.117313):0.043347)',
    '(1:0.064662,2:0.038006,(((3:0.079095,6:0.879746):0.015347,5:0.150154):0.001972,4:0.117155):0.043899)'
]
if jobid > 105:
    mcmc.starting_tree_source = TreeCollection(newick=tree_descriptions[95])
    ss.refdist_definition_file = 'rokas6.txt'
else:
    mcmc.starting_tree_source = TreeCollection(newick=tree_descriptions[jobid-1])
if testing_sample:
    ss.nbetavals = 2
else:
    ss.nbetavals = 21
if regular_mcmc:
    mcmc()
    '''sump.out.log = 'mcmc.sump.%d.txt' % (jobid,)
    sump.out.log.mode = REPLACE
    sump.file = 'mcmc.params.%d.p' % (jobid,)
    sump()'''
    sumt.burnin = 5
    sumt.out.refdistfile.prefix = "refdist"
    sumt.out.refdistfile.mode = REPLACE
    sumt.trees = 'trees.106.t'
    sumt.out.log = 'mcmc.sumt.%d.txt' % (jobid,)
    sumt.out.log.mode = REPLACE
    sumt()
else:
    ss()

    if testing_sample:
        sumt.trees = 'trees.106.t'
        sumt.out.splits = False
        if larget_simon:
            sumt.out.log = 'trees_from_refdist_larget_simon.%d.txt' % (jobid,)
        else:
            sumt.out.log = 'trees_from_refdist_direct_draw.%d.txt' % (jobid,)
        sumt.out.log.mode = REPLACE
        sumt.tree_credible_prob = 1.0
        sumt()
    else:
        sump.out.log = 'ss.sump.%d.txt' % (jobid,)
        sump.out.log.mode = REPLACE
        sump.file = 'params.%d.p' % (jobid,)
        sump()

