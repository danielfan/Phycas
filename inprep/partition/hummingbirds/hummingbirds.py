from phycas import *

scheme          = 1
just_testing    = False
fixed_topology  = False
user_tree_def   = '(1:0.04640953,((((2:0.10083239,(((74:0.06573636,80:0.06079629):0.02499439,(((76:0.04211028,(77:0.02028882,78:0.04075987):0.02200350):0.00258667,79:0.04988411):0.03466390,81:0.08275208):0.00469364):0.01183470,75:0.08269005):0.00647305):0.00737666,(((((((((3:0.14241130,89:0.11968377):0.01588842,(((9:0.12023159,(127:0.12199255,(160:0.13578537,(161:0.01966805,162:0.02135534):0.08645868):0.00834312):0.08307470):0.04293166,71:0.19213915):0.12868863,(((11:0.34624168,(23:0.39285592,50:0.21062105):0.03549638):0.02132359,(109:0.69280942,129:0.34361994):0.01055400):0.04012362,163:0.34133378):0.06635304):0.36777659):0.00515860,((((4:0.04823198,5:0.04727534):0.02489479,(87:0.03615125,88:0.03853889):0.02747223):0.04491776,(((((25:0.04150487,35:0.03548123):0.00886628,(29:0.03521314,(31:0.00107159,32:0.00342959):0.04180869):0.00470227):0.01482175,26:0.07611199):0.00814396,(((27:0.01760838,39:0.01542277):0.01168872,(30:0.03456101,36:0.02628006):0.00447441):0.02636653,40:0.04315070):0.01481749):0.00567593,(((28:0.02236402,33:0.01780789):0.01365659,(34:0.00000277,37:0.00325805):0.04108633):0.03908936,38:0.05660866):0.00473664):0.03095190):0.06414104,(90:0.12859580,91:0.11412686):0.05521239):0.02150577):0.02742153,(((17:0.11455889,(110:0.05852276,(111:0.04430812,112:0.03594057):0.00628836):0.06535513):0.00785913,(104:0.04681275,105:0.04360069):0.05409576):0.03696037,(((21:0.05194478,22:0.06449211):0.03573269,(((94:0.00815936,95:0.01052913):0.02188843,(143:0.03631264,144:0.04108205):0.00567946):0.02430230,119:0.06117142):0.05177471):0.04846384,((72:0.14996305,84:0.13109170):0.01080090,145:0.13842759):0.03028277):0.00692950):0.00895385):0.03284032,((((((((6:0.01045312,7:0.02166827):0.02988570,(123:0.03641825,124:0.02919843):0.01128157):0.01304598,((52:0.00366851,(102:0.00407939,103:0.00661140):0.00437190):0.02478178,(92:0.01447651,93:0.01470186):0.01916399):0.01529024):0.04341607,((8:0.05544777,142:0.05339483):0.02665006,((((((61:0.01626012,68:0.03126787):0.00835518,(118:0.01677689,147:0.01422693):0.00833725):0.00509386,151:0.02913949):0.01611813,((((67:0.02391212,154:0.02757016):0.00321202,152:0.02409423):0.02120587,(69:0.02774941,(70:0.03176441,108:0.01934726):0.00758293):0.02254480):0.00188634,(153:0.02194580,155:0.01502917):0.02249719):0.00442383):0.00876045,(148:0.04859432,150:0.04963204):0.00751926):0.01147687,(146:0.01641112,149:0.01311977):0.03974895):0.02509894):0.02229126):0.01433085,((44:0.05234383,66:0.05792028):0.02230819,(((130:0.04565417,132:0.03955999):0.00647512,131:0.04771428):0.00081837,133:0.06391111):0.02030140):0.02139708):0.00694373,((120:0.02009009,122:0.01752456):0.01563402,121:0.02422447):0.08473443):0.01159870,(((((((10:0.01067187,14:0.00615424):0.00748994,15:0.00284290):0.00181056,16:0.01856089):0.00544178,(134:0.01819503,135:0.02531030):0.00834799):0.02108478,(140:0.00850415,141:0.01403871):0.02628805):0.01095648,((18:0.01609805,51:0.01604822):0.00159213,(128:0.01251850,(136:0.01340645,137:0.01520677):0.00430941):0.00276681):0.02793322):0.05492404,(((42:0.06318201,73:0.09676024):0.01163173,96:0.08365538):0.01178319,(62:0.08766184,(63:0.00036942,64:0.00000288):0.08445080):0.00793953):0.00799359):0.01671256):0.00500329,41:0.13628057):0.01929874):0.00759159,((((12:0.02828814,13:0.03660278):0.09099250,((57:0.03022991,58:0.02734952):0.04691793,(106:0.02799337,107:0.03392233):0.05782938):0.05524737):0.00546959,((((((19:0.04622391,(59:0.02149172,60:0.03668707):0.01391867):0.02547258,48:0.03876160):0.00399554,(((43:0.03624041,(47:0.04386983,125:0.03938019):0.00575737):0.00519697,126:0.04290819):0.00543736,(((53:0.02726146,56:0.02713151):0.01541637,54:0.03356577):0.00531865,55:0.03816579):0.00837725):0.01315062):0.00338499,(45:0.00958633,46:0.01173891):0.05286142):0.01301322,((156:0.02836968,157:0.02861064):0.03919395,164:0.06518590):0.00854136):0.01321227,24:0.13824418):0.00912199):0.00367894,(82:0.04829424,83:0.03217149):0.07438489):0.01328676):0.01599822,((85:0.01278588,86:0.01138009):0.07440001,(((97:0.02527642,98:0.02829281):0.03410401,100:0.05435568):0.01353992,99:0.05523894):0.02211434):0.02294889):0.00454520,(158:0.01061004,159:0.01626976):0.09371373):0.00187559,(65:0.09113909,((113:0.04322415,117:0.03926446):0.04090617,((114:0.04382225,116:0.04255993):0.01355817,115:0.04626452):0.02645982):0.03229439):0.00555639):0.00428728):0.00308157,(20:0.10139505,101:0.08525093):0.01160755):0.01234313,(138:0.02869866,139:0.02283562):0.05044033):0.03078831,49:0.05630108);'

fnprefix = 'hbirds%d' % scheme

# McGuire, JA, CC Witt, DL Altshuler and JV Remsen. 2007.
# Phylogenetic Systematics and Biogeography of Hummingbirds: 
# Bayesian and Maximum Likelihood Analyses of Partitioned Data 
# and Selection of an Appropriate Partitioning Strategy.
# Systematic Biology 56(5): 837-856.

def abort(message):
    print 'Aborting...'
    print 'Reason:',message
    import sys
    sys.exit(0)

# scheme  subsets  description from Table 2 of McGuire et al. (2007)
#   1        1     Unpartitioned: (All nucleotide positions: GTR+I+G)
#   2        2     Partitions: (mtDNA, nuc: GTR+I+G) 
#   3        4     Partitions: (mtpos1, mtpos2, mtpos3, nuc: GTR+I+G)
#   4        5     Partitions A: (AK1, BFib: GTR+G), (mtpos1, mtpos2, mtpos3: GTR+I+G)
#   5        5     Partitions B: (tRNA, mtpos1, mtpos2, mtpos3, nuc: GTR+I+G) 
#   6        6     Partitions A: (tRNA, mtpos1, mtpos3, ND2pos2, ND4pos2, nuc: GTR+I+G) 
#   7        6     Partitions B: (AK1, BFib: GTR+G), (mtpos1, mtpos3, ND2pos2, ND4pos2: GTR+I+G) 
#   8        8     Partitions: (tRNA, ND2pos1, ND2pos2, ND2pos3, ND4pos1,ND4pos2, ND4pos3, nuc: GTR+I+G) 
#   9        9     Partitions: (tRNA, ND2pos1, ND2pos2, ND2pos3, ND4pos1, ND4pos2, ND4pos3: GTR+I+G), (BFib, AK1: GTR+G)
def createModelForScheme(scheme):
    # mtDNA sites
    tRNA = subset(1,6) + subset(1048,1074) + subset(3903,4122) 
    ND2 = subset(7,1047)
    ND4 = subset(3208,3902)
    mtDNA = tRNA + ND2 + ND4
    
    # nuclear sites
    AK1 = subset(1075,1815)
    bfib = subset(1816,3207) 
    nuc = AK1 + bfib
    
    # MLEs based on published tree
    relrate_MLEs    = [0.20301, 6.91377, 0.42937, 0.41175, 4.24005, 1.0]
    state_freq_MLEs = [0.343455, 0.385235, 0.082384, 0.188926]
    shape_MLE       = 0.357591
    pinvar_MLE      = 0.116526

    if scheme == 1:
        model.type = 'gtr'
        
        # relative rates
        model.update_relrates_separately = False
        model.relrates = relrate_MLEs
        model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
        
        # state frequencies
        model.update_freqs_separately = False
        model.state_freqs = state_freq_MLEs
        model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))
        
        # discrete gamma rate heterogeneity
        model.num_rates = 4
        model.gamma_shape = shape_MLE
        model.gamma_shape_prior = Exponential(1.0)
        
        # proportion of invariable sites
        model.pinvar_model = True
        model.pinvar = pinvar_MLE
        model.pinvar_prior = Beta(1.0, 1.0)
        
        # edge lengths
        model.edgelen_hyperprior = None    
        model.edgelen_prior = Exponential(10.0)
        
    elif scheme == 2:
        model.type = 'gtr'
        
        # relative rates
        model.update_relrates_separately = False
        model.relrates = relrate_MLEs
        model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
        
        # state frequencies
        model.update_freqs_separately = False
        model.state_freqs = state_freq_MLEs
        model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))
        
        # discrete gamma rate heterogeneity
        model.num_rates = 4
        model.gamma_shape = shape_MLE
        model.gamma_shape_prior = Exponential(1.0)
        
        # proportion of invariable sites
        model.pinvar_model = True
        model.pinvar = pinvar_MLE
        model.pinvar_prior = Beta(1.0, 1.0)
        
        # edge lengths
        model.edgelen_hyperprior = None    
        model.edgelen_prior = Exponential(10.0)
        
        m1 = model()
        m2 = model()

        partition.addSubset(mtDNA, m1, 'mtDNA')
        partition.addSubset(nuc, m2, 'nuc')
        partition()
        
    elif scheme == 3:
        abort('Sorry, not currently able to implement scheme 3')
        
    elif scheme == 4:
        abort('Sorry, not currently able to implement scheme 4')
        
    elif scheme == 5:
        abort('Sorry, not currently able to implement scheme 5')
        
    elif scheme == 6:
        abort('Sorry, not currently able to implement scheme 6')
        
    elif scheme == 7:
        abort('Sorry, not currently able to implement scheme 7')
        
    elif scheme == 8:
        abort('Sorry, not currently able to implement scheme 8')
        
    elif scheme == 9:
        abort('Sorry, not currently able to implement scheme 9')
        
    else:
        abort('Sorry, unrecognized scheme (%s)' % str(scheme))
        
setMasterSeed(13579)
data_file_name = 'M3354.nex'
blob = readFile(data_file_name)

createModelForScheme(scheme)

mcmc.data_source = blob.characters
mcmc.out.log.prefix = fnprefix
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = fnprefix
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = fnprefix
mcmc.out.params.mode = REPLACE
mcmc.nchains = 1

if just_testing:
    mcmc.ncycles = 2
    mcmc.adapt_first = 3
    mcmc.sample_every = 1
    mcmc.report_every = 1
else:
    mcmc.ncycles = 20000
    mcmc.adapt_first = 2
    mcmc.sample_every = 10
    mcmc.report_every = 100
    
mcmc.edge_move_weight = 1
mcmc.ls_move_weight = 100
mcmc.state_freq_weight = 10
mcmc.state_freq_psi = 500.0
mcmc.rel_rate_weight = 10
mcmc.rel_rate_psi = 500.0
mcmc.tree_scaler_weight = 1
mcmc.slice_weight = 1
mcmc.rel_rate_psi = 1000.0
mcmc.state_freq_psi = 1000.0
#mcmc.starting_tree_source = randomtree(n_taxa=len(blob.taxon_labels))
mcmc.starting_tree_source = TreeCollection(newick=user_tree_def)
mcmc.debugging = False
mcmc.fix_topology = fixed_topology
mcmc()

#ss.nbetavals = 11
#ss()

if not just_testing:
    sump.file = fnprefix + '.p'
    sump.out.log = fnprefix + '.sump.txt'
    sump.out.log.mode = REPLACE
    sump()


