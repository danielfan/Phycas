from phycas import *

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
filename = 'M3354.nex'
blob = readFile(filename)

scheme = 1
createModelForScheme(scheme)

mcmc.data_source = blob.characters
mcmc.out.log = 'hbirds%d.txt' % scheme
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = 'hbirds%d' % scheme
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = 'hbirds%d' % scheme
mcmc.out.params.mode = REPLACE
mcmc.nchains = 1
mcmc.ncycles = 20000
mcmc.adapt_first = 2
mcmc.sample_every = 10
mcmc.report_every = 100
mcmc.ls_move_weight = 100
mcmc.tree_scaler_weight = 1
mcmc.slice_weight = 1
mcmc.rel_rate_psi = 1000.0
mcmc.state_freq_psi = 1000.0
#mcmc.starting_tree_source = randomtree(n_taxa=len(blob.taxon_labels))
mcmc.starting_tree_source = TreeCollection(newick='(1:0.04768660,(((2:0.10098358,(((74:0.06589900,80:0.06184227):0.02528591,(((76:0.04215109,(77:0.02151332,78:0.04143849):0.02269567):0.00233907,79:0.04919162):0.03324410,81:0.08278496):0.00470941):0.01163869,75:0.08327564):0.00612124):0.00813234,(((((((((3:0.13940404,89:0.12526779):0.01522875,(((9:0.13475320,(127:0.11840046,(160:0.13339025,(161:0.02006331,162:0.02142480):0.09086069):0.00941754):0.08392620):0.03956939,71:0.20635832):0.12396491,(((11:0.33283895,(23:0.40750144,50:0.20732170):0.03803436):0.02550846,109:0.69105665,129:0.34958372):0.03969887,163:0.34069059):0.06542661):0.37240874):0.00537249,((((4:0.04899118,5:0.04759248):0.02494154,(87:0.03583240,88:0.03849916):0.02854235):0.04399677,(((((25:0.04194178,35:0.03524009):0.00932649,(29:0.03543673,(31:0.00162983,32:0.00324578):0.04227134):0.00432925):0.01499588,26:0.07493367):0.00897908,(((27:0.01748666,39:0.01537163):0.01223755,(30:0.03565502,36:0.02582077):0.00443332):0.02561120,40:0.04374404):0.01453110):0.00542783,(((28:0.02202643,33:0.01832647):0.01421800,(34:0.00000269,37:0.00325012):0.04163762):0.03866394,38:0.05760168):0.00459669):0.03332763):0.06317082,(90:0.12844150,91:0.11445175):0.05903818):0.02079016):0.02745374,(((17:0.11692930,(110:0.05790205,(111:0.04485819,112:0.03573642):0.00656574):0.06608928):0.00785209,(104:0.04705200,105:0.04309560):0.05436958):0.03719303,(((21:0.05125532,22:0.06407948):0.03547400,(((94:0.00819456,95:0.01050483):0.02277788,(143:0.03767667,144:0.04086343):0.00605241):0.02413990,119:0.06207209):0.05283197):0.04772769,((72:0.15155109,84:0.13418848):0.01064449,145:0.13678191):0.03043255):0.00670630):0.00895404):0.03192576,((((((((6:0.01048439,7:0.02139101):0.03055663,(123:0.03587248,124:0.02931768):0.01093434):0.01268665,((52:0.00368660,(102:0.00405488,103:0.00659333):0.00431879):0.02458680,(92:0.01532053,93:0.01499198):0.01882948):0.01663087):0.04370688,((8:0.05521258,142:0.05405667):0.02741971,((((((61:0.01630204,68:0.03150828):0.00822851,(118:0.01759894,147:0.01385940):0.00838555):0.00546748,151:0.02945145):0.01624864,((((67:0.02423844,154:0.02836473):0.00302549,152:0.02582856):0.02167713,(69:0.02883840,(70:0.03197003,108:0.02036098):0.00721819):0.02300618):0.00173692,(153:0.02266568,155:0.01526556):0.02278339):0.00487221):0.01085886,(148:0.04892461,150:0.05093126):0.00678132):0.01168358,(146:0.01659856,149:0.01400996):0.03917163):0.02516607):0.02174773):0.01424202,((44:0.05220464,66:0.05647598):0.02156968,((130:0.04528318,132:0.04013132):0.00638956,131:0.04697534,133:0.06525193):0.01997858):0.02140898):0.00683174,((120:0.01978299,122:0.01965305):0.01634931,121:0.02416101):0.08842266):0.01114488,(((((((10:0.01061174,14:0.00611948):0.00837220,15:0.00280368):0.00182956,16:0.01838305):0.00538963,(134:0.01851026,135:0.02509363):0.00825915):0.02075760,(140:0.00836236,141:0.01403857):0.02580512):0.01100333,((18:0.01601416,51:0.01620088):0.00159429,(128:0.01247464,(136:0.01377010,137:0.01503784):0.00460270):0.00271367):0.02870943):0.05317429,(((42:0.06223725,73:0.09681727):0.01125489,96:0.08384717):0.01135838,(62:0.08759078,(63:0.00036989,64:0.00000256):0.08285430):0.00783253):0.00846028):0.01666963):0.00477906,41:0.13954814):0.01930400):0.00752049,((((12:0.02846935,13:0.03804053):0.08963056,((57:0.02984950,58:0.02747685):0.04714500,(106:0.02781025,107:0.03466439):0.05663632):0.05578246):0.00529355,(((((19:0.04420769,(59:0.02004292,60:0.03507684):0.01367395):0.02432107,48:0.04219743):0.00423469,(((43:0.03721925,47:0.04560292,125:0.04117641):0.00745851,126:0.04224107):0.00631124,(((53:0.02809543,56:0.02671735):0.01505298,54:0.03294409):0.00527076,55:0.03742727):0.00974751):0.01324152,(45:0.01031152,46:0.01145121):0.05496183):0.01464784,((156:0.02852333,157:0.02801544):0.04028629,164:0.06542955):0.00811920):0.01276752,24:0.13648654):0.00867270):0.00352657,(82:0.04955706,83:0.03279034):0.07486111):0.01277608):0.01602549,((85:0.01415294,86:0.01073645):0.07496437,(((97:0.02456780,98:0.02941261):0.03441133,100:0.05388739):0.01333286,99:0.05544148):0.02239536):0.02219923):0.00489347,(158:0.01123808,159:0.01606124):0.09140013):0.00190004,(65:0.09268878,((113:0.04475343,117:0.04063688):0.03942415,((114:0.04432701,116:0.04328489):0.01251020,115:0.04662176):0.02818692):0.03152127):0.00581219):0.00442037,(20:0.10046226,101:0.08441685):0.01417396):0.01400496,(138:0.02831079,139:0.02252448):0.04935209):0.03010906,49:0.05599866);')
mcmc.debugging = False
mcmc()



