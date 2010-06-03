import os,sys,math,random
from phycas import *
from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from MCMCManager import MCMCManager
from phycas.ProbDist import StopWatch
from phycas.ReadNexus import NexusReader


DEBUGGING_OUTPUT = False

def check(msg = 'check'):
	raw_input(msg)

class MCMCImpl(CommonFunctions):
	#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
	"""
	Needs to be written.
	
	"""
	def __init__(self, opts):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Initializes MCMCImpl object by assigning supplied phycas object
		to a data member variable.
		
		"""
		CommonFunctions.__init__(self, opts)
		self.models					= None		# This variable is set in setup
		
		# These copied over from Phycas.py - many are not used and should be weeded out
		self.data_matrix			= None
		self.file_name_trees_stored = None
		self.do_marginal_like		= False
		self.mcmc_manager			= MCMCManager(self)
		self.heat_vector			= None		# Leave set to None unless you are implementing some ad hoc heating scheme. This vector ordinarily computed using nchains and self.heating_lambda
		self.stopwatch				= StopWatch()
		self.sim_model_tree			= None		# Will hold the model tree used by simulateDNA 
		self.starting_tree			= []		# starting_tree[i] will contain reference to Phylogeny.Tree object for chain i
		self.warn_tip_numbers		= False		# True only if tip numbers were not able to be created using the tip names in the tree description (always False if starting_tree_source == 'random' because BuildTreeFromString is not called in this case)
		self.ntax					= 0			# Will hold the actual number of taxa after data file read
		self.nchar					= 0			# Will hold the actual number of characters after data file has been read
		self.npatterns			= []		# Will hold the actual number of patterns for each subset after data file has been read
		self.taxon_labels			= []		# Will hold taxon labels from data file or default names if self.data_source equals None
		#self.sssf					= None
		self.paramf					= None
		self.treef					= None
		self.sitelikef				= None
		#self.tree_file_name		 = ''		 # Will hold tree file name (see openParameterAndTreeFiles)
		#self.param_file_name		 = ''		 # Will hold parameter file name (see openParameterAndTreeFiles)
		#self.tmp_simdata			 = SimData()
		self.gg_Pm					= 0.0		# Penalty component (same for all k)
		self.gg_Gm					= []		# Vector of goodness-of-fit components (one for each k in gg_kvect)
		self.gg_Dm					= []		# Vector of overall measures (one for each k in gg_kvect)
		self.reader					= NexusReader()
		#self._logFileName			 = None
		self.addition_sequence		= []		# List of taxon numbers for addition sequence
		self.ref_tree				= None		# a reference tree against which the current MCMC tree can be compared (user specifies with, e.g., "mcmc.reference_tree_source = TreeCollection(file='hkyml.tre')")
		self.stored_tree_defs		= None
		self.psf					= None
		self.pdf_splits_to_plot		= None
		#self.param_file_name		 = None	 
		#self.tree_file_name		 = None
		self.nsamples				= None
		self.unimap_manager			= None
		self.nsamples				= 0
		self.burnin					= 0		# same as self.opts.burnin except for path sampling, when it drops to 0 after first beta value
		self.ncycles				= 0		# same as self.ncycles except for steppingstone sampling, when ss.xcycles are added to mcmc.ncycles for the first (beta = 1) step
		self.cycle_start			= 0		# used in path sampling to avoid starting over the cycle count for each beta value
		self.cycle_stop				= 0		# total number of cycles (used for computing time remaining) 
		self.last_adaptation		= 0
		self.next_adaptation		= 0
		self.ss_beta				= 1.0
		self.ss_beta_index			= 0
		self.ss_sampled_betas		= None
		self.ss_sampled_likes		= None
		self.siteIndicesForPatternIndex = None
		
	def setSiteLikeFile(self, sitelikef):
		if sitelikef is not None:
			self.sitelikef = sitelikef
			
	def siteLikeFileSetup(self, coldchain):
		if self.sitelikef is not None:
			# Set up the siteIndicesForPatternIndex, which holds a list of sites for each pattern index
			# This allows us to spit out site likelihoods for each site, even though site likelihoods
			# are stored for patterns, many of which represent numerous sites
			v = coldchain.likelihood.getCharIndexToPatternIndex()
			self.phycassert(len(v) == self.nchar,'Number of sites returned by coldchain.likelihood.getCharIndexToPatternIndex differs from MCMCImpl.nchar in MCMCImpl.siteLikeFileSetup()')
			npatterns = coldchain.likelihood.getNPatterns()
			self.siteIndicesForPatternIndex = []
			for i in range(npatterns):
				self.siteIndicesForPatternIndex.append([])
			for i,p in enumerate(v):
				self.siteIndicesForPatternIndex[p].append(i)

	def unsetSiteLikeFile(self):
		self.sitelikef = None
		self.siteIndicesForPatternIndex = None

	def adaptSliceSamplers(self):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Cycles through all slice samplers and adapts each one. Adaptation of
		a slice sampler involves changing the unit width of segments used to
		construct a slice. If the slice unit width is too small, too many
		likelihood function evaluations are needed to span the full
		conditional posterior density. If the slice unit width is too large,
		too many failed sampling attempts are needed before a valid sample
		can be obtained. Adaptation adjusts the slice unit width of each
		slice sampler in an attempt to bring it closer to the optimum width
		using experience from past sampling attempts.
		
		"""
		summary = ''
		# need to adapt all chains, not just the cold one!
		cold_chain_manager = self.mcmc_manager.getColdChainManager()
		for p in cold_chain_manager.getAllUpdaters():
			p_summary = ''
			nm = p.getName()
			if p.hasSliceSampler():
				s = p.getSliceSampler()
				if s.getNumSamples() > 0:
					s.adaptSimple(self.opts.adapt_simple_param)
					mode = s.getMode()
					accept_pct = 100.0*float(s.getNumSamples())/float(s.getNumFuncEvals())
					p_summary += ' * efficiency = %.1f%%, mode=%.5f (%s)\n' % (accept_pct, mode, nm)
					s.resetDiagnostics()
			else:
				naccepts = p.getNumAccepts()
				nattempts = p.getNumAttempts()
				if nattempts > 0:
					accept_pct = 100.0*float(naccepts)/float(nattempts)
					if nattempts == 1:
						p_summary += '   accepted %.1f%% of 1 attempt (%s)\n' % (accept_pct, nm)
					else:
						p_summary += '   accepted %.1f%% of %d attempts (%s)\n' % (accept_pct, nattempts, nm)
				else:
					accept_pct = 0.0
				p.resetDiagnostics()
			summary += p_summary
		
		if self.opts.verbose and summary != '':
			self.output('\nUpdater diagnostics (* = slice sampler):')
			self.output(summary)
			
	def obsoleteUpdateAllUpdaters(self, chain, chain_index, cycle):
		# This function abandoned; functionality moved to C++ side 
		# for speed reasons: see MCMCChainManager::updateAllUpdaters
		#import gc
		#gc.enable()
		if self.opts.debugging:
			tmpf = file('debug_info.txt', 'a')
			tmpf.write('************** cycle=%d, chain=%d\n' % (cycle,chain_index))
		for p in chain.chain_manager.getAllUpdaters():
			w = p.getWeight()
			#print "param = %s (weight = %f), chain = %d" % (p.getName(), w, chain_index)
			for x in range(w):
				if self.opts.debugging:
					p.setSaveDebugInfo(True)
				p.update()
				#if p.getName() == 'Bush move':
				#	 print '  counts after bush move =',gc.get_count()
				#	 print '  thresholds =',gc.get_threshold()
				#	 #raw_input('debug stop')
				if self.opts.debugging:
					tmpf.write('%s | %s\n' % (p.getName(), p.getDebugInfo()))
		
		if self.opts.debugging:
			tmpf.close()

	def showTopoPriorInfo(self):
		m = self.mcmc_manager.getColdChain()
		self.output('\nTopology prior:')
		if not self.opts.allow_polytomies:
			self.output('  flat across all fully-resolved tree topologies (polytomies not allowed)')
		else:			 
			if m.topo_prior_calculator.isPolytomyPrior():
				self.output('  Prior type: polytomy prior')
			else:
				self.output('  Prior type: resolution class prior')
			self.output('  Prior strength (C): %s' % m.topo_prior_calculator.getC())
			self.output('  Prior probability for each resolution class:')
			self.output('  Note: 0.00000000 does *not* mean that the prior is zero! It simply')
			self.output('		 indicates that the prior is less than 0.000000005\n')
			self.output('%20s %20s' % ('internal nodes', 'prior probability'))
			self.output('%20s %20s' % ('--------------', '-----------------'))
			topo_priors = m.topo_prior_calculator.getRealizedResClassPriorsVect()
			for i,v in enumerate(topo_priors):
				if i == 0:
					denom = v	# first element of vector is log of normalization constant (sum of all other elements)
				else:
					topo_prior = math.exp(v - denom)
					self.output('%20d %20.8f' % (i,topo_prior))
			self.output()

	def showParamInfo(self, p):
		if p.computesUnivariatePrior() or p.computesMultivariatePrior():
			self.output('  Parameter name:	   %s' % p.getName())
			self.output('  Prior distribution: %s' % p.getPriorDescr())
			if p.isMasterParameter():
				self.output('  Master parameter (no current value)')
			else:
				if p.computesUnivariatePrior():
					v = p.getCurrValueFromModel()
					self.output('  Current value:	   %s' % v)
				else:
					v = p.listCurrValuesFromModel()
					self.output('  Current value:	   %s' % ','.join(['%.5f' % x for x in v]))				
			self.output('  Prior log-density:  %s' % p.getLnPrior())
			self.output()
				
	def treeFileOpen(self):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Opens the tree file and writes a translate table.
		
		"""
		#self.tree_file_name = self.opts.out.trees
		
		tree_file_spec = self.opts.out.trees
		self.treef = None
		try:
			self.treef = tree_file_spec.open(self.stdout)
		except:
			print '*** Attempt to open tree file (%s) failed.' % self.opts.out.trees.filename

		if self.treef:
			self.mcmc_manager.treeFileHeader(self.treef)

	def treeFileClose(self):
		self.treef.write('end;\n')
		self.treef.close()

	def paramFileOpen(self):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Opens the parameter file and writes a header line. Additionally, if
		performing smoothed steppingstone sampling, opens a file with 
		extension .sss to hold the likelihood, prior and working prior samples
		needed to estimate the marginal likelihood.
		
		"""
		#self.param_file_name = self.opts.out.params
		
		param_file_spec = self.opts.out.params
		self.paramf = None
		try:
			self.paramf = param_file_spec.open(self.stdout)
		except:
			print '*** Attempt to open parameter file (%s) failed.' % self.opts.out.params.filename

		if self.paramf:
			self.mcmc_manager.paramFileHeader(self.paramf)
			self.paramf.write('\n')
			
	def paramFileClose(self):
		if self.paramf is not None:
			self.paramf.close()
		#if self.sssf is not None:
		#	self.sssf.close()

	def openParameterAndTreeFiles(self):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Creates parameter and tree file names based on the data file name or the
		user-supplied prefix and opens the files
		
		"""
		#prefix = self.getPrefix()
		#self.param_file_name = prefix + '.p'
		#self.tree_file_name = prefix + '.t'

		self.paramFileOpen()
		self.treeFileOpen()
		
	def _loadData(self, matrix):
		self.phycassert(matrix is not None, 'Tried to load data from a non-existant matrix')
		self.data_matrix = matrix
		self.taxon_labels = matrix.getTaxLabels()
		self.ntax = self.data_matrix.getNTax()
		self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only
		self.phycassert(len(self.taxon_labels) == self.ntax, "Number of taxon labels does not match number of taxa.")
		
	def cumLagrange(self, which, x, y):
		xx = x[which]
		term1 = (xx - x[0])*(y[0] + (xx - x[0])*(y[1] - y[0])/(2.0*(x[1] - x[0])))
		term2a = 2.0*(xx**2.0) - xx*x[0] - (x[0]**2.0) + 2.0*x[0]*x[1] - 3.0*xx*x[1]
		term2b = (y[2] - y[1])/(x[2] - x[1]) - (y[1] - y[0])/(x[1] - x[0])
		term2c = (xx - x[0])/(x[2] - x[0])
		term2 = term2a*term2b*term2c/6.0
		cum = term1 + term2
		return cum
		
	def getStartingTree(self):
		#		  if self.starting_tree is None:
		#			  if False:
		#				  if self.opts.starting_tree_source == 'file':
		#					  self.phycassert(self.data_source, "Specified starting_tree_source to be 'file' when data_source was None (file was not read)")
		#					  tree_defs = self.reader.getTrees()
		#					  self.phycassert(len(tree_defs) > 0, 'a trees block defining at least one tree must be stored in the nexus data file')
		#					  # Grab first tree description in the data file
		#					  # TODO allow some other tree than the first
		#					  self.starting_tree = tree_defs[0]
		#				  elif self.opts.starting_tree_source == 'usertree':
		#					  self.starting_tree = Newick(self.opts.tree_topology)
		#				  elif self.opts.starting_tree_source == 'random':
		#					  self.phycassert(self.ntax > 0, 'expecting ntax to be greater than 0')
		#					  self.starting_tree = None
		#				  else:
		#					  self.phycassert(False, "starting_tree_source should equal 'random', 'file', or 'usertree', but instead it was this: %s" % self.starting_tree_source)
		#			  else:
		# If user failed to specify starting_tree_source, get starting tree from randomtree object
		# as it is currently configured
		tr_source = self.opts.starting_tree_source
		if tr_source is None:
			tr_source = randomtree()
		try:
			tr_source.setActiveTaxonLabels(self.taxon_labels)
			i = iter(tr_source)
			# self.starting_tree = i.next()
			self.starting_tree.append(i.next())
		except:
			self.stdout.error("A starting tree could not be obtained from the starting_tree_source")
			raise
		t = self.starting_tree[-1]
		num_degree_two_nodes = t.deroot()
		if num_degree_two_nodes > 0:
			self.stdout.warning("A total of %d degree-2 nodes were removed from tree defined in starting_tree_source" % num_degree_two_nodes)
		return t
		
	def storeRefTreeIfSupplied(self):
		cold_chain = self.mcmc_manager.getColdChain()
		
		# If a reference tree was specified, create that tree now
		tr_source = self.opts.reference_tree_source
		if tr_source is not None:
			try:
				tr_source.setActiveTaxonLabels(self.taxon_labels)
				i = iter(tr_source)
				self.ref_tree = i.next()
			except:
				self.stdout.error("A reference tree could not be obtained from the specified reference_tree_source")
		if self.ref_tree is not None:
			cold_chain.chain_manager.setRefTree(self.ref_tree)

	def setup(self):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		This function is for parts of the setup that should occur right before
		run() is called. Setup is deferred until this point to give the
		user a chance to change the default settings before the call to
		run(). This function does these things: 
		1) reads the data and ensures that the taxon_labels list is filled
		with the correct number of taxon labels; 
		2) creates a starting tree description; 
		3) creates an appropriate heat_vector
		4) calls MCMCManager's createChains function to handle setup for each
		individual chain; 
		5) opens the parameter and tree files; and
		6) establishes an output log file name if requested
		
		"""
		ds = self.opts.data_source
		if ds is None:
			# User apparently wants to run without data
			self.data_matrix = None
			self.ntax = self.opts.ntax
			self.nchar = 0 # used for Gelfand-Ghosh simulations only
			self.phycassert(self.ntax > 0, 'Number of taxa (mcmc.ntax) should be > 0 if mcmc.data_source is None')
			self.taxon_labels = ['taxon%d' % (i+1,) for i in range(self.ntax)]
		else:
			mat = ds.getMatrix()
			self.phycassert(mat is not None, 'Data matrix could not be input')
			self._loadData(mat)
			self.phycassert(self.ntax > 0, 'Number of taxa in data matrix was 0')
			
		
		#print 'In MCMCImpl.py, function setup():'
		#print '  mat =', mat
		#print '  self.nchar = %d' % self.nchar
		#raw_input('debug check')
		
		# Next line creates a default partition if a partition was not defined by the user
		self.opts.partition.validate(self.nchar)

		self.models			= [m for (n,s,m) in self.opts.partition.subset]
		self.model_names	= [n for (n,s,m) in self.opts.partition.subset]
		
		for m,n in zip(self.models, self.model_names):
			#print '==> checking model %s' % n
			if m.edgelen_prior is not None:
				# set both internal and external edge length priors to edgelen_prior
				m.internal_edgelen_prior = m.edgelen_prior
				m.external_edgelen_prior = m.edgelen_prior
			else:
				# Ensure that user has specified both internal and external edge length priors
				self.phycassert(m.internal_edgelen_prior is not None, 'In model %s, internal_edgelen_prior cannot be None if edgelen_prior is None' % n)
				self.phycassert(m.external_edgelen_prior is not None, 'In model %s, external_edgelen_prior cannot be None if edgelen_prior is None' % n)
				
			if m.edgelen_hyperprior is not None:
				# Ensure that both internal and external edgelen priors are Exponential
				if m.internal_edgelen_prior.getDistName() != 'Exponential':
					m.internal_edgelen_prior = Exponential(1.0)
					self.warning('In model %s, internal_edgelen_prior reset to Exponential because edgelen_hyperprior was specified' % n)
				if m.external_edgelen_prior.getDistName() != 'Exponential':
					m.external_edgelen_prior = Exponential(1.0)
					self.warning('In model %s, external_edgelen_prior reset to Exponential because edgelen_hyperprior was specified' % n)
			
		# Determine heating levels if multiple chains
		if self.opts.heat_vector == None:
			if self.opts.nchains == 1:
				self.heat_vector = [1.0]
			else:
				# Determine vector of powers for each chain
				self.heat_vector = []
				for i in range(self.opts.nchains):
					# For n=5 chains (1 cold, 4 heated), min_heat_power = 0.5, we have:
					# lambda = (1 - min_heat_power)/(min_heat_power*(n-1))
					#		 = (1 - 0.5)/(0.5*4)
					#		 = 0.25
					# 0 1.000 = 1/1.0 cold chain explores posterior
					# 1 0.800 = 1/1.25
					# 2 0.667 = 1/1.50
					# 3 0.571 = 1/1.75
					# 4 0.500 = 1/2.00
					z = self.opts.min_heat_power
					n = self.opts.nchains
					lamda = (1.0 - z)/(z*(n - 1.0))
					temp = 1.0/(1.0 + float(i)*lamda)
					self.heat_vector.append(temp)
		else:
			# User supplied his/her own heat_vector; perform sanity checks
			self.heat_vector = self.opts.heat_vector
			#self.opts.nchains = len(self.heat_vector)
			self.phycassert(self.heat_vector[0] == 1.0, 'user-supplied heat_vector does not allow for a cold chain (one power must be 1.0)')
			h = list(self.heat_vector)
			h.sort(reverse=True)
			if h != self.heat_vector:
				self.phycassert(False, 'chain heating powers must be in decreasing order')
			self.phycassert(h[-1] > 0.0, 'all chain heating powers must be positive')

		self.mcmc_manager.createChains()
		self.storeRefTreeIfSupplied()
		self.openParameterAndTreeFiles()
		
		if self.opts.doing_steppingstone_sampling:
			# start with posterior (ss_beta = 1) and work toward the prior (ss_beta = 0)
			self.ss_beta = self.opts.ssobj.maxbeta		
			cc = self.mcmc_manager.getColdChain()
			cc.setPower(self.ss_beta)
		self.siteLikeFileSetup(self.mcmc_manager.getColdChain())
		
	def beyondBurnin(self, cycle):
		c = cycle + 1
		return (c > self.burnin)
		
	def doThisCycle(self, cycle, mod):
		c = cycle + 1
		return ((c % mod) == 0)
		
	def getModelIndex(self, name):
		"""
		If an updater has name like 'rAC_k', this function will return the model index k-1.
		If an updater has name 'rAC', returns model index 0.
		"""
		import re
		mo = re.search('_([0-9]+)$',name)
		if mo is not None:
			return int(mo.group(1)) - 1
		else:
			return 0
		
	############################ exploreWorkingPrior ############################
	def exploreWorkingPrior(self, cycle):
		#raw_input('entering exploreWorkingPrior')
		chain_index = 0
		cold_chain = self.mcmc_manager.getColdChain()
		#tm = Phylogeny.TreeManip(cold_chain.tree)
		
		nmodels = cold_chain.partition_model.getNumSubsets()
		unpartitioned = (nmodels == 1)

		new_edge_lens = []
		new_internal_edge_lens = None
		new_external_edge_lens = None
		
		all_updaters = cold_chain.chain_manager.getAllUpdaters() 
		for u in all_updaters:		# good candidate for moving into C++
			if u.isFixed() or not u.isPriorSteward():
				continue
			name = u.getName()
			if name.find('relrates') == 0:						# C++ class RelRatesMove
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				rate_vector = u.sampleMultivariateWorkingPrior()
				m.setRelRates(rate_vector)
			elif name.find('state_freqs') == 0:					# C++ class StateFreqMove
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				freq_vector = u.sampleMultivariateWorkingPrior()
				m.setStateFreqsUnnorm(freq_vector)
			elif name.find('subset_relrates') == 0:					# C++ class SubsetRelRatesMove
				rates_vector = u.sampleMultivariateWorkingPrior()
				while min(rates_vector) == 0.0:
					self.output('\nWarning: resampling working prior for subset relative rates because first draw contained\n  an element equal to zero: (%s)' % ','.join(['%g' % r for r in rates_vector]))
					rates_vector = u.sampleMultivariateWorkingPrior()
				cold_chain.partition_model.setSubsetRelRatesVect(rates_vector)
			elif name.find('gamma_shape') == 0:						# C++ class DiscreteGammaShapeParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				new_gamma_shape = u.sampleWorkingPrior()
				m.setShape(new_gamma_shape)
			elif name.find('pinvar') == 0:							# C++ class PinvarParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				new_pinvar = u.sampleWorkingPrior()
				m.setPinvar(new_pinvar)
			elif name.find('kappa') == 0:							# C++ class KappaParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				new_kappa = u.sampleWorkingPrior()
				m.setKappa(new_kappa)
			elif name.find('rAC') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				rr = list(m.getRelRates())
				new_rAC = u.sampleWorkingPrior()
				rr[0] = new_rAC
				m.setRelRates(rr)
			elif name.find('rAG') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				rr = list(m.getRelRates())
				new_rAG = u.sampleWorkingPrior()
				rr[1] = new_rAG
				m.setRelRates(rr)
			elif name.find('rAT') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				rr = list(m.getRelRates())
				new_rAT = u.sampleWorkingPrior()
				rr[2] = new_rAT
				m.setRelRates(rr)
			elif name.find('rCG') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				rr = list(m.getRelRates())
				new_rCG = u.sampleWorkingPrior()
				rr[3] = new_rCG
				m.setRelRates(rr)
			elif name.find('rCT') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				rr = list(m.getRelRates())
				new_rCT = u.sampleWorkingPrior()
				rr[4] = new_rCT
				m.setRelRates(rr)
			elif name.find('rGT') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				rr = list(m.getRelRates())
				new_rGT = u.sampleWorkingPrior()
				rr[5] = new_rGT
				m.setRelRates(rr)
			elif name.find('freqA') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				new_freqA = u.sampleWorkingPrior()
				m.setStateFreqUnnorm(0, new_freqA)
			elif name.find('freqC') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				new_freqC = u.sampleWorkingPrior()
				m.setStateFreqUnnorm(1, new_freqC)
			elif name.find('freqG') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				new_freqG = u.sampleWorkingPrior()
				m.setStateFreqUnnorm(2, new_freqG)
			elif name.find('freqT') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = cold_chain.partition_model.getModel(i)
				new_freqT = u.sampleWorkingPrior()
				m.setStateFreqUnnorm(3, new_freqT)
			elif name.find('edgelen') == 0:							# C++ class EdgeLenParam
				# depends on edge length updaters being in preorder sequence
				edgelen = u.sampleWorkingPrior()
				new_edge_lens.append(edgelen)
			elif name.find('master_edgelen') == 0:					# C++ class EdgeLenMasterParam
				num_edge_lens = cold_chain.tree.getNNodes() - 1
				new_edge_lens = [u.sampleWorkingPrior() for j in range(num_edge_lens)]
			elif name.find('external_edgelen') == 0:				# C++ class EdgeLenMasterParam
				num_edge_lens = cold_chain.tree.getNTips() - 1
				new_external_edge_lens = [u.sampleWorkingPrior() for j in range(num_edge_lens)]
			elif name.find('internal_edgelen') == 0:				# C++ class EdgeLenMasterParam
				num_edge_lens = cold_chain.tree.getNInternals()
				new_internal_edge_lens = [u.sampleWorkingPrior() for j in range(num_edge_lens)]
			else:
				self.phycassert(0, 'model uses an updater (%s) that has not yet been added to MCMCImpl.exploreWorkingPrior (workaround: specify mcmc.draw_directly_from_prior = False)' % name)
		
		if new_internal_edge_lens is not None:
			# Case of a separate master edge length parameter for internals and tips
			self.phycassert(new_external_edge_lens is not None, 'not expecting new_external_edge_lens to be None in MCMCImpl.exploreWorkingPrior') 
			i = 0	# indexes internals
			j = 0	# indexes tips
			for nd in cold_chain.tree:
				if nd.isRoot():
					continue
				elif nd.isInternal():
					nd.setEdgeLen(new_internal_edge_lens[i])
					i += 1
				elif nd.isTip():
					nd.setEdgeLen(new_external_edge_lens[j])
					j += 1
				else:
					self.phycassert(0, 'nd is neither a tip nor an internal node in MCMCImpl.exploreWorkingPrior')
		else:
			# Case of fixed topology (and hence each edge has its own length parameter) or there is a single master edge length parameter governing both tips and internals
			self.phycassert(len(new_edge_lens) == cold_chain.tree.getNNodes() - 1, 'new_edge_lens has %d elements but expecting %d in MCMCImpl.exploreWorkingPrior' % (len(new_edge_lens), cold_chain.tree.getNNodes() - 1)) 
			i = 0
			for nd in cold_chain.tree:
				if nd.isRoot():
					continue
				nd.setEdgeLen(new_edge_lens[i])
				i += 1
					
		# replace the model
		cold_chain.prepareForLikelihood()
		cold_chain.likelihood.replaceModel(cold_chain.partition_model)
		
		# recalculate the likelihood
		cold_chain_manager = self.mcmc_manager.getColdChainManager()
		cold_chain_manager.refreshLastLnLike()
		cold_chain_manager.refreshLastLnPrior()

		
		# what about polytomies?
											
		############################ end exploreWorkingPrior ############################
		
	def explorePrior(self, cycle):
		chain_index = 0
		chain = self.mcmc_manager.getColdChain()
		tm = Phylogeny.TreeManip(chain.tree)
		if self.opts.debugging:
			tmpf = file('debug_info.txt', 'a')
			tmpf.write('************** cycle=%d, chain=%d\n' % (cycle,chain_index))
		edgelens_generated = False
		
		nmodels = chain.partition_model.getNumSubsets()
		unpartitioned = (nmodels == 1)

		for p in chain.chain_manager.getAllUpdaters():
			if p.isFixed():
				continue
			w = p.getWeight()
			name = p.getName()
			if (name.find('edgelen_hyper') == 0) or (name.find('external_hyper') == 0) or (name.find('internal_hyper') == 0):	# C++ class HyperPriorParam
				if not edgelens_generated:
					# Choose hyperparam, then use it to choose new edge lengths for a newly-created tree
					m = chain.partition_model.getModel(0)
					if m.isSeparateInternalExternalEdgeLenPriors():
						# draw an edge length hyperparameter value for external edges
						edgelen_hyperparam = m.getEdgeLenHyperPrior().sample()
						chain.chain_manager.setEdgeLenHyperparam(0, edgelen_hyperparam)
						m.getExternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
						
						# draw an edge length hyperparameter value for internal edges
						edgelen_hyperparam = m.getEdgeLenHyperPrior().sample()
						chain.chain_manager.setEdgeLenHyperparam(1, edgelen_hyperparam)
						m.getInternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
					else:
						# draw an edge length hyperparameter value that applies to all edges
						edgelen_hyperparam = m.getEdgeLenHyperPrior().sample()
						chain.chain_manager.setEdgeLenHyperparam(0, edgelen_hyperparam)
						m.getInternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
						m.getExternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
					if self.opts.fix_topology:
						tm.setRandomInternalExternalEdgeLengths(m.getInternalEdgeLenPrior(), m.getExternalEdgeLenPrior()) 
					else:
						tm.equiprobTree(chain.tree.getNTips(), chain.r, m.getInternalEdgeLenPrior(), m.getExternalEdgeLenPrior())
					edgelens_generated = True
			elif name.find('master_edgelen') == 0:
				pass
			elif name.find('external_edgelen') == 0:
				pass
			elif name.find('internal_edgelen') == 0:
				pass
			elif name.find('subset_relrates') == 0:
				new_subset_relrate_vector = chain.partition_model.getSubsetRelRatePrior().sample()
				chain.partition_model.setSubsetRelRatesVect(new_subset_relrate_vector)
			elif name.find('kappa') == 0:							# C++ class KappaParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_kappa = m.getKappaPrior().sample()
				m.setKappa(new_kappa)
			elif name.find('omega') == 0:							# C++ class OmegaParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_omega = m.getOmegaPrior().sample()
				m.setOmega(new_omega)
			elif name.find('rAC') == 0:								# C++ class GTRRateParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_rAC = m.getRelRateParamPrior().sample()
				m.setRelRateUnnorm(0, new_rAC) 
			elif name.find('rAG') == 0:								# C++ class GTRRateParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_rAG = m.getRelRateParamPrior().sample()
				m.setRelRateUnnorm(1, new_rAG) 
			elif name.find('rAT') == 0:								# C++ class GTRRateParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_rAT = m.getRelRateParamPrior().sample()
				m.setRelRateUnnorm(2, new_rAT) 
			elif name.find('rCG') == 0:								# C++ class GTRRateParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_rCG = m.getRelRateParamPrior().sample()
				m.setRelRateUnnorm(3, new_rCG) 
			elif name.find('rCT') == 0:								# C++ class GTRRateParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_rCT = m.getRelRateParamPrior().sample()
				m.setRelRateUnnorm(4, new_rCT) 
			elif name.find('rGT') == 0:								# C++ class GTRRateParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_rGT = m.getRelRateParamPrior().sample()
				m.setRelRateUnnorm(5, new_rGT) 
			elif name.find('freqA') == 0:							# C++ class StateFreqParam
				new_freq_param_A = m.getStateFreqParamPrior().sample()
				m.setStateFreqParam(0, new_freq_param_A) 
			elif name.find('freqC') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_freq_param_C = m.getStateFreqParamPrior().sample()
				m.setStateFreqParam(1, new_freq_param_C) 
			elif name.find('freqG') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_freq_param_G = m.getStateFreqParamPrior().sample()
				m.setStateFreqParam(2, new_freq_param_G) 
			elif name.find('freqT') == 0:							# C++ class StateFreqParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_freq_param_T = m.getStateFreqParamPrior().sample()
				m.setStateFreqParam(3, new_freq_param_T) 
			elif name.find('gamma_shape') == 0:						# C++ class DiscreteGammaShapeParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_gamma_shape = m.getDiscreteGammaShapePrior().sample()
				m.setShape(new_gamma_shape)
			elif name.find('pinvar') == 0:							# C++ class PinvarParam
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				new_pinvar = m.getPinvarPrior().sample()
				m.setPinvar(new_pinvar)
			elif name.find('relrates') == 0:						# C++ class RelRatesMove
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				rate_vector = m.getRelRatePrior().sample()
				# Drawing values from a Dirichlet prior, but relative rates should have mean 1, not sum 1,
				# so multiply each by the value nmodels to correct this. 
				m.setRelRates([6.0*x for x in rate_vector])
			elif name.find('state_freqs') == 0:						# C++ class StateFreqMove
				i = unpartitioned and 0 or self.getModelIndex(name)
				m = chain.partition_model.getModel(i)
				freq_vector = m.getStateFreqPrior().sample()
				m.setStateFreqsUnnorm(freq_vector)
			elif name.find('edge_move') == 0:						# C++ class EdgeMove
				pass
			elif name.find('larget_simon_local') == 0:				# C++ class LargetSimonMove
				pass
			elif name.find('tree_scaler') == 0:						# C++ class TreeScalerMove
				pass
			elif name.find('bush_move') == 0:						# C++ class BushMove
				pass	# polytomies handled further down (by randomly pruning fully-resolved equiprobable tree)
			else:
				self.phycassert(0, 'model uses an updater (%s) that has not yet been added to MCMCImpl.explorePrior (workaround: specify mcmc.draw_directly_from_prior = False)' % name)

		# If no edge length hyperprior was specified, then build the tree with edge lengths now

		if not edgelens_generated:
			m = chain.partition_model.getModel(0)
			if self.opts.fix_topology:
				tm.setRandomInternalExternalEdgeLengths(m.getInternalEdgeLenPrior(), m.getExternalEdgeLenPrior()) 
			else:
				tm.equiprobTree(chain.tree.getNTips(), chain.r, m.getInternalEdgeLenPrior(), m.getExternalEdgeLenPrior())
			
		if self.opts.allow_polytomies:
			# Choose number of internal nodes
			num_internal_nodes = chain.topo_prior_calculator.sample(chain.r)
					
			# Delete nodes at random from tree to achieve chosen number of internal nodes
			orig_num_internal_nodes = chain.tree.getNInternals()
			num_internals_to_delete = orig_num_internal_nodes - num_internal_nodes
			for i in range(num_internals_to_delete):
				#open('doofus.tre','w').write('%s;\n' % chain.tree.makeNumberedNewick())
				tm.deleteRandomInternalEdge(chain.r)
				
		chain.prepareForLikelihood()
		chain.likelihood.replaceModel(chain.partition_model)
				
		if False:
			# debugging code
			chain.likelihood.storeSiteLikelihoods(True)
			from phycas.Utilities.kappa2tratio import convert
			f = chain.model.getStateFreqs()
			k = convert(chain.model.getKappa(), f[0], f[1], f[2], f[3])
			print 'cycle = %d, model = %s' % (cycle + 1, chain.model.getModelName())
			print '	 lset tratio=%.5f basefreq=(%.5f %.5f %.5f) rates=gamma ncat=4 shape=%.5f;' % (k, f[0], f[1], f[2], chain.model.getShape())
			print 'taxon names:', self.opts.data_source.taxon_labels
			chain.tree.rectifyNames(self.opts.data_source.taxon_labels)
			print '	 utree one = %s;' % chain.tree.makeNewick()
			print '	 sum of edge lengths = %.5f' % chain.tree.edgeLenSum()
			raw_input('stopped before computing likelihood')
		
		# recalculate the likelihood
		cold_chain_manager = self.mcmc_manager.getColdChainManager()
		cold_chain_manager.refreshLastLnLike()
		
		if False:
			# debugging code
			counts = chain.likelihood.getPatternCounts()
			sitelikes = chain.likelihood.getSiteLikelihoods()
			print '	 lnL = %.6f' % cold_chain_manager.getLastLnLike()
			sumlikes = 0.0
			for sitelike,count in zip(sitelikes, counts):
				if count > 100:
					print '%6.0f -> %12.5f' % (count, sitelike)
				sumlikes += count*sitelike
			raw_input('check: sum = %.5f' % sumlikes)
		
		if self.opts.debugging:
			tmpf.close()
			
	def computeTimeRemaining(self, secs, ndone, ntotal):
		if ndone < 1:
			return ''
		days_remaining = 0
		hours_remaining = 0
		secs_remaining = float(secs)*(float(ntotal)/float(ndone) - 1.0)
		time_left = []
		if secs_remaining > 86400.0:
			days_remaining = math.floor(secs_remaining/86400.0)
			secs_remaining -= 86400.0*days_remaining
			if days_remaining > 0:
				if days_remaining == 1:
					time_left.append('1 day')
				else:
					time_left.append('%d days' % days_remaining)
		if secs_remaining > 3600.0:
			hours_remaining = math.floor(secs_remaining/3600.0)
			secs_remaining -= 3600.0*hours_remaining
			if hours_remaining > 0:
				if hours_remaining == 1:
					time_left.append('1 hour')
				else:
					time_left.append('%d hours' % hours_remaining)
		if secs_remaining > 60.0:
			minutes_remaining = math.floor(secs_remaining/60.0)
			secs_remaining -= 60.0*minutes_remaining
			if minutes_remaining > 0:
				if minutes_remaining == 1 and (days_remaining + hours_remaining) == 0:
					time_left.append('less than 2 minutes')
				else:
					time_left.append('%d minutes' % minutes_remaining)
		if len(time_left) > 0:
			return ', '.join(time_left) + ' remaining'
		elif math.floor(secs_remaining) == 1:
			return '1 second remaining'
		else:
			return '%d seconds remaining' % math.floor(secs_remaining)
		
	def mainMCMCLoop(self, explore_prior = False):
		levels_file_created = False	#temp!
		nchains = len(self.mcmc_manager.chains)
		# print '******** nchains =',nchains
		self.last_adaptation = 0
		self.next_adaptation = self.opts.adapt_first
		
		for cycle in xrange(self.burnin + self.ncycles):
			# Update all updaters
			if explore_prior and self.opts.draw_directly_from_prior:
				if self.opts.doing_steppingstone_sampling and self.opts.ssobj.scubed:
					self.exploreWorkingPrior(cycle)
				else:
					self.explorePrior(cycle)
			else:
				for i,c in enumerate(self.mcmc_manager.chains):
					if DEBUGGING_OUTPUT:
						if cycle == 0:
							c.r.setSeed(364665646)
						print "seed=", c.r.getSeed()
						print "	  model = " + c.model.paramReport(self.mcmc_manager.parent.opts.ndecimals)
						print '	  tree rep.%d = %s;' % (cycle + 1, c.tree.makeNewick(self.mcmc_manager.parent.opts.ndecimals))
			
					c.chain_manager.updateAllUpdaters()
					
			# print '******** time_stamp =',self.mcmc_manager.getColdChain().model.getTimeStamp()
					
			# Attempt to swap two random chains
			if nchains > 1:
				self.mcmc_manager.attemptChainSwap(cycle)
	
			# Provide progress report to user if it is time
			if self.opts.verbose and self.doThisCycle(cycle, self.opts.report_every):
				# Refresh log-likelihood of cold chain if necessary
				if self.ss_beta == 0.0:
					self.mcmc_manager.getColdChainManager().refreshLastLnLike()
						
				self.stopwatch.normalize()
				secs = self.stopwatch.elapsedSeconds()
				time_remaining = self.computeTimeRemaining(secs, self.cycle_start + cycle + 1, self.cycle_stop)
				if time_remaining != '':
					time_remaining = '(' + time_remaining + ')'
				if self.opts.doing_steppingstone_sampling:
					cold_chain_manager = self.mcmc_manager.getColdChainManager()
					msg = 'beta = %.5f, cycle = %d, lnL = %.5f %s' % (self.ss_beta, cycle + 1, cold_chain_manager.getLastLnLike(), time_remaining)
				else:
					if nchains == 1:
						cold_chain_manager = self.mcmc_manager.getColdChainManager()
						msg = 'cycle = %d, lnL = %.5f %s' % (cycle + 1, cold_chain_manager.getLastLnLike(), time_remaining)
					else:
						msg = 'cycle = %d, ' % (cycle + 1)
						for k in range(nchains):
							c = self.mcmc_manager.chains[k]
							msg += 'lnL(%.3f) = %.5f, ' % (c.heating_power, c.chain_manager.getLastLnLike())
						msg += '%s' % time_remaining
				self.output(msg)

			# Sample chain if it is time
			if self.beyondBurnin(cycle) and self.doThisCycle(cycle - self.burnin, self.opts.sample_every):
				# Refresh log-likelihood(s) if necessary
				if self.ss_beta == 0.0:
					for i,c in enumerate(self.mcmc_manager.chains):
						# is this necessary?
						c.chain_manager.refreshLastLnLike()
						c.chain_manager.refreshLastLnPrior()
						
				if self.opts.doing_steppingstone_sampling and self.opts.ssobj.scubed and self.ss_beta_index == 0:
					self.mcmc_manager.recordSample(True, self.cycle_start + cycle)	# dofit = True (i.e. educate the working prior if doing SS and currently exploring the posterior)
				else:
					self.mcmc_manager.recordSample(False, self.cycle_start + cycle)	# dofit = False
				cold_chain_manager = self.mcmc_manager.getColdChainManager()
				sampled_lnL = cold_chain_manager.getLastLnLike()
				self.ss_sampled_likes[self.ss_beta_index].append(sampled_lnL)
				self.stopwatch.normalize()

			# Adapt slice samplers if it is time
			if self.doThisCycle(cycle, self.next_adaptation):
				self.adaptSliceSamplers()
				self.next_adaptation += 2*(self.next_adaptation - self.last_adaptation)
				self.last_adaptation = cycle + 1
		self.cycle_start += self.burnin + self.ncycles
		
	def run(self):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Performs the MCMC analysis. 
		
		"""		   
		self.setup()
		
		# If user has set quiet to True, then phycas.output calls will have no effect
		self.quiet = self.opts.quiet
		
		# Tell TreeLikelihood object if user wants to run with no data
		#if not self.data_source:
		#	 self.likelihood.setNoData()

		self.nsamples = self.opts.ncycles//self.opts.sample_every
		nchains = len(self.mcmc_manager.chains)
		
		cold_chain = self.mcmc_manager.getColdChain()
		
		if self.opts.verbose:
			if self.data_matrix == None:
				self.output('Data source:	 None (running MCMC with no data to explore prior)')
			else:
				self.output('Data source:	 %s' % str_value_for_user(self.opts.data_source))
				self.output('  No. taxa:                %d' % self.ntax)
				self.output('  No. included characters: %d' % cold_chain.likelihood.sumPatternCounts())
				all_missing = cold_chain.likelihood.getListOfAllMissingSites()
				num_excluded = len(all_missing)
				if num_excluded > 0:
					self.output('  *** Note: the following %d sites were automatically excluded because' % num_excluded)
					self.output('  *** they exhibited completely missing data for all taxa:')
					while len(all_missing) > 0:
						tmp = all_missing[:10]
						all_missing = all_missing[10:]
						self.output('  ***   '+','.join([str(i+1) for i in tmp]))
						
			if nchains > 1:
				for c in range(nchains):
					self.output('Starting tree for chain %d:  %s' % (c, self.starting_tree[c]))
			else:
				self.output('Starting tree:	 %s' % str(self.starting_tree[0]))
			
			if self.opts.fix_topology:
				self.output('\nTree topology fixed.\n')
			else:
				self.output('\nTree topology free to vary.\n')
			
			if self.opts.doing_steppingstone_sampling:
				self.output('\nPerforming steppingstone sampling to estimate marginal likelihood.')
				self.output('Likelihood will be raised to the power beta, and beta will be')
				self.output('decremented from 1.0 to 0.0 in a series of steps.')
				self.output('  No. steps:				%s' % self.opts.ssobj.nbetavals)
				self.output('  No. cycles per step:		%s' % self.opts.ncycles)
				self.output('  Sample every:			%s' % self.opts.sample_every)
				self.output('  No. samples per step:	%s' % self.nsamples)
				self.output('\n')
			else:
				self.output('No. cycles:	 %s' % self.opts.ncycles)
				self.output('Sample every:	 %s' % self.opts.sample_every)
				self.output('No. samples:	 %s' % self.nsamples)
			self.output('Sampled trees will be saved in %s' % str_value_for_user(self.opts.out.trees))
			self.output('Sampled parameters will be saved in %s' % str_value_for_user(self.opts.out.params))
			if self.opts.use_unimap:
				self.output('Using uniformized mapping MCMC')
			else:
				self.output('Using standard MCMC (i.e. no uniformized mapping)')

			if not self.warn_tip_numbers:
				self.output('Tip node numbers were set using the names in the tree description')
			else:
				self.output('Warning: tip node numbers were NOT set using the names in the tree description')

		if nchains == 1:
			self.output('Creating one chain (i.e. not using heated chains to improve mixing)')
		else:
			self.output('Creating %d chains with these temperatures:' % (nchains))
			for t in self.heat_vector:
				self.output('  %.5f %s' % (t, t == 1.0 and '(cold chain)' or ''))
			
		# Compute the current log-likelihood and log-prior in case first updater 
		# is a move and will thus depend on these quantities being accurate
		for c in self.mcmc_manager.chains:
			c.chain_manager.refreshLastLnLike()
			c.chain_manager.refreshLastLnPrior()
			if c.heating_power == 1.0:
				self.output('Starting log-likelihood = %s' % c.chain_manager.getLastLnLike())
				self.output('Starting log-prior = %s' % c.chain_manager.getLastLnPrior())
		
		# Show starting parameter info 
		self.output('\nParameter starting values and prior densities:')
		cold_chain_manager = self.mcmc_manager.getColdChainManager()
		for p in cold_chain_manager.getAllUpdaters():
			self.showParamInfo(p)
		#for p in cold_chain_manager.getEdgeLenParams():
		#	self.showParamInfo(p)
		#for p in cold_chain_manager.getEdgeLenHyperparams():
		#	self.showParamInfo(p)
		#for p in cold_chain_manager.getModelParams():
		#	self.showParamInfo(p)
			
		# Show updater names
		self.output('\nHere is a list of all updaters that will be used for this analysis:')
		for p in cold_chain_manager.getAllUpdaters():
			if p.getWeight() > 0:
				p.setUseWorkingPrior(False)
				if p.isMove():
					self.output('  %s (Metropolis-Hastings)' % p.getName())
				elif p.hasSliceSampler():
					self.output('  %s (slice sampler)' % p.getName())
				else:
					self.output('  %s (computes prior but does not update)' % p.getName())

		# Debugging: show data patterns
		if self.opts.debugging:
			#cold_chain = self.mcmc_manager.getColdChain()
			s = cold_chain.likelihood.listPatterns()
			print '\nDebug Info: List of data patterns and their frequencies (see TreeLikelihood::listPatterns):'
			print s

		# Show information about topology prior to be used
		self.showTopoPriorInfo()
		
		self.stopwatch.start()
		self.mcmc_manager.resetNumLikelihoodEvals()
		
		if self.opts.doing_steppingstone_sampling:
			self.output('\nSampling (%d cycles for each of the %d values of beta)...' % (self.opts.ncycles, self.opts.ssobj.nbetavals))
		else:
			self.output('\nSampling (%d cycles)...' % self.opts.ncycles)
		if self.opts.verbose:
			print
			
		# POL moved these lines to beginning of mainMCMCLoop so that adaptation cycle
		# starts again each time path sampling beta value is changed
		#self.last_adaptation = 0
		#self.next_adaptation = self.opts.adapt_first

		# Lay down first line in params file (recorded as cycle 0) containing starting values of parameters
		self.mcmc_manager.recordSample(False)
		
		if self.opts.doing_steppingstone_sampling:
			self.phycassert(self.data_matrix is not None, 'path sampling requires data')
			self.phycassert(nchains == 1, 'path sampling requires nchains to be 1')
			chain = self.mcmc_manager.getColdChain()
			if self.opts.ssobj.nbetavals > 1:
				# Set up the list ss_sampled_betas
				# Beta distribution will be divided into ss_nbetavals intervals, each of which has an equal area
				segment_area = 1.0/float(self.opts.ssobj.nbetavals - 1)
				cum_area = 0.0
				lower_boundary = 0.0
				self.ss_sampled_betas = [self.opts.ssobj.minbeta]
				total_extent = float(self.opts.ssobj.maxbeta - self.opts.ssobj.minbeta)
				betadist = ProbDist.Beta(self.opts.ssobj.shape1, self.opts.ssobj.shape2)
				for i in range(self.opts.ssobj.nbetavals - 1):
					cum_area += segment_area
					upper_boundary = betadist.getQuantile(cum_area)
					scaled_upper_boundary = self.opts.ssobj.minbeta + total_extent*upper_boundary
					self.ss_sampled_betas.append(scaled_upper_boundary)
					lower_boundary = upper_boundary
					
				# Reverse ss_sampled_betas so that sampled beta values start at 1.0 and decrease toward 0.0
				self.ss_sampled_betas.reverse()
				
				# Output the beta values that will be used
				self.output('%d %s chosen from a discrete\nBeta(%.5f, %.5f) distribution:' % (self.opts.ssobj.nbetavals, (self.opts.ssobj.nbetavals == 1 and 'value was' or 'values were'), self.opts.ssobj.shape1, self.opts.ssobj.shape2))
				for i,x in enumerate(self.ss_sampled_betas):
					self.output('%6d %12.5f' % (i+1,x))
				self.output('An MCMC analysis will be performed exploring each of the')
				self.output('power posteriors defined by these values.')
				self.output()
			else:
				self.ss_sampled_betas = [self.opts.ssobj.minbeta]
			
			# Run the main MCMC loop for each beta value in ss_sampled_betas
			self.ss_sampled_likes = []
			working_priors_calculated = False
			for self.ss_beta_index, self.ss_beta in enumerate(self.ss_sampled_betas):
				if self.ss_beta_index > 0 and self.opts.ssobj.scubed and not working_priors_calculated:
					# if using working prior with steppingstone sampling, it is now time to 
					# parameterize the working prior for all updaters so that this working prior
					# can be used in the sequel
					self.output('\nWorking prior details:')
					all_updaters = cold_chain.chain_manager.getAllUpdaters() 
					for u in all_updaters:		# good candidate for moving into C++
						if not u.isFixed():
							u.setUseWorkingPrior(True)
							if u.computesUnivariatePrior() or u.computesMultivariatePrior():
								self.output('  Finalizing working prior for %s...' % u.getName())
								u.finalizeWorkingPrior()
								self.output('  %s --> %s' % (u.getName(), u.getWorkingPriorDescr()))
					self.output()
					working_priors_calculated = True
					
				self.ss_sampled_likes.append([])
				chain.setPower(self.ss_beta)
				boldness = 100.0*(1.0 - self.ss_beta)
				chain.setBoldness(boldness)
				print 'Setting chain boldness to %g based on beta = %g' % (boldness,self.ss_beta)
				self.cycle_stop = self.opts.burnin + len(self.ss_sampled_betas)*self.opts.ncycles + self.opts.ssobj.xcycles
				if self.ss_beta_index > 0:
					self.burnin = 0
					self.ncycles = self.opts.ncycles
				else:
					self.burnin = self.opts.burnin
					self.ncycles = self.opts.ncycles + self.opts.ssobj.xcycles
					self.cycle_start = 0
				if self.ss_beta == 0.0:
					self.mainMCMCLoop(explore_prior=True)
				else:
					self.mainMCMCLoop()
		else:	# not doing steppingstone sampling
			self.ss_sampled_likes = []
			self.ss_sampled_likes.append([])
			self.ss_beta_index = 0
			self.cycle_start = 0
			self.cycle_stop = self.opts.burnin + self.opts.ncycles
			self.burnin = self.opts.burnin
			self.ncycles = self.opts.ncycles
			if self.data_matrix is None:
				self.mainMCMCLoop(explore_prior=True)
			else:
				self.mainMCMCLoop()

		self.adaptSliceSamplers()
		total_evals = self.mcmc_manager.getTotalEvals() #self.likelihood.getNumLikelihoodEvals()
		total_secs = self.stopwatch.elapsedSeconds()
		self.output('%d likelihood evaluations in %.5f seconds' % (total_evals, total_secs))
		if (total_secs > 0.0):
			self.output('  = %.5f likelihood evaluations/sec' % (total_evals/total_secs))

		if self.treef:
			self.treeFileClose()
		if self.paramf:
			self.paramFileClose()
			
