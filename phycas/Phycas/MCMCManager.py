import os,sys,math
from phycas import *
import phycas.Phylogeny as Phylogeny
import phycas.ProbDist as ProbDist
import phycas.Likelihood as Likelihood
from LikelihoodCore import LikelihoodCore
from MarkovChain import MarkovChain

class MCMCManager:
	#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
	"""
	This class manages one to many Markov chains participating in a
	Bayesian MCMC analysis. It creates the MarkovChain objects and manages
	swapping and sampling the chains. If likelihood heating is employed,
	it also has the ability to estimate (albeit crudely) the marginal
	likelihood of the model. A single instance of MCMCManager is created
	by parent in the parent contructor.
	
	"""
	def __init__(self, parent):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Stores the parent object passed into the constructor and creates an
		empty self.chains list.
		
		"""
		self.parent = parent	# parent is MCMCImpl object
		self.chains = []
		self.swap_table = None

	def paramFileHeader(self, paramf):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Simply passes the file object paramf on to the paramFileHeader method
		of the MarkovChain object representing the cold chain.
		
		"""
		m = self.getColdChain()
		m.paramFileHeader(paramf)

	def treeFileHeader(self, treef):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Simply passes the file object treef on to the treeFileHeader method
		of the MarkovChain object representing the cold chain.
				
		"""
		m = self.getColdChain()
		m.treeFileHeader(treef)

	def createChains(self):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Creates a separate MarkovChain for every element in parent.heat_vector
		and adds it to self.chains. This function is also responsible for 
		performing several sanity checks (last chance to abort before creating
		the MCMC machinery).
		
		"""
		# Note: self.parent is the MCMCImpl object
		
		# Sanity checks
		unimap_and_flex			   = (self.parent.opts.use_unimap and self.parent.opts.model.use_flex_model)
		unimap_and_ratehet		   = (self.parent.opts.use_unimap and self.parent.opts.model.num_rates > 1)
		unimap_and_polytomies	   = (self.parent.opts.use_unimap and self.parent.opts.allow_polytomies)
		unimap_and_multiple_chains = (self.parent.opts.use_unimap and self.parent.opts.nchains > 1)
		#unimap_and_samc			= (self.parent.opts.use_unimap and self.parent.opts.doing_samc)
		#self.parent.phycassert(not unimap_and_samc, 'SAMC cannot (yet) be used in conjunction with use_unimap')
		self.parent.phycassert(not unimap_and_polytomies, 'Allowing polytomies cannot (yet) be used in conjunction with use_unimap')
		self.parent.phycassert(not unimap_and_flex, 'Flex model cannot (yet) be used in conjunction with use_unimap')
		self.parent.phycassert(not unimap_and_ratehet, 'Rate heterogeneity cannot (yet) be used in conjunction with use_unimap')
	
		# Create the chains
		# print '@@@@@ heat vector =',self.parent.heat_vector
		for i, heating_power in enumerate(self.parent.heat_vector):
			markov_chain = MarkovChain(self.parent, heating_power)
			self.chains.append(markov_chain)
			
		n = len(self.parent.heat_vector)
		self.swap_table = [[0]*n for i in range(n)]

	def getNumChains(self):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Returns the current number of chains being managed by this object
		(i.e. returns the length of the self.chains list).
		
		"""
		return len(self.chains)

	def getColdChain(self):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Returns the cold chain MarkovChain object. The cold chain is the
		first chain in the data member chains whose power equals 1.0.
		
		"""
		# Note: self.parent is the MCMCImpl object
		i = self.parent.heat_vector.index(1.0)
		return self.chains[0]

	def getColdChainManager(self):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Returns the chain manager for the cold chain. The cold chain is the
		first chain in the data member chains whose power equals 1.0.
		
		"""
		return self.chains[0].chain_manager

	def setChainPower(self, chain_index, power):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Loops through all updaters in the chain having the specified 
		chain_index, setting the power for each updater to the supplied power.
		
		"""
		# Note: self.parent is the MCMCImpl object
		self.parent.phycassert(len(self.chains) > chain_index, 'chain index specified (%d) too large for number of chains (%d)' % (chain_index, len(self.chains)))
		self.parent.phycassert(len(self.parent.heat_vector) == len(self.chains), 'length of heat vector (%d) not equal to number of chains (%d)' % (len(self.parent.heat_vector), len(self.chains)))
		self.chains[chain_index].setPower(power)

	def resetNumLikelihoodEvals(self):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Calls resetNumLikelihoodEvals method for every MarkovChain in the
		self.chains list.
		
		"""
		for c in self.chains:
			c.resetNumLikelihoodEvals()

	def setRandomSeedAllChains(self, rnseed):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Sets the master random number seed to rnseed, then passes rnseed to
		the setSeed method of every MarkovChain in the self.chains list.
		
		"""
		# Note: self.parent is the MCMCImpl object
		self.parent.opts.random_seed = rnseed
		for c in self.chains:
			c.r.setSeed(int(rnseed))

	def getTotalEvals(self):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Returns the total number of likelihood evaluations over all 
		MarkovChain objects in the self.chains list by calling the 
		getNumLikelihoodEvals method for each chain.
		
		"""
		total = 0
		for c in self.chains:
			total += c.getNumLikelihoodEvals()
		return total

	def recordSample(self, cycle = -1):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Records the current tree topology and edge lengths by adding a line to
		the tree file, and records tree length and substitution parameters
		by adding a line to the parameter file.
		
		"""
		# Note: self.parent is the MCMCImpl object
		float_format_str = '%%.%df\t' % self.parent.opts.ndecimals
		float_format_notab_str = '%%.%df' % self.parent.opts.ndecimals
		
		# Gather log-likelihoods, and if path sampling save in path_sample list for later
		lnLikes = []
		for i,c in enumerate(self.chains):
			lnLi = c.chain_manager.getLastLnLike()
			lnLikes.append(lnLi)
		
		# Only record samples from the current cold chain
		cold_chain = self.parent.mcmc_manager.getColdChain()
		
		# Add line to parameter file if it exists
		if self.parent.paramf:
			self.parent.paramf.write('%d\t' % (cycle + 1))
			if self.parent.opts.doing_steppingstone_sampling:
				self.parent.paramf.write(float_format_str % (cold_chain.heating_power))
				
			# lnL for each chain
			for lnl in lnLikes:
				self.parent.paramf.write(float_format_str % lnl)
			
			# tree length
			if self.parent.opts.fix_topology:
				self.parent.paramf.write('%s\t' % '\t'.join([float_format_notab_str % brlen for brlen in cold_chain.tree.edgeLens()]))
			else:
				self.parent.paramf.write(float_format_str % cold_chain.tree.edgeLenSum())
			
			if partitioning:
				# record parameter values for each model in partition
				nmodels = cold_chain.partition_model.getNumSubsets()
				if nmodels > 1:
					for i in range(nmodels):
						ssrr = cold_chain.partition_model.getSubsetRelRate(i)
						self.parent.paramf.write(float_format_str % ssrr)					
				for i in range(nmodels):
					m = cold_chain.partition_model.getModel(i)
					self.parent.paramf.write(m.paramReport(self.parent.opts.ndecimals))
					if m.hasEdgeLenHyperPrior():
						if m.isSeparateInternalExternalEdgeLenPriors():
							self.parent.paramf.write(float_format_str % m.getExternalEdgelenHyperparam())
							self.parent.paramf.write(float_format_str % m.getInternalEdgelenHyperparam())
						else:
							self.parent.paramf.write(float_format_str % m.getInternalEdgelenHyperparam())
					if m.isFlexModel():
						rates_vector = cold_chain.likelihood.getRateMeans()
						for rr in rates_vector:
							self.parent.paramf.write(float_format_str % rr)
						probs_vector = cold_chain.likelihood.getRateProbs()
						for rp in probs_vector:
							self.parent.paramf.write(float_format_str % rp)
				self.parent.paramf.write('\n')
				self.parent.paramf.flush()
			else:
				self.parent.paramf.write(cold_chain.model.paramReport(self.parent.opts.ndecimals))
				if self.parent.opts.model.edgelen_hyperprior is not None:
					for p in cold_chain.chain_manager.getEdgeLenHyperparams():
						self.parent.paramf.write(float_format_str % p.getCurrValue())
				if self.parent.opts.model.use_flex_model:
					rates_vector = cold_chain.likelihood.getRateMeans()
					for rr in rates_vector:
						self.parent.paramf.write(float_format_str % rr)
					probs_vector = cold_chain.likelihood.getRateProbs()
					for rp in probs_vector:
						self.parent.paramf.write(float_format_str % rp)
				self.parent.paramf.write('\n')
				self.parent.paramf.flush()
			
		# Add line to tree file if it exists
		if self.parent.treef:
			self.parent.treef.write('\ttree rep.%d = %s;\n' % (cycle + 1, cold_chain.tree.makeNewick(self.parent.opts.ndecimals)))
			self.parent.treef.flush()

		# Add line to sitelike file if it exists and if we are saving site-likelihoods
		if self.parent.opts.saving_sitelikes:
			if self.parent.sitelikef:
				cold_chain.likelihood.storeSiteLikelihoods(True)
				cold_chain.likelihood.calcLnL(cold_chain.tree)
				cold_chain.likelihood.storeSiteLikelihoods(False)
				
				patternLnLikes = cold_chain.likelihood.getSiteLikelihoods()
				siteloglikes = [0.0]*self.parent.nchar
				for i,patternLnL in enumerate(patternLnLikes):
					for site in self.parent.siteIndicesForPatternIndex[i]:
						siteloglikes[site] = patternLnL
				for siteLnL in siteloglikes:
					self.parent.sitelikef.write(float_format_str % siteLnL)
				self.parent.sitelikef.write('\n')
				self.parent.sitelikef.flush()
								
	def attemptChainSwap(self, cycle):
		#---+----|----+----|----+----|----+----|----+----|----+----|----+----|
		"""
		Attempts to swap two randomly-chosen chains. Letting
		p_i be posterior density for chain i
		p_j be posterior density for chain j
		theta_i be the parameter vector for chain i
		theta_j be the parameter vector for chain j
		
					   p_i(theta_j)	 p_j(theta_i)
		accept_ratio = -----------	 ------------
					   p_i(theta_i)	 p_j(theta_j)
		
		  L(theta_j)^power_i f(theta_j)	 L(theta_i)^power_j f(theta_i) 
		= -----------------------------	 -----------------------------
		  L(theta_i)^power_i f(theta_i)	 L(theta_j)^power_j f(theta_j) 
		
		= L(theta_i)^(power_j - power_i) L(theta_j)^(power_i - power_j)
		
		"""
		n = len(self.chains)
		if n < 2:
			return
			
		# Figure out which chains to swap
		
		# Note: self.parent is the MCMCImpl object
		r = self.parent._getLot()
		i = r.sampleUInt(n)
		j = r.sampleUInt(n - 1)
		if j >= i:
			j += 1
		else:
			i,j = j,i	# make sure i < j so attempted swaps are consistently in upper triangle of swap table
		
		# Compute acceptance ratio
		cmi = self.chains[i].chain_manager
		power_i = self.chains[i].heating_power
		cmi.refreshLastLnLike()
		lnLi = cmi.getLastLnLike()
		cmi.refreshLastLnPrior()
		lnPriori = cmi.getLastLnPrior()

		cmj = self.chains[j].chain_manager
		power_j = self.chains[j].heating_power
		cmj.refreshLastLnLike()
		lnLj = cmj.getLastLnLike()
		cmj.refreshLastLnPrior()
		lnPriorj = cmj.getLastLnPrior()

		log_accept_ratio = (power_j - power_i)*(lnLi + lnPriori - lnLj - lnPriorj)
		log_u = math.log(r.uniform())
		swap_accepted = False
		self.swap_table[i][j] += 1	# upper triangle
		if log_u < log_accept_ratio:
			swap_accepted = True
			self.setChainPower(i, power_j)
			self.setChainPower(j, power_i)
			self.chains[i], self.chains[j] = self.chains[j], self.chains[i]
			
			self.swap_table[j][i] += 1	# lower triangle
			
		if cycle % 1 == 0:
			print '\n********** Attempting chain swap **********'
			print 'n =',n
			print 'i =',i
			print 'j =',j
			print 'log_accept_ratio =',log_accept_ratio
			print 'log_u			=',log_u
			if swap_accepted:
				print 'swap_accepted	= True'
			else:
				print 'swap_accepted	= False'
			print 'swap_table:'
			for ii in range(n):
				for jj in range(n):
					print '%12d' % self.swap_table[ii][jj],
				print
			print '*********************************************\n'
