from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas import mcmc

class SS(PhycasCommand):
	def __init__(self):
		args = (   ("nbetavals", 21, "The number of values beta will take on during the run; for example, if this value is 4, then beta will take on these values: 1, 2/3, 1/3, 0", IntArgValidate(min=1)),
				   ("ti", False, "If True, the marginal likelihood will be estimated using thermodynamic integration and the stepping stone method with reference distribution equal to the prior; if False (the default), the stepping stone method with reference distribution approximating the posterior will be used (this greatly improves the accuracy of the stepping stone method and is strongly recommended).", BoolArgValidate),
				   ("xcycles", 0, "The number of extra cycles (above and beyond mcmc.ncycles) that will be spent exploring the posterior (additional posterior cycles help stepping stone analyses formulate an effective reference distribution).", IntArgValidate()),
				   ("maxbeta", 1.0, "The first beta value that will be sampled.", FloatArgValidate(min=0.0, max=1.0)),
				   ("minbeta", 0.0, "The last beta value that will be sampled.", FloatArgValidate(min=0.0, max=1.0)),
				   ("minsample", 10, "Minimum sample size needed to create a split-specific edge length reference distribution.", IntArgValidate(min=0)),
				   ("shape1", 1.0, "The first shape parameter of the distribution used to determine the beta values to be sampled. This distribution is, confusingly, a Beta distribution. Thus, if both shape1 and shape2 are set to 1, beta values will be chosen at uniform intervals from minbeta to maxbeta.", FloatArgValidate(greaterthan=0.0)),
				   ("shape2", 1.0, "The second shape parameter of the distribution used to determine the beta values to be sampled. This distribution is, confusingly, a Beta distribution. Thus, if both shape1 and shape2 are set to 1, beta values will be chosen at uniform intervals from minbeta to maxbeta.", FloatArgValidate(greaterthan=0.0)),
				)
		# Specify output options
		#self.__dict__["hidden"] = True # hide from main phycas help list of commands until working
		#o = PhycasCommandOutputOptions()
		#o.__dict__["_help_order"] = ["sss"]
		PhycasCommand.__init__(self, args, "ss", "Performs stepping stone method (and also thermodynamic integration) for purposes of estimating the marginal likelihood of the current model.")

		# The data member below is hidden from the user because it overrides something that users should not be able to override 
		self.__dict__["override_fixed_topology_restriction"] = False
		
		# The data members added below are hidden from the user because they are set when the mcmc command runs
		self.__dict__["sampled_likes"] = None
		self.__dict__["sampled_betas"] = None
		
		# The data members added below are hidden from the user because they are for developer use only
		self.__dict__["refdist_definition_file"] = None	# specify file name of file containing referenced distributions (one per line in same order that they are outut in the log file)

	def hidden():
		""" 
		Overrides the PhycasCommand.hidden method to keep SS's name from being displayed 
		in the list of classes displayed when users type help. Delete this function, or
		change its return value to False, when it is ready to be advertised.
		"""
		return False
		
	hidden = staticmethod(hidden)

	def checkSanity(self):
		"""
		Place asserts in this function that should be checked before anything substantive
		is done during a call of a SS object.
		"""
		cf = CommonFunctions(self)
		cf.phycassert(mcmc.ncycles > 0, 'mcmc.ncycles cannot be less than 1 for the stepping-stone method')
		if not self.override_fixed_topology_restriction:
			cf.phycassert(mcmc.fix_topology == True, "mcmc.fix_topology must be True to use the stepping-stone method (we're working on relaxing this requirement)")
		
	def __call__(self, **kwargs):
		self.set(**kwargs)
		self.checkSanity()
		mcmc.doing_steppingstone_sampling = True
		mcmc.ssobj = self
		if self.ti:
			mcmc.ss_heating_likelihood = True
		else:
			mcmc.ss_heating_likelihood = False
		mcmc()
		mcmc.ss_heating_likelihood = False
		mcmc.doing_steppingstone_sampling = False
		mcmc.ssobj = None
		self.sampled_betas = mcmc.ss_sampled_betas
		self.sampled_likes = mcmc.ss_sampled_likes
