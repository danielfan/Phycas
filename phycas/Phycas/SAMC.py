from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas import mcmc

class SAMC(PhycasCommand):
	def __init__(self):
		args = (   ("nlevels", 5, "The number of energy levels used to partition the parameter space", IntArgValidate(min=3)),
				   ("hilog", 0.0, "The natural logarithm of the posterior for a relatively good tree under the model to be used for the SAMC analysis.", FloatArgValidate()),
				   ("lolog", 0.0, "The natural logarithm of the posterior for a relatively bad (e.g. random) tree under the model to be used for the SAMC analysis.", FloatArgValidate()),
				)
		# Specify output options
		#self.__dict__["hidden"] = True # hide from main phycas help list of commands until working
		#o = PhycasCommandOutputOptions()
		#o.__dict__["_help_order"] = ["sss"]
		#sss = TextOutputSpec(prefix='marglike', suffix=".sss", help_str="The text file in which all sampled likelihood, prior and working prior values are saved if performing smoothed steppingstone sampling. This file is needed by the sump command to estimate the marginal likelihood.")
		#o.__dict__["sss"] = sss
		PhycasCommand.__init__(self, args, "samc", "Performs a Bayesian analysis using Stochastic Approximation Monte Carlo (SAMC), the goal of which is to search for tree topologies with high marginal posterior probability while avoiding getting trapped in local optima.")

		# The data members added below are hidden from the user because they are set automatically given the settings the user has provided
		self.__dict__["energy_levels"] = None

	def hidden():
		""" 
		Overrides the PhycasCommand.hidden method to keep SS's name from being displayed 
		in the list of classes displayed when users type help. Delete this function, or
		change its return value to False, when it is ready to be advertised.
		"""
		return True
		
	hidden = staticmethod(hidden)

	def checkSanity(self):
		"""
		Place asserts in this function that should be checked before anything substantive
		is done during a call of a SS object.
		"""
		cf = CommonFunctions(self)
		cf.phycassert(mcmc.ncycles > 0, 'mcmc.ncycles cannot be less than 1 for SAMC analyses')
		
	def initLevels(self):
		"""
		Establish nlevels energy levels for the analysis, with level 0 being lowest (below the level set by lolog)
		and level nlevels-1 being highest (above the level set by hilog). The remaining nlevels-3 boundaries
		will be equally spaced (on the log scale) between hilog and lolog.
		"""
		n = self.nlevels
		self.energy_levels = [0.0]*n
		self.energy_levels[0] = lolog
		self.energy_levels[n-1] = hilog
		incr = (hilog - lolog)/float(n-2)
		prevlog = lolog
		for i in range(1, n-1):
			self.energy_levels[i] = prevlog + incr
		
	def __call__(self, **kwargs):
		self.set(**kwargs)
		self.checkSanity()
		self.initLevels()
		mcmc.samcobj = self
		mcmc.doing_samc = True
		mcmc()
		mcmc.doing_samc = False
