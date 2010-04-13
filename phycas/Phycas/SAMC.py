from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas import mcmc

class SAMC(PhycasCommand):
	def __init__(self):
		args = (   ("nlevels",                 50, "The number of energy levels used to partition the parameter space", IntArgValidate(min=2)),
				   ("hilog",                  0.0, "The natural logarithm of the posterior for a relatively good tree under the model to be used for the SAMC analysis.", FloatArgValidate()),
				   ("lolog",                  0.0, "The natural logarithm of the posterior for a relatively bad (e.g. random) tree under the model to be used for the SAMC analysis.", FloatArgValidate()),
				   ("gain_t0",               50.0, "The numerator in the gain factor update formula. Larger values maintain larger gain factors longer.", FloatArgValidate(greaterthan=0.0)),
				   ("gain_eta",               0.6, "The power in the gain factor update formula. Must be greater than 0.5 and less than or equal to 1.0.", FloatArgValidate(greaterthan=0.5,max=1.0)),
				   ("likelihood_only",      False, "If True, SAMC analysis will use an improper prior such that the posterior is equivalent to the normalized likelihood surface. If False, SAMC will be based on the usual posterior distribution (with a proper joint prior).", BoolArgValidate),
				   ("reference_tree_source", None, "A TreeCollection that will serve as the source of the reference tree topology. The reference tree should represent the best tree topology known for the current model and data set. Specifying a reference tree makes it possible to determine if and when SAMC finds that tree topology. If a string is passed in, it is interpreted as a the path to a file with trees.", TreeSourceValidate),
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
		cf.phycassert(self.hilog >= self.lolog, 'samc.hilog must be greater than or equal to samc.lolog')
		
	def initLevels(self):
		"""
		Establish nlevels energy levels for the analysis, with level 0 being lowest (below the level set by lolog)
		and level nlevels-1 being highest (above the level set by hilog). The remaining nlevels-3 boundaries
		will be equally spaced (on the log scale) between hilog and lolog.
		"""
		wider_at_bottom  = 0
		wider_at_top     = 1
		equally_spaced   = 2
		up_from_bottom   = 3
		
		choice = equally_spaced
		
		if choice == wider_at_bottom:
			m = self.nlevels - 2
			self.energy_levels = [0.0]*(m+1)
			a = 2.0*(self.hilog - self.lolog - float(m))/float(m*(m-1))
			self.energy_levels[m] = self.hilog
			prevlog = self.energy_levels[m]
			for i in range(m):
				currlog = prevlog - (1.0 + a*float(i))
				prevlog = self.energy_levels[m-i-1] = currlog
		elif choice == wider_at_top:
			n = self.nlevels
			self.energy_levels = [0.0]*(n-1)
			incr = (self.hilog - self.lolog)/float(n-2)
			#raw_input('incr = %g' % incr)
			prevlog = self.lolog
			for i in range(n-1):
				self.energy_levels[i] = prevlog + incr*float(i)
		elif choice == equally_spaced:
			m = self.nlevels - 2
			incr = (self.hilog - self.lolog)/float(m)
			print 'incr = %g' % incr
			self.energy_levels = [0.0]*(m + 1)
			for i in range(m + 1):
				self.energy_levels[i] = self.lolog + incr*float(i)
		elif choice == up_from_bottom:
			n = self.nlevels
			self.energy_levels = [0.0]*(n - 1)
			for i in range(n - 1):
				self.energy_levels[i] = self.lolog + float(i)			
		
	def __call__(self, **kwargs):
		self.set(**kwargs)
		self.checkSanity()
		self.initLevels()
		mcmc.samcobj = self
		mcmc.doing_samc = True
		mcmc()
		mcmc.doing_samc = False
