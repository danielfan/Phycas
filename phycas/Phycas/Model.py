import copy
from phycas.Utilities.PhycasCommand import *
from phycas.ProbDist import Beta, Exponential, InverseGamma, Dirichlet
class Model(PhycasCommand):
	def __init__(self):
		args = ( 
				("type",					   'hky',								 "Can be 'jc', 'hky', 'gtr' or 'codon'", EnumArgValidate(['jc', 'hky', 'gtr', 'codon'])),
				("update_relrates_separately", True,								 "If True, GTR relative rates will be individually updated using slice sampling; if False, they will be updated jointly using a Metropolis-Hastings move (generally both faster and better)."),
				("relrate_prior",			   Dirichlet([1.0,1.0,1.0,1.0,1.0,1.0]), "The joint prior distribution for all six GTR relative rate parameters. Used only if update_relrates_separately is False."),
				("relrate_param_prior",		   Exponential(1.0),					 "The prior distribution for individual GTR relative rate parameters.  Used only if update_relrates_separately is true."),
				("relrates",				   [1.0, 4.0, 1.0, 1.0, 4.0, 1.0] ,		 "The current values for GTR relative rates. These should be specified in this order: A<->C, A<->G, A<->T, C<->G, C<->T, G<->T."),
				("fix_relrates",			   False,								 "If True, GTR relative rates will not be modified during the course of an MCMC analysis", BoolArgValidate),
				("kappa_prior",				   Exponential(1.0),					 "The prior distribution for the kappa parameter in an HKY model"),
				("kappa",					   4.0,									 "The current value for the kappa parameter in an HKY model", FloatArgValidate(greaterthan=0.0)),
				("fix_kappa",				   False,								 "If True, the HKY kappa parameter will not be modified during the course of an MCMC analysis", BoolArgValidate),
				("omega_prior",				   Exponential(20.0),					 "The prior distribution for the omega parameter in a codon model"),
				("omega",					   0.05,								 "The current value for the omega parameter in a codon model", FloatArgValidate(greaterthan=0.0)),
				("fix_omega",				   False,								 "If True, the omega parameter will not be modified during the course of an MCMC analysis", BoolArgValidate),
				("num_rates",				   1,									 "The number of relative rates used for the discrete gamma rate heterogeneity submodel; default is rate homogeneity (i.e. 1 rate)", IntArgValidate(min=1)),
				("gamma_shape_prior",		   Exponential(1.0),					 "The prior distribution for the shape parameter of the gamma among-site rate distribution"),
				("gamma_shape",				   0.5,									 "The current value for the gamma shape parameter", FloatArgValidate(greaterthan=0.0)),
				("fix_shape",				   False,								 "If True, the gamma shape parameter will not be modified during the course of an MCMC analysis", BoolArgValidate),
				("use_inverse_shape",		   False,								 "If True, gamma_shape_prior is applied to 1/shape rather than shape", BoolArgValidate),
				("pinvar_model",			   False,								 "If True, an invariable sites submodel will be applied and the parameter representing the proportion of invariable sites will be estimated", BoolArgValidate),
				("pinvar_prior",			   Beta(1.0, 1.0),						 "The prior distribution for pinvar, the proportion of invariable sites parameter"),
				("pinvar",					   0.2,									 "The current value of pinvar, the proportion of invariable sites parameter", ProbArgValidate()),
				("fix_pinvar",				   False,								 "If True, the proportion of invariable sites parameter (pinvar) will not be modified during the course of an MCMC analysis", BoolArgValidate),
				("update_freqs_separately",	   True,								 "If True, state frequencies will be individually updated using slice sampling; if False, they will be updated jointly using a Metropolis-Hastings move (generally both faster and better)."),
				("state_freq_prior",		   Dirichlet([1.0, 1.0, 1.0, 1.0]),		 "The joint prior distribution for the relative state frequency parameters. Used only if update_freqs_separately is False."),
				("state_freq_param_prior",	   Exponential(1.0),					 "The prior distribution for the individual base frequency parameters; these parameters, when normalized to sum to 1, represent the equilibrium proportions of the nucleotide states. Used only if update_freqs_separately is True."),
				("state_freqs",				   [1.0, 1.0, 1.0, 1.0],				 "The current values for the four base frequency parameters"),
				("fix_freqs",				   False,								 "If True, the base frequencies will not be modified during the course of an MCMC analysis", BoolArgValidate),
				("edgelen_hyperprior",		   InverseGamma(2.1,1.0/1.1),			 "The prior distribution for the hyperparameter that serves as the mean of an Exponential edge length prior. If set to None, a non-hierarchical model will be used with respect to edge lengths. Note that specifying an edge length hyperprior will cause internal and external edge length priors to be Exponential distributions (regardless of what you assign to internal_edgelen_prior, external_edgelen_prior or edgelen_prior)."),
				("fix_edgelen_hyperparam",	   False,								 "If True, the hyperparameter that governs the mean of the Exponential edge length prior will be fixed at the value edgelen_hyperparam.", BoolArgValidate),
				("edgelen_hyperparam",		   0.05,								 "The current value of the edge length hyperparameter - setting this currently has no effect", FloatArgValidate(greaterthan=0.0)),
				("internal_edgelen_prior",	   Exponential(2.0),					 "Can be used to set a prior distribution for internal edges that differs from that applied to external edges. If this is set to something besides None, you should also set external_edgelen_prior appropriately. Setting the edgelen_prior option sets both external_edgelen_prior and internal_edgelen_prior to the same value"),
				("external_edgelen_prior",	   Exponential(2.0),					 "Can be used to set a prior distribution for external edges that differs from that applied to internal edges. If this is set to something besides None, you should also set internal_edgelen_prior appropriately. Setting the edgelen_prior option sets both external_edgelen_prior and internal_edgelen_prior to the same value"),
				("edgelen_prior",			   None,								 "Sets both internal_edgelen_prior and external_edgelen_prior to the supplied value. Use this setting if you want all edges in the tree to have the same prior distribution. Using this setting will overwrite any values previously supplied for internal_edgelen_prior and external_edgelen_prior"),
				("fix_edgelens",			   False,								 "not yet documented", BoolArgValidate),
				)
		PhycasCommand.__init__(self, args, "model", "Defines a substitution model.")

	def saveas(self):
		return copy.deepcopy(self)

	def __deepcopy__(self, memo):
		"""
		Used by the copy module to make a deep copy of a Model object. The memo parameter
		is a dictionary holding objects already copied.
		"""
		# If this object has already been copied (i.e. it can be found in the memo
		# dictionary), simply return it
		c = memo.get(self)
		if c:
			return c

		# Only get here if the Model object needs to be copied
		new_model = Model()
		new_model.type							= copy.deepcopy(self.type, memo)

		new_model.update_relrates_separately	= copy.deepcopy(self.update_relrates_separately, memo)
		new_model.update_freqs_separately		= copy.deepcopy(self.update_freqs_separately, memo)

		new_model.num_rates						= copy.deepcopy(self.num_rates, memo)
		new_model.use_inverse_shape				= copy.deepcopy(self.use_inverse_shape, memo)
		new_model.pinvar_model					= copy.deepcopy(self.pinvar_model, memo)

		new_model.relrate_param_prior			= copy.deepcopy(self.relrate_param_prior, memo)
		new_model.relrate_prior					= copy.deepcopy(self.relrate_prior, memo)
		new_model.omega_prior					= copy.deepcopy(self.omega_prior, memo)
		new_model.kappa_prior					= copy.deepcopy(self.kappa_prior, memo)
		new_model.gamma_shape_prior				= copy.deepcopy(self.gamma_shape_prior, memo)
		new_model.pinvar_prior					= copy.deepcopy(self.pinvar_prior, memo)
		new_model.state_freq_prior				= copy.deepcopy(self.state_freq_prior, memo)
		new_model.state_freq_param_prior		= copy.deepcopy(self.state_freq_param_prior, memo)

		new_model.edgelen_prior					= copy.deepcopy(self.edgelen_prior, memo)
		new_model.internal_edgelen_prior		= copy.deepcopy(self.internal_edgelen_prior, memo)
		new_model.external_edgelen_prior		= copy.deepcopy(self.external_edgelen_prior, memo)

		new_model.edgelen_hyperprior			= copy.deepcopy(self.edgelen_hyperprior, memo)
		new_model.edgelen_hyperparam			= copy.deepcopy(self.edgelen_hyperparam, memo)

		new_model.relrates						= copy.deepcopy(self.relrates, memo)
		new_model.omega							= copy.deepcopy(self.omega, memo)
		new_model.kappa							= copy.deepcopy(self.kappa, memo)
		new_model.gamma_shape					= copy.deepcopy(self.gamma_shape, memo)
		new_model.pinvar						= copy.deepcopy(self.pinvar, memo)
		new_model.state_freqs					= copy.deepcopy(self.state_freqs, memo)

		new_model.fix_relrates					= copy.deepcopy(self.fix_relrates, memo)
		new_model.fix_omega						= copy.deepcopy(self.fix_omega, memo)
		new_model.fix_kappa						= copy.deepcopy(self.fix_kappa, memo)
		new_model.fix_shape						= copy.deepcopy(self.fix_shape, memo)
		new_model.fix_pinvar					= copy.deepcopy(self.fix_pinvar, memo)
		new_model.fix_freqs						= copy.deepcopy(self.fix_freqs, memo)
		
		new_model.fix_edgelens					= copy.deepcopy(self.fix_edgelens, memo)
		new_model.fix_edgelen_hyperparam		= copy.deepcopy(self.fix_edgelen_hyperparam, memo)

		memo[self] = new_model
		return new_model
		
	def __call__(self, **kwargs):
		self.set(**kwargs)
		return self.saveas()
		








