from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas import model

class Subset(object):
	def __init__(self):
		self.start = 0
		self.stop = 0
		self.incr = 0
		
	def __call__(self, first_site, last_site, step_size = 1):
		"""
		Creates a range using the supplied first site, last site and step size.
		"""
		self.start = int(first_site)
		self.stop = int(last_site) + 1
		self.incr = int(step_size)
		return range(self.start, self.stop, self.incr)

class Partition(PhycasCommand):
	"""
	The Python class Partition handles the user interface for defining a partition.
	It has no exact counterpart on the C++ side. The class PartitionModel 
	stores partitioning information on the C++ side.
	"""
	def __init__(self):
		args =	(
				('name', 					'mypart', 	'Choose a name for this partitioning scheme'),
				('fix_subset_relrates',		False,		'If True, the vector of subset relative rates will not be modified during the course of an MCMC analysis', BoolArgValidate),
				('subset_relrates',			None,		'The vector of subset relative rates. If None, the vector of subset relative rates be set to a vector of the appropriate dimension consisting of all 1.0', BoolArgValidate),
				('subset_relrates_prior',	None,		'The joint prior distribution for the relative rates of partition subsets. If None, a flat Dirichlet prior of the appropriate dimension will be used.'),
				)
		
		# Specify output options
		PhycasCommand.__init__(self, args, "partition", "Sets up a data partitioning scheme, allowing different models to be applied to different data subsets.")

		# a list containing (name, sitelist, model) tuples
		self.__dict__['subset'] = []
		
		# a list of integers representing the index of the model for each site
		self.__dict__['sitemodel'] = []
		
		# needed for using phycassert
		self.__dict__['cf'] = CommonFunctions(self)
		
	def hidden():
		""" 
		Overrides the PhycasCommand.hidden method to keep Partitions's name from being displayed 
		in the list of classes displayed when users type help. Delete this function, or
		change its return value to False, when it is ready to be advertised.
		"""
		return True
		
	hidden = staticmethod(hidden)
	
	def flatten(self, alist):
		flattened = []
		for s in alist:
			if s.__class__.__name__ == 'list' or s.__class__.__name__ == 'tuple':
				# an element of the supplied list is itself a list or tuple, so recurse
				s = self.flatten(s)
				flattened.extend(s)
			else:
				if s.__class__.__name__ == 'int':
					self.cf.phycassert(s > 0, 'elements of site lists must be greater than 0; you supplied the value %d' % s)
					flattened.append(s)
				else:
					self.cf.phycassert(False, 'elements of site lists must be of type int; you supplied a value of type %s' % s.__class__)
			 
		return flattened
		
	def getSiteModelVector(self):
		if len(self.subset) < 2:
			# user has not specified a partition, so return an empty list (which will signify
			# the default partition)
			self.sitemodel = []
		return self.sitemodel
		
	def getModels(self):
		"""
		Creates and returns a list of just the models stored in subset.
		"""
		return [m for n,s,m in self.subset]
		
	def validate(self, nsites):
		"""
		If user has not specified a partition, create a default partition now using the current model.
		"""
		if len(self.subset) < 2:
			self.sitemodel = []
			self.addSubset(range(1, nsites + 1), model, 'default')
	
	def addSubset(self, sites, model_for_sites, name = None):
		"""
		Applies a model (model_for_sites) to a list of sites (sites). For example,
		
		model.type = 'hky'
		hky = model()
		partition.addSubset(range(1,100,3), hky, 'first')
		
		"""
		# the index of the new subset is simply the length of self.subset
		subset_index = len(self.subset)
		if name is None:
			name = 'subset%d' % (1 + subset_index)
			
		self.cf.output('Processing subset %s...' % name)
			
		# make a copy of the supplied list in order to flatten it (if necessary)
		sitelist = self.flatten(sites)
		sitelist.sort()
		
		# expand sitemodel list if last site in sorted `sites' list is larger
		# than the last site in self.sitemodel
		curr_size = len(self.sitemodel)
		needed_size = sitelist[-1]
		if curr_size < needed_size:
			xtra = [-1]*(needed_size - curr_size)
			self.sitemodel.extend(xtra)
			
		#print 'adding subset named ',name
		#print 'curr_size	=',curr_size
		#print 'needed_size =',needed_size
		#print 'size		=',len(self.sitemodel)

		# add sites in sitelist to sitemodel (the master list)
		for s in sitelist:
			assigned_index = self.sitemodel[s - 1]
			if assigned_index == subset_index:
				self.cf.phycassert(False, 'site %d has already been assigned a model by the current subset (%s)' % (s,name))
			elif assigned_index >= 0:
				assigned_subset = self.subset[assigned_index]
				self.cf.phycassert(False, 'site %d has already been assigned a model by subset %s' % (s,assigned_subset[0]))
			self.sitemodel[s - 1] = subset_index
			
		#raw_input('length of sitemodel = %d' % len(self.sitemodel))
		self.subset.append((name, sitelist, model_for_sites))
		
	def save(self):
		return copy.deepcopy(self)
	
	def __deepcopy__(self, memo):
		# memo is a dictionary used to avoid infinite recursion
		# memo keeps track of objects already copied
		c = memo.get(self)
		if c:
			return c
			
		# self not already copied, so copy it now
		new_partition = Partition()
		new_partition.name = copy.deepcopy(self.name, memo)
		memo[self] = new_partition
		return new_partition
		
	def partitionReport(self):
		"""
		Ensures that no site has more than one model and issues a summary of the 
		partitioning scheme.
		"""
		print 'Assignment of subsets to sites:'
		for i,x in enumerate(self.sitemodel):
			if x < 0:
				nm = 'default'
			else:
				nm = self.subset[x][0]
			#print '--> %6d %6d %s' % (i,x,nm)

	def __call__(self, **kwargs):
		"""
		Calling a Partition object results in partitionReport() being called, which in turn
		does some sanity checks and issues a summary of the partitioning scheme.
		"""
		self.set(**kwargs)
		return self.partitionReport()
