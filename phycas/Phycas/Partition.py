from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions

class Partition(PhycasCommand):
    def __init__(self):
        args =  (
                ("name", 'mypart', 'Choose a name for this partitioning scheme'),
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
    
    def addSubset(self, sites, model_for_sites, name = None):
        """
        Applies a model (model_for_sites) to a list of sites (sites). 
        """
        subset_index = len(self.subset)
        if name is None:
            name = 'subset%d' % (1 + subset_index)
            
        self.cf.output('Processing subset %s...' % name)
            
        # make a copy of the supplied list in order to flatten it (if necessary)
        sitelist = self.flatten(sites)
        sitelist.sort()
        
        # expand sitemodel list if necessary
        curr_size = len(self.sitemodel)
        needed_size = sitelist[-1]
        if curr_size < needed_size:
            xtra = [-1]*(needed_size - curr_size)
            self.sitemodel.extend(xtra)
            
        #print 'curr_size   =',curr_size
        #print 'needed_size =',needed_size
        #print 'size        =',len(self.sitemodel)

        # add sites in sitelist to sitemodel (the master list)
        for s in sitelist:
            assigned_index = self.sitemodel[s - 1]
            if assigned_index == subset_index:
                self.cf.phycassert(False, 'site %d has already been assigned a model by the current subset (%s)' % (s,name))
            elif assigned_index >= 0:
                assigned_subset = self.subset[assigned_index]
                self.cf.phycassert(False, 'site %d has already been assigned a model by subset %s' % (s,assigned_subset[0]))
            self.sitemodel[s - 1] = subset_index
            
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
        new_partition.name             = copy.deepcopy(self.name, memo)
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
            print '%6d %12d %s' % (i,x,nm)

    def __call__(self, **kwargs):
        """
        Calling a Partition object results in partitionReport() being called, which in turn
        does some sanity checks and issues a summary of the partitioning scheme.
        """
        self.set(**kwargs)
        return self.partitionReport()
