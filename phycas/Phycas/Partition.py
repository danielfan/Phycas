from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions

class Partition(PhycasCommand):
    def __init__(self):
        args =  (
                ("name", 'mypart', 'Choose a name for this partitioning scheme'),
                )
        
        # Specify output options
        PhycasCommand.__init__(self, args, "partition", "Sets up a data partitioning scheme, allowing different models to be applied to different data subsets.")
        

    def hidden():
        """ 
        Overrides the PhycasCommand.hidden method to keep Partitions's name from being displayed 
        in the list of classes displayed when users type help. Delete this function, or
        change its return value to False, when it is ready to be advertised.
        """
        return False
        
    hidden = staticmethod(hidden)
    
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

    def __call__(self, **kwargs):
        # calling a Partition object returns a deep copy of the partition object
        # Not sure this is the best use of the () operator, however.
        self.set(**kwargs)
        return copy.deepcopy(self)
