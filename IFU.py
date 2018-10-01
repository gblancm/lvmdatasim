class IFU(object):
    """Read an existing IFU model stored in the data directory as a pickle
    'ID'= lens id
    'x'= x position in focal plane
    'y'= y position in focal plane
    'r'= lens radius defined as ???
    """
    def __init__ (self, ifuname):
        """ Initialize the IFU class"""
        self.lensKernel = None
        self.lensx = None
        self.lensy = None
        self.lensr = None
        self.lensID= None
