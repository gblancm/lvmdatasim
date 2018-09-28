class IFUmodel(object):
    """Read an existing IFU model stored in the data directory as a pickle
    'x'= x position in focal plane
    'y'= y position in focal plane
    'r'= lens radius defined as ???
    'ID'= lens id
    """
    def __init__ (self, ifuname):
        self.lensKernel = None
        self.lensx = None
        self.lensy = None
        self.lensr = None
        self.lensID= None