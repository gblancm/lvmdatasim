import hexagonlib as hexlib
import astropy.io.ascii as ascii
import sys

class IFU(object):
    """Read an existing IFU model stored in the data directory as a pickle
    Internal Parameters
    ----------
    ID: array
        array of lens IDs in the focal plane
    blockID: array 
        IFU block number 
    blockN: array
        Fiber ID in the block
    x: array
        x-position in focal plane
    y: array
        y-position in focal plane
    r: array
        lens radius defined as average radius from center to corner
    lensKernel: 2d-array
        Place holder to store a 2d-array lens kernel 
        calculated from the analytic hexagon and 
        the datacube PIXSCALE
    cuber: array
        Sudo radius in cubic coordinates
    ring: array
        array containing index of the hexagonal ring
    Transmission:array


    """
    def __init__ (self, ifuModel):
        """ Initialize the IFU class """
        self.ifuModel = ifuModel
        self.lensKernel = None
        self.lensID = None
        self.blockID= None
        self.blockN = None
        self.CurrentID = 0

        if isinstance(self.ifuModel, (float, int)):
            if int(self.ifuModel) == self.ifuModel:
                self.ifuModel = int(self.ifuModel)            
            self.lensx = None
            self.lensy = None
            self.lensr = None
            self.cubex = None # this will be blank if the ifu model is read from a table
            self.cubey = None # this will be blank if the ifu model is read from a table
            self.ring = None # this MAY be blank if the ifu model is read from a table
            self.Trans = 1.0
            self.hexLayout = hexlib.Layout(hexlib.layout_flat, 1.0, hexlib.Point(0,0))
            self.mkIfu(N=self.ifuModel)
        elif isinstance(self.ifuModel, str):
            try:
                ascii.read(self.ifuModel, format="commented_header", delimiter="\t")
            sys.exit("Reading an ifuModel from a file is not yet defined.")


    def mkIfu(self, N):
        self.mkCubeIfu(N=N)
        self.cubeIfuToXY()

    def mkCubeIfu(self, N=30):
        """
        An IFU of N rings above the central fiber is returned as a list of hexlib.Hex objects.
        Coordinates are in Cube coordinates (3 numbers per hexagon).
        """
        self.lensx = []
        self.lensy = []
        self.lensr = []
        self.lensID= []
        self.cubex = [] 
        self.cubey = []
        self.ring = []
        self.cubeIfu = []

        self.CurrentID = 0
        for dx in range(-N,N+1):
            for dy in range(max(-N, -dx-N), min(N, -dx+N)+1):
                dz = -dx-dy
                self.cubeIfu.append(hexlib.Hex(dx, dy, dz))
                self.ring.append( int(abs(dx)+abs(dy) + abs(dz))/2 )
                self.lensID.append(self.CurrentID)
                self.CurrentID += 1

    def cubeIfuToXY(self):
        """Modify this to calculate and return 2 arrays, not just one.
        xy = xy position in the focal plane stored as Point(x,y)
        cubexy = cube xy position stored as Point(cubex,cubey)
        Then these can be used to calculate arrays of ring number and theta. Those get passed to fiber-ID assinging function
        """
        #self.xyIfu = [hexlib.hex_to_pixel(self.hexLayout, cubeCoord)  for cubeCoord in self.cubeIfu]
        self.xyIfu = []
        for cubeCoord in self.cubeIfu:
            self.xyIfu.append(hexlib.hex_to_pixel(self.hexLayout, cubeCoord))
            (x,y) = hexlib.hex_to_pixel(self.hexLayout, cubeCoord)
            self.lensx.append(x)
            self.lensy.append(y)
            self.lensr.append((x**2+y**2)**0.5)


def main():
    tel = "LVM160-SCI-S"
    ifu = IFU(3)
    print(ifu.ifuModel)


if __name__ == '__main__':
    main()