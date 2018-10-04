import lvmdatasim

input='./testInput/small.sii.fits'

mysim=lvmdatasim.LVMSimulator(input=input, telescopeName='LVM160-SCI-S', psfModel=1.0, inputType='fitscube', fluxType='flux', saveConvCube=True)

mysim.simulate()

print("done")

