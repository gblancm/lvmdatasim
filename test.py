import lvmdatasim

input='./testInput/small.sii.fits'

# Define LVMSimulator object
mysim=lvmdatasim.LVMSimulator(input=input, telescopeName='LVM160-SCI-S', psfModel=1.0, inputType='fitscube', fluxType='flux', saveConvCube=True)

# Set simulation parameters that are different from defaults
mysim.simparam['ra']=77.4873
mysim.simparam['dec']=-68.898
mysim.simparam['theta']=0.0

# # Set wavelength range for image
# mysim.simparam['minWave']=4500.
# mysim.simparam['maxWave']=7100.
# mysim.simparam['waveStep']=3.

# Run simulatioon

mysim.simulate()

print("done")

