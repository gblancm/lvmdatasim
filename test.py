import lvmdatasim
import numpy as np

input='./testInput/small.sii.fits'

# Define LVMSimulator object
mysim=lvmdatasim.LVMSimulator(input=input, telescopeName='LVM160-SCI-S', psfModel=1.0, inputType='fitscube', fluxType='flux', saveConvCube=True)

# Set simulation parameters that are different from defaults
mysim.simparam['ra']=77.4873
mysim.simparam['dec']=-68.898
mysim.simparam['theta']=0.0

# Run simulatioon

mysim.simulate()

def print_median_snr(fiber=0, sim=mysim.simulator):
    for output in sim.camera_output:
        name = output.meta['name']
        pixel_size = output.meta['pixel_size']
        snr = (output['num_source_electrons'][:, fiber] /
               np.sqrt(output['variance_electrons'][:, fiber]))
        print('{0} median SNR = {1:.3f} / {2:.1f}'.format(name, np.median(snr), pixel_size))
        
print_median_snr(fiber=0)

