
from brian import *
import numpy as np
from matplotlib import pyplot as plt

resp = np.load('/Users/tpfeffer/Dropbox/projects/phd/spiking/proc/pmod_spiking.npy')

colors = np.matlib.repmat(np.linspace(0,0.5,11),3,1)

fig = plt.figure()
ax = plt.subplot()

ax.set_color_cycle(np.transpose(colors))

plt.plot(resp)

plt.xlabel("Excitatory input [Hz]")
plt.xticks(np.arange(0,11),int(np.linspace(0,2000,2000/200+1)))
plt.ylabel("Population firing rate [Hz]")

plt.tick_params(direction='out',length=5)
fig.set_facecolor((1,1,1))

plt.show()

# -------------------------------------
# PLOT RESULTS OF FFT ANALYSIS
# -------------------------------------
pow = np.load('/Users/tpfeffer/Dropbox/projects/phd/spiking/proc/pmod_spiking_fft_inp1_inh1.npy.npz')
plt.plot(np.log10(pow['freqs']),pow['pow'])
pow = np.load('/Users/tpfeffer/Dropbox/projects/phd/spiking/proc/pmod_spiking_fft_inp5_inh5.npy.npz')
plt.plot(np.log10(pow['freqs']),pow['pow'])
pow = np.load('/Users/tpfeffer/Dropbox/projects/phd/spiking/proc/pmod_spiking_fft_inp10_inh10.npy.npz')
plt.plot(np.log10(pow['freqs']),pow['pow'])

plt.show()