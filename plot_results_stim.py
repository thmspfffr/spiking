import scipy
import neo
import elephant
from brian import *
import numpy 
import h5py
from neo.core import SpikeTrain
from quantities import *
from elephant.conversion import BinnedSpikeTrain

inputs      = np.linspace(0,1,1/0.1+1) # note that input 200 is baseline
GABA_mods   = np.linspace(2.7,3.5,0.8/0.1+1)
trial_len   = 2.5
neuron_num  = 640
mean_corr   = np.zeros([inputs.size,GABA_mods.size])
spike_rate  = np.zeros([inputs.size,GABA_mods.size])

for iinp in range(0,inputs.size):
    for igaba in range(0,GABA_mods.size):
        print("Processing input #%d, inh #%d ...") % (iinp, igaba)
        fn = '/Users/tpfeffer/Dropbox/projects/phd/spiking/proc/pmod_stim_spiketimes_igaba%d_inp%d_v3.h5' % (igaba, iinp)
        f = h5py.File(fn, 'r') 
        sp_E = np.array(f[list(f.keys())[0]])
        sp_I = np.array(f[list(f.keys())[1]])
        f.close()   
        spike_rate[iinp,igaba] = float(len(sp_E[0]))/(640.0*trial_len)
        spikes = dict()
        first_spike = 0.25 # start analysis at 1s
        for ineuron in range(0,neuron_num):
            spikes[ineuron] = SpikeTrain(sp_E[0][(sp_E[1]==ineuron) & (sp_E[0]>first_spike)]*s, t_start = 0.0, t_stop = trial_len)
        st = []
        subsamp = 5
        matidx = np.triu_indices(len(range(0,len(spikes),subsamp)),1)
        for isp in range(0,len(spikes),subsamp):
            st.append(spikes[isp])
        sts=BinnedSpikeTrain(st, binsize=10*ms)
        print("Compute correlation...")
        corr=elephant.spike_train_correlation.corrcoef(sts)
        mean_corr[iinp,igaba] = corr[matidx].mean()


def fun(x,spike_rate):
    return np.array(spike_rate - (x[3]* (inputs**x[0] / (inputs**x[0] + x[1]**x[0])) + x[2]) )
def fun_plot(x):
    return np.array((x[3]* (np.linspace(0,1,1/0.001+1)**x[0] / (np.linspace(0,1,1/0.001+1)**x[0] + x[1]**x[0])) + x[2]) )

xscale1 = scipy.optimize.leastsq(fun, [1, 1, 1, 1], args=(spike_rate[:,0]))
xscale2 = scipy.optimize.leastsq(fun, [0.3, 1, 1, 50], args=(spike_rate[:,1]))
xscale3 = scipy.optimize.leastsq(fun, [1, 1, 1, 50], args=(spike_rate[:,2]))
xscale4 = scipy.optimize.leastsq(fun, [1, 1, 1, 50], args=(spike_rate[:,3]))
xscale5 = scipy.optimize.leastsq(fun, [1, 1, 1, 50], args=(spike_rate[:,4]))
xscale6 = scipy.optimize.leastsq(fun, [1, 1, 1, 50], args=(spike_rate[:,5]))

fig = plt.figure(1,figsize=(4,6))

plt.plot(np.linspace(0,1,1/0.1+1), spike_rate[:,0],'ro')
plt.plot(np.linspace(0,1,1/0.001+1), fun_plot(xscale1[0][:]),'r')
plt.plot(np.linspace(0,1,1/0.1+1), spike_rate[:,1],'ro')
plt.plot(np.linspace(0,1,1/0.001+1), fun_plot(xscale2[0][:]),'r')
plt.plot(np.linspace(0,1,1/0.1+1), spike_rate[:,2],'ro')
plt.plot(np.linspace(0,1,1/0.001+1), fun_plot(xscale3[0][:]),'r')
plt.plot(np.linspace(0,1,1/0.1+1), spike_rate[:,3],'bs')
plt.plot(np.linspace(0,1,1/0.001+1), fun_plot(xscale4[0][:]),'b')
plt.plot(np.linspace(0,1,1/0.1+1), spike_rate[:,4],'bs')
plt.plot(np.linspace(0,1,1/0.001+1), fun_plot(xscale5[0][:]),'b')
plt.plot(np.linspace(0,1,1/0.1+1), spike_rate[:,5],'bs')
plt.plot(np.linspace(0,1,1/0.001+1), fun_plot(xscale6[0][:]),'b')

plt.xlabel("Contrast (%)")
plt.ylabel("Firing rate")
plt.show()


