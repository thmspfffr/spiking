'''
SPIKING MODEL
----------------------
Run simulation of recurrent spiking network across wide range of parameters.
The goal of this exhaustive search is to identify parameter combinations that yield
a low spontaenous firing rate (1-5 Hz) and low spike correlations (r < 0.10). 
----------------------
Based on code from Wimmer et al. (2014), initial parameters from Wang (2002).
Adapted by thmspfffr, 10/2018
'''

from brian import *
import numpy 
from numpy.random import rand as rand
from numpy.random import randn as randn
import h5py
import os.path
from subprocess import call
from sys import platform
import elephant
from neo.core import SpikeTrain
from elephant.conversion import BinnedSpikeTrain
import quantities as pq

if platform == 'darwin':
    root_dir = '~/Dropbox/projects/phd/spiking/'
else:
    root_dir = '/home/tpfeffer/spiking/'

from integration_circuit_mod import make_integration_circuit

# Version
v = 6

if __name__ == '__main__':

    #------------------------------------------------------------------------------ 
    # Simulation parameters 
    #------------------------------------------------------------------------------ 
    # Timing 
    runtime = 10000.0 * ms 
    # Inputs: stimululus, AMPA, NMDA, GABA
    inputs      = np.linspace(0,1.2,1.2/0.1+1) 
    AMPA_mods   = np.linspace(0.25,5,3.75/0.1+1)
    NMDA_mods   = np.linspace(0.25,5,3.75/0.1+1)
    GABA_mods   = np.linspace(0.25,5,3.75/0.1+1)
    # preallocate
    resp = np.zeros([len(AMPA_mods), len(NMDA_mods), len(GABA_mods), len(inputs)])
    mean_corr = np.zeros([len(AMPA_mods), len(NMDA_mods), len(GABA_mods), len(inputs)])

    # Loop through exp parameters
    for igaba in range(0,len(GABA_mods)): 
        for iinp in range(0,inputs.size):
            for iampa in range(0,len(AMPA_mods)):
                for inmda in range(0,len(NMDA_mods)):

                    fn = os.path.expanduser(root_dir + 'proc/pmod_spiketimes_iampa%d_inmda%d_gaba%d_inp%d_v%d_processing.txt') % (iampa, inmda, igaba, iinp,v)
                    if os.path.isfile(fn)==False:
                        call(['touch', fn])
                    else:
                        continue
                    #  initialize  
                    defaultclock.reinit()
                    clear(True) 

                    print("Computing INPUT%d, GABA%d, AMPA%d and NMDA%d ...") % (iinp, igaba,iampa,inmda)

                    inp = 2000 * (inputs[iinp]**1.2 / (inputs[iinp]**1.2 + 0.133**1.2)) 
                    inh = 1
                    AMPA_mod = AMPA_mods[iampa]
                    NMDA_mod = NMDA_mods[inmda]
                    GABA_mod = GABA_mods[igaba]
    
                    Dgroups, Dconnections, Dnetfunctions, subgroups = make_integration_circuit(inp,GABA_mod,AMPA_mod,NMDA_mod)

                    # get populations from the integrations circuit
                    decisionE = Dgroups['DE']
                    decisionI = Dgroups['DI']

                    # ---- set initial conditions (random)
                    decisionE.gen = decisionE.gen * (1 + 0.2 * rand(decisionE.__len__()))
                    decisionI.gen = decisionI.gen * (1 + 0.2 * rand(decisionI.__len__()))
                    decisionE.V = decisionE.V + rand(decisionE.__len__()) * 2 * mV
                    decisionI.V = decisionI.V + rand(decisionI.__len__()) * 2 * mV

                    # record spikes of excitatory neurons
                    Sp_E = SpikeMonitor(decisionE, record=True)
                    # record spikes of inhibitory neurons
                    Sp_I = SpikeMonitor(decisionI, record=True)
                    # record instantaneous excitatory populations activity
                    R_E = PopulationRateMonitor(decisionE, bin=5*ms)
                    # record instantaneous inhibitory populations activity
                    R_I = PopulationRateMonitor(decisionI, bin=5*ms)
                    # record voltage
                    Vm_E = StateMonitor(decisionE, 'V', record=True)
                    Vm_I = StateMonitor(decisionI, 'V', record=True)

                    #------------------------------------------------------------------------------
                    # Run the simulation
                    #------------------------------------------------------------------------------
                    print("Running simulation...")
                    net = Network(Dgroups.values(),  Dconnections.values(), Dnetfunctions, Sp_E, Sp_I, R_E, R_I, Vm_E, Vm_I)
                    net.prepare()
                    net.run(runtime) 
                    
                    # convert to array and save output has .h5            
                    spt_E = []; spt_E_idx = []
                    spt_I = []; spt_I_idx = []

                    for ineuron in range(0,len(Sp_E.spiketimes)):
                        spt_E = np.append(spt_E,Sp_E.spiketimes.values()[ineuron], axis=None)
                        spt_E_idx = np.append(spt_E_idx,np.matlib.repmat(ineuron,1,len(Sp_E.spiketimes.values()[ineuron])), axis=None)

                    for ineuron in range(0,len(Sp_I.spiketimes)):
                        spt_I = np.append(spt_I,Sp_I.spiketimes.values()[ineuron], axis=None)
                        spt_I_idx = np.append(spt_I_idx,np.matlib.repmat(ineuron,1,len(Sp_I.spiketimes.values()[ineuron])), axis=None)

                    spt_E = np.vstack((spt_E,spt_E_idx))
                    spt_I = np.vstack((spt_I,spt_I_idx))

                    # COMPUTE SPIKE COUNT CORRELATIONS
                    print("Computing spike count correlations...")
                    spikes = dict()
                    first_spike = 1 # start analysis at 1s
                    for ineuron in range(0,640):
                        spikes[ineuron] = SpikeTrain(spt_E[0][(spt_E[1]==ineuron) & (spt_E[0]>first_spike)]*pq.s, t_start = 1.0, t_stop = 10.0)
                    st = []
                    subsamp = 1
                    matidx = np.triu_indices(len(range(0,len(spikes),subsamp)),1)
                    for isp in range(0,len(spikes),subsamp):
                        st.append(spikes[isp])
                    sts=BinnedSpikeTrain(st, binsize=10*pq.ms)
                    corr=elephant.spike_train_correlation.corrcoef(sts)
                    mean_corr = corr[matidx].mean()
                    resp = float(len(spt_E[0]))/(640.0*10)
                    print("Saving output...")
                    hf = h5py.File(os.path.expanduser(root_dir + 'proc/pmod_spiketimes_iinp%d_ampa%d_nmda%d_gaba%d_v%d.h5') % (iinp, iampa, inmda, igaba, v), 'w')
                    hf.create_dataset('spt_E', data=spt_E)
                    hf.create_dataset('spt_I', data=spt_I)
                    hf.create_dataset('spt_E_r', data=mean_corr)
                    hf.create_dataset('spt_E_fr', data=resp)
                    hf.close()
                    

