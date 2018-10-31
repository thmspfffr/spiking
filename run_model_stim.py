'''
Model with current parameter settings gives
roughly a population firing rate of 3.9 Hz
as well as low spike correlations of 0.13
'''

from brian import *
import numpy 
from numpy.random import rand as rand
from numpy.random import randn as randn
import h5py
import os.path
from subprocess import call
from sys import platform

if platform == 'darwin':
    root_dir = '~/Dropbox/projects/phd/spiking/'
else:
    root_dir = '/home/tpfeffer/spiking/'

from integration_circuit_mod import make_integration_circuit

v = 5

if __name__ == '__main__':

    #------------------------------------------------------------------------------ 
    # Simulation parameters 
    #------------------------------------------------------------------------------ 

    # Timing 
    runtime = 2500.0 * ms                        # total simulation time
    inputs = 200 + np.linspace(0,1,1/0.1+1) # note that input 200 is baseline
    
    AMPA_mods = np.linspace(0.25,4,3.75/0.25+1)[12]
    NMDA_mods = np.linspace(0.25,4,3.75/0.25+1)[9]
    GABA_mods = np.linspace(2.7,3.5,0.8/0.1+1)

    # Integration circuit
    resp = np.zeros([AMPA_mods.size, NMDA_mods.size, GABA_mods.size])
    mean_corr = np.zeros([AMPA_mods.size, NMDA_mods.size, GABA_mods.size])
    for iinp in range(0,inputs.size):
        
        for igaba in range(0,GABA_mods.size):
            fn = os.path.expanduser(root_dir + 'proc/pmod_stim_spiketimes_inp%d_igaba%d_v%d_processing.txt') % (iinp,igaba,v)
            if os.path.isfile(fn)==False:
                call(['touch', fn])
            else:
                continue

            #  initialize  
            defaultclock.reinit()
            clear(True) 

            print("Computing INPUT #%d and GABA #%d ...") % (iinp, igaba)
            # see murphy & miller (2003) j neurosci for contrast-response curve
            inp = 2000 * (inputs[iinp]**1.2 / (inputs[iinp]**1.2 + 0.133**1.2)) 
            inh = 1
            AMPA_mod = AMPA_mods
            NMDA_mod = NMDA_mods
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

            print("Saving output...")
            hf = h5py.File(os.path.expanduser(root_dir + 'proc/pmod_stim_spiketimes_igaba%d_inp%d_v%d.h5') % (igaba, iinp, v), 'w')
            hf.create_dataset('spt_E', data=spt_E)
            hf.create_dataset('spt_I', data=spt_I)
            hf.close()
        


