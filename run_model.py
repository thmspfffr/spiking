'''
run_model.py
---------------------
Simulate leaky integrate-and-fire network with excitatory and inhibitory neurons
Parameters are identical (where possible) to Wang (2002) Neuron
Code modified from Wimmer et al. (2015) Nat Comm

See model code in integration_circuit_mod.py

Last updated 04/10/2018, thmspfffr
---------------------
'''


from brian import *
import numpy
import random as pyrandom
from numpy.random import rand as rand
from numpy.random import randn as randn
from scipy.signal import lfilter
from scipy.signal import welch

from integration_circuit_mod import make_integration_circuit

if __name__ == '__main__':

    #  initialize  
    defaultclock.reinit()
    clear(True) 

    #------------------------------------------------------------------------------ 
    # Simulation parameters 
    #------------------------------------------------------------------------------ 

    # Timing 
    runtime = 3000.0 * ms                        # total simulation time
    inputs = np.linspace(0,2000,2000/200+1)
    inhibition = np.linspace(0.75,1.25,0.5/0.05+1)

    AMPA_mods = np.linspace(0.75,1.25,0.5/0.05+1)
    NMDA_mods = np.linspace(0.75,1.25,0.5/0.05+1)

    # Integration circuit
    resp = np.zeros([len(inputs),len(inhibition)])

    for iinp in range(0,len(inputs)):
        for i_inhibition in range(0,len(inhibition)):

            print "Computing input #%d..." % (iinp)

            #inp = inputs[iinp]
            inp = inputs[1];
            inh = inhibition[10]
            AMPA_mod = AMPA_mods[iinp]
            NMDA_mod = NMDA_mods[i_inhibition]

            Dgroups, Dconnections, Dnetfunctions, subgroups = make_integration_circuit(inp,inh,AMPA_mod,NMDA_mod)

            # get populations from the integrations circuit
            decisionE = Dgroups['DE']
            decisionI = Dgroups['DI']

            # ---- set initial conditions (random)
            decisionE.gen = decisionE.gen * (1 + 0.2 * rand(decisionE.__len__()))
            decisionI.gen = decisionI.gen * (1 + 0.2 * rand(decisionI.__len__()))
            decisionE.V = decisionE.V + rand(decisionE.__len__()) * 2 * mV
            decisionI.V = decisionI.V + rand(decisionI.__len__()) * 2 * mV

            # record spikes of excitatory neurons
            S_DE = SpikeMonitor(decisionE, record=True)
            # record spikes of inhibitory neurons
            S_DI = SpikeMonitor(decisionI, record=True)
            # record instantaneous excitatory populations activity
            R_DE = PopulationRateMonitor(decisionE, bin=5*ms)
            # record instantaneous inhibitory populations activity
            R_DI = PopulationRateMonitor(decisionI, bin=5*ms)

            #------------------------------------------------------------------------------
            # Run the simulation
            #------------------------------------------------------------------------------
            net = Network(Dgroups.values(),  Dconnections.values(),
                            Dnetfunctions, S_DE, S_DI, R_DE, R_DI)
            net.prepare()
            net.run(runtime)      

            resp[iinp,i_inhibition] = R_DE.rate.mean()
            
            # Obtain population firing rate frequency
            f, Pxx = welch(R_DE.rate, fs=200, window='hanning', nperseg=100,noverlap=0.5, nfft=None, detrend=False)

            fn = '/Users/tpfeffer/Dropbox/projects/phd/spiking/proc/pmod_spiking_fft_inp%d_inh%d.npy' % (iinp,i_inhibition)
            savez(fn,freqs=f,pow=Pxx)

    save('/Users/tpfeffer/Dropbox/projects/phd/spiking/proc/pmod_spiking.npy',resp)

