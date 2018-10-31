
import numpy as np
import matplotlib
matplotlib.use('Qt4Agg')
from matplotlib import pyplot as plt
import h5py    
import os
import elephant
from neo.core import SpikeTrain
from quantities import *
from elephant.conversion import BinnedSpikeTrain

v= 4
igaba = 5
iampa = 4
trial_len = 10
bin_length = 0.2 # in seconds (i.e., 100 ms)
overlap = 0
neuron_num = 640
bin_num = int(np.floor((trial_len - bin_length) /(bin_length*(1-overlap)) + 1))
shift = bin_length*(1-overlap)

# mean_corr = np.zeros([7])
# mean_cov = np.zeros([7])
# mean_corr_avg = np.zeros([7])
# k = 0

# for i in range(0,14,2):
#     fn = os.path.expanduser('~/cluster_home/spiking/proc/pmod_spiketimes_ampa%d_nmda%d_gaba%d_v%d.h5') % (iampa,i,igaba,v)
#     f = h5py.File(fn, 'r') 
#     sp_E = np.array(f[list(f.keys())[0]])
#     sp_I = np.array(f[list(f.keys())[1]])
#     f.close()
#     spikes = np.zeros([neuron_num,bin_num])
#     time = [0, bin_length]
#     for ibin in range(0,bin_num):    
#         print("Bin #%d / %d ...") % (ibin, bin_num)
#         for ineuron in range(0,neuron_num):
#             spikes[ineuron,ibin] = sum((sp_E[0][sp_E[1]==ineuron] >= time[0]) & (sp_E[0][sp_E[1]==ineuron] <= time[1]))   
#         time = [time[0]+shift, time[0]+shift+bin_length]
#     mean_corr[k] = sum(sum(np.triu(np.corrcoef(spikes),1)))/((neuron_num*neuron_num-neuron_num)/2)
#     k = k + 1
    
    


# -------------------------------------
# RASTER PLOTS
# -------------------------------------
v= 4

for igaba in range(7,16):
    for iampa in range(0,16):
        print "Processing gaba%d, ampa%d..." % (igaba,iampa)
        fig, axs = plt.subplots(7,1, figsize=(9,7))
        k = 0
        spike_rate = np.zeros([7])
        for i in range(0,14,2):
            fn = os.path.expanduser('~/cluster_home/spiking/proc/pmod_spiketimes_ampa%d_nmda%d_gaba%d_v%d.h5') % (iampa,i,igaba,v)
            f = h5py.File(fn, 'r') 
            sp_E = np.array(f[list(f.keys())[0]])
            sp_I = np.array(f[list(f.keys())[1]])
            spike_rate[k] = float(len(sp_E[0]))/(640.0*trial_len)
            if spike_rate[k] > 60:
                break
            f.close()
            # COMPUTE SPIKE COUNT CORRELATIONS
            print("Computing spike count correlations...")
            spikes = dict()
            first_spike = 1 # start analysis at 1s
            for ineuron in range(0,neuron_num):
                spikes[ineuron] = SpikeTrain(sp_E[0][(sp_E[1]==ineuron) & (sp_E[0]>first_spike)]*s, t_start = 1.0, t_stop = 10.0)
            st = []
            subsamp = 10
            matidx = np.triu_indices(len(range(0,len(spikes),subsamp)),1)
            for isp in range(0,len(spikes),subsamp):
                st.append(spikes[isp])
            sts=BinnedSpikeTrain(st, binsize=10*ms)
            corr=elephant.spike_train_correlation.corrcoef(sts)
            mean_corr = corr[matidx].mean()
            print("Computing spike count correlations... Done")
            print mean_corr
            # time = [0, bin_length]
            # spikes = np.zeros([neuron_num,bin_num])
            # for ibin in range(0,bin_num):    
            #     print("Bin #%d / %d ...") % (ibin, bin_num)
            #     for ineuron in range(0,neuron_num):
            #         spikes[ineuron,ibin] = sum((sp_E[0][sp_E[1]==ineuron] >= time[0]) & (sp_E[0][sp_E[1]==ineuron] <= time[1]))   
            #     time = [time[0]+shift, time[0]+shift+bin_length]    
            # mean_corr = sum(sum(np.triu(np.corrcoef(spikes),1)))/((neuron_num*neuron_num-neuron_num)/2)          

            fig.add_axes(axs[k])
            plt.plot(sp_E[0][(sp_E[0]<1.8) & (sp_E[0]>1.2)],sp_E[1][(sp_E[0]<1.8) & (sp_E[0]>1.2)],'ro',markersize=1)
            plt.xlabel('Time [in s]'); plt.ylabel('%.1f Hz' %(spike_rate[k]))
            plt.text(1, 1, 'r = %.2f' % mean_corr, fontsize=10, transform=axs[k].transAxes, rotation=270)
            k = k + 1     

        plt.savefig('pmod_spiketimes_ampa%d_nmda%d_gaba%d_v%d.png' % (iampa,i,igaba,v), dpi=300) 
        plt.close()



# for igaba in range(6,16):
#     for iampa in range(0,16):
#         print "Processing gaba%d, ampa%d..." % (igaba,iampa)
#         mean_corr = np.zeros([7])
#         mean_cov = np.zeros([7])
#         mean_corr_avg = np.zeros([7])
#         k = 0

#         for i in range(0,14,2):
#             fn = os.path.expanduser('~/cluster_home/spiking/proc/pmod_spiketimes_ampa%d_nmda%d_gaba%d_v%d.h5') % (iampa,i,igaba,v)
#             f = h5py.File(fn, 'r') 
#             sp_E = np.array(f[list(f.keys())[0]])
#             sp_I = np.array(f[list(f.keys())[1]])
#             f.close()
#             spikes = np.zeros([neuron_num,bin_num])
#             time = [0, bin_length]
#             for ibin in range(0,bin_num):    
#                 print("Bin #%d / %d ...") % (ibin, bin_num)
#                 for ineuron in range(0,neuron_num):
#                     spikes[ineuron,ibin] = sum((sp_E[0][sp_E[1]==ineuron] >= time[0]) & (sp_E[0][sp_E[1]==ineuron] <= time[1]))   
#                 time = [time[0]+shift, time[0]+shift+bin_length]
#             mean_corr[k] = sum(sum(np.triu(np.corrcoef(spikes),1)))/((neuron_num*neuron_num-neuron_num)/2)
#             k = k + 1
            
#             np.savetxt('test.txt', mean_corr, delimiter=',')
            
            


# # -------------------------------------
# # SPIKE TIME TILING COEFFICIENT
# # -------------------------------------
# k = 0
# neuron_num = 1600
# mean_corr = np.zeros([5])
# mean_corr_avg = np.zeros([5])

# for i in range(0,5,1):
#     fn = os.path.expanduser('~/Dropbox/projects/phd/spiking/proc/pmod_spiketimes_ampa%d_nmda%d_gaba%d_v%d.h5') % (iampa,i,igaba,v)
#     f = h5py.File(fn, 'r') 
#     sp_E = np.array(f[list(f.keys())[0]])
#     sp_I = np.array(f[list(f.keys())[1]])
#     f.close()
#     trial_len = 3
#     bin_length = 0.1 # in seconds (i.e., 100 ms)
#     overlap = 0.9
#     bin_num = int(np.floor((trial_len - bin_length) /(bin_length*(1-overlap)) + 1))
#     time = [0, bin_length]
#     shift = bin_length*(1-overlap)
#     spikes = dict()
#     for ineuron in range(0,neuron_num):
#         spikes[ineuron] = SpikeTrain(sp_E[0][sp_E[1]==ineuron]*s, t_stop = 3.0)
#     spike_corr = np.zeros([50,50])
#     for ii in range(0,50):
#         print ii
#         for jj in range(0,50):
#             spike_corr[ii,jj] = elephant.spike_train_correlation.spike_time_tiling_coefficient(spikes[ii],spikes[jj],dt=np.array(0.005) * s)
#     mean_corr[i] = spike_corr.mean()
#     del(spikes)





