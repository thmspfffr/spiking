clear fr r

v = 6;
for iinp = 1
  for iampa = 0 : 23
    for igaba = 0 : 23
      
      fr(iampa+1,igaba+1) = h5read(sprintf('~/spiking/proc/pmod_spiketimes_iinp%d_ampa%d_nmda0_gaba%d_v%d.h5',iinp,iampa,igaba,v),'/spt_E_fr');
       r(iampa+1,igaba+1) = h5read(sprintf('~/spiking/proc/pmod_spiketimes_iinp%d_ampa%d_nmda0_gaba%d_v%d.h5',iinp,iampa,igaba,v),'/spt_E_r');

    end
  end
end

%%
addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/

inhib = [17 11 5];
ampa = [5 10 15];

% 
inhib = [17 17 17];
ampa = [5 14 23];

figure; set(gcf,'color','w');

subplot(1,2,1);
imagesc(fr,[0 20]);
axis square
tp_editplots
ylabel('Excitation (g_{EE,AMPA})')
title('Firing rate')
hold on
scatter(inhib(1),ampa(1),20,'markerfacecolor','w','markeredgecolor','k')
scatter(inhib(2),ampa(2),20,'markerfacecolor','w','markeredgecolor','k')
scatter(inhib(3),ampa(3),20,'markerfacecolor','w','markeredgecolor','k')
xlabel('Inhibition (g_{EI,GABA})')

subplot(1,2,2);
imagesc(r,[0 1]);
colormap(plasma)
axis square
title('Spike correlations')
tp_editplots
hold on
scatter(inhib(1),ampa(1),20,'markerfacecolor','w','markeredgecolor','k')
scatter(inhib(2),ampa(2),20,'markerfacecolor','w','markeredgecolor','k')
scatter(inhib(3),ampa(3),20,'markerfacecolor','w','markeredgecolor','k')
xlabel('Inhibition (g_{EI,GABA})')

print(gcf,'-dpdf',sprintf('~/spiking/plots/spiking_ei_gain_v2.pdf'))
%%

clear frinp rinp

v = 6;
for iinp = 0 : 10
  for iampa = 0:23
    for igaba = 0:23
      
      frinp(iinp+1,iampa+1,igaba+1) = h5read(sprintf('~/spiking/proc/pmod_spiketimes_iinp%d_ampa%d_nmda0_gaba%d_v%d.h5',iinp,iampa,igaba,v),'/spt_E_fr');
       rinp(iinp+1,iampa+1,igaba+1) = h5read(sprintf('~/spiking/proc/pmod_spiketimes_iinp%d_ampa%d_nmda0_gaba%d_v%d.h5',iinp,iampa,igaba,v),'/spt_E_r');

    end
  end
end

%%
figure;  set(gcf,'color','w'); hold on
subplot(4,3,1); hold on
plot(frinp(:,ampa(1),inhib(1)),'k-')
plot(frinp(:,ampa(1),inhib(1)-1),'r-')
plot(frinp(:,ampa(1),inhib(1)+1),'b-')
xlabel('Input'); ylabel('Firing rate [Hz]'); tp_editplots
text(5,4,sprintf('r=%.3f',nanmean(rinp(:,ampa(1),inhib(1)))),'color','k')
text(5,5,sprintf('r=%.3f',nanmean(rinp(:,ampa(1),inhib(1)-1))),'color','r')
text(5,3,sprintf('r=%.3f',nanmean(rinp(:,ampa(1),inhib(1)+1))),'color','b')
axis([1 12 0 15])

subplot(4,3,2); hold on
plot(frinp(:,ampa(2),inhib(2)),'k-')
plot(frinp(:,ampa(2),inhib(2)-1),'r-')
plot(frinp(:,ampa(2),inhib(2)+1),'b-')
xlabel('Input'); ylabel('Firing rate [Hz]'); tp_editplots
text(5,4,sprintf('r=%.3f',nanmean(rinp(:,ampa(2),inhib(2)))),'color','k')
text(5,5,sprintf('r=%.3f',nanmean(rinp(:,ampa(2),inhib(2)-1))),'color','r')
text(5,3,sprintf('r=%.3f',nanmean(rinp(:,ampa(2),inhib(2)+1))),'color','b')
axis([1 12 0 15])

subplot(4,3,3); hold on
plot(frinp(:,ampa(3),inhib(3)),'k-')
plot(frinp(:,ampa(3),inhib(3)-1),'r-')
plot(frinp(:,ampa(3),inhib(3)+1),'b-')
xlabel('Input'); ylabel('Firing rate [Hz]'); tp_editplots
text(5,4,sprintf('r=%.3f',nanmean(rinp(:,ampa(3),inhib(3)))),'color','k')
text(5,5,sprintf('r=%.3f',nanmean(rinp(:,ampa(3),inhib(3)-1))),'color','r')
text(5,3,sprintf('r=%.3f',nanmean(rinp(:,ampa(3),inhib(3)+1))),'color','b')
axis([1 12 0 15])
print(gcf,'-dpdf',sprintf('~/spiking/plots/spiking_ei_gain_io_v2.pdf'))

%% FITTING

figure;  set(gcf,'color','w'); hold on
subplot(1,3,1); hold on

pars_init = [1 1 1 50];
y = 0:0.1:1;

x = frinp(1:end,ampa(1),inhib(1));
pars = tp_fitnakarushton(x,y,pars_init);
plot(y,x,'k.','markersize',20); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'k')

x = frinp(1:end,ampa(1),inhib(1)-1);
pars = tp_fitnakarushton(x,y,pars);
plot(y,x,'r.','markersize',20); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'r')

x = frinp(1:end,ampa(1),inhib(1)+1);
pars = tp_fitnakarushton(x,y,pars);
plot(y,x,'b.','markersize',20); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'b')
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1],'xticklabel',num2cell([0 0.2 0.4 0.6 0.8 1]))
xlabel('Contrast'); ylabel('Firing rate [Hz]'); tp_editplots
axis([-0.1 1 0 15])

% second subplot
subplot(1,3,2); hold on

pars_init = [5 0.1 5 50];
y = 0:0.1:1;

x = frinp(1:end,ampa(2),inhib(2));
pars = tp_fitnakarushton(x,y,pars_init);
plot(y,x,'k.','markersize',20); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'k')

x = frinp(1:end,ampa(2),inhib(2)-1);
pars = tp_fitnakarushton(x,y,pars);
plot(y,x,'r.','markersize',20); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'r')

x = frinp(1:end,ampa(2),inhib(2)+1);
pars = tp_fitnakarushton(x,y,pars);
plot(y,x,'b.','markersize',20); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'b')

xlabel('Contrast'); ylabel('Firing rate [Hz]'); tp_editplots
axis([-0.1 1 0 15])
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1],'xticklabel',num2cell([0 0.2 0.4 0.6 0.8 1]))
% third subplot
subplot(1,3,3); hold on

pars_init = [200 2 0.5 50];
y = 0:0.1:1;

x = frinp(1:end,ampa(3),inhib(3));
pars = tp_fitnakarushton(x,y,pars_init);
plot(y,x,'k.','markersize',20); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'k')

x = frinp(1:end,ampa(3),inhib(3)-1);
pars = tp_fitnakarushton(x,y,pars_init);
plot(y,x,'r.','markersize',20); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'r')

x = frinp(1:end,ampa(3),inhib(3)+1);
pars = tp_fitnakarushton(x,y,pars_init);
plot(y,x,'b.','markersize',20); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'b')

xlabel('Contrast'); ylabel('Firing rate [Hz]'); tp_editplots
axis([-0.1 1 0 15])
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1],'xticklabel',num2cell([0 0.2 0.4 0.6 0.8 1]))

print(gcf,'-dpdf',sprintf('~/spiking/plots/spiking_ei_gain_io_v2.pdf'))
