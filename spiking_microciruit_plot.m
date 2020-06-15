clear fr r

% From python script
inputs    = linspace(0,1.2,13) 
AMPA_mods   = linspace(0.2,5,24);
NMDA_mods   = 1;
GABA_mods   = linspace(0.2,5,24) 

v = 1;
for itr = 1
  for iinp = 0 : 12
    
    iinp
    for iampa = 0 : 23
      for igaba = 0 : 23

        fr(iinp+1,iampa+1,igaba+1,itr+1) = h5read(sprintf('~/spiking/proc/spiking_circuitmodel_iinp%d_ampa%d_nmda0_gaba%d_tr%d_v%d.h5',iinp,iampa,igaba,itr,v),'/spt_E_fr');
         r(iinp+1,iampa+1,igaba+1,itr+1) = h5read(sprintf('~/spiking/proc/spiking_circuitmodel_iinp%d_ampa%d_nmda0_gaba%d_tr%d_v%d.h5',iinp,iampa,igaba,itr,v),'/spt_E_r');

      end
    end
  end

end

%%
addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/

inhib = [17 17];
ampa = [5 23];

figure; set(gcf,'color','w');

subplot(1,2,1);
imagesc(squeeze(fr(4,:,:)),[0 3]);
% axis square
tp_editplots
ylabel('Excitation (g_{EE,AMPA})')
title('Firing rate')
hold on
for i = 1 : 2
scatter(inhib(i),ampa(i),20,'markerfacecolor','w','markeredgecolor','k')
end
xlabel('Inhibition (g_{EI,GABA})')
axis([12.5 24.5 0.5 24.5])
set(gca,'XTick',1:2:length(GABA_mods),'XTickLabels',num2cell(round(GABA_mods(1:2:end)*10)/10),'ydir','normal')
set(gca,'YTick',1:2:length(AMPA_mods),'YTickLabels',num2cell(round(AMPA_mods(1:2:end)*10)/10),'ydir','normal')

subplot(1,2,2);
imagesc(squeeze(r(4,:,:)),[0 0.2]);
colormap(plasma)
% axis square
title('Spike correlations')
tp_editplots
hold on
for i = 1 : 2
scatter(inhib(i),ampa(i),20,'markerfacecolor','w','markeredgecolor','k')
end
xlabel('Inhibition (g_{EI,GABA})')
axis([12.5 24.5 0.5 24.5])
set(gca,'XTick',1:2:length(GABA_mods),'XTickLabels',num2cell(round(GABA_mods(1:2:end)*10)/10),'ydir','normal')
set(gca,'YTick',1:2:length(AMPA_mods),'YTickLabels',num2cell(round(AMPA_mods(1:2:end)*10)/10),'ydir','normal')

print(gcf,'-dpdf',sprintf('~/spiking/plots/spiking_circuitmodel_ei_gain_v%d.pdf',v))

%% FITTING
% Fit Naka-Rushton function to simulated data, using linear least squares
% Naka-Rushton: R = Rmax * (C^n / (C^n + C50^n)) + S
% pars(1) = Rmax; pars(2) = C50; pars(3) = n; pars(4) = S
    
marker = 8;
figure;  set(gcf,'color','w'); hold on
subplot(1,5,1); hold on

pars_init = [10 1 1 1];
y = 0:0.1:1;
start_val = 1;

x = fr(start_val:end,ampa(1),inhib(1));
pars = tp_fitnakarushton(x,y,pars_init);
plot(y,x,'.','color',[0.7 0.7 0.7],'markersize',marker); hold on
val = pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3)));
slp(1) = val(end)-val(2);
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'color',[0.7 0.7 0.7])
Rmax(1,1) = pars(1);
C50(1,1) = pars(2);

x = fr(start_val:end,ampa(1),inhib(1)-1);
pars = tp_fitnakarushton(x,y,pars);
plot(y,x,'r.','markersize',marker); hold on
val = pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3)));
slp(2) = val(end)-val(2);
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'r')
Rmax(1,2) = pars(1);
C50(1,2) = pars(2);

x = fr(start_val:end,ampa(1),inhib(1)+1);
pars = tp_fitnakarushton(x,y,pars);
plot(y,x,'b.','markersize',marker); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'b')
val = pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3)));
slp(3) = val(end)-val(2);
set(gca,'xtick',[0 0.5 1],'xticklabel',num2cell([0 0.5 1]))
xlabel('Contrast'); ylabel('Firing rate [Hz]'); tp_editplots
axis([-0.1 1 0 15]); axis square 
Rmax(1,3) = pars(1);
C50(1,3) = pars(2);

% second subplot
subplot(1,5,2); hold on
pars_init = [10 0.1 5 50];

x = fr(start_val:end,ampa(2),inhib(2));
pars = tp_fitnakarushton(x,y,pars_init);
plot(y,x,'.','color',[0.7 0.7 0.7],'markersize',marker); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'color',[0.7 0.7 0.7])
axis square; tp_editplots
Rmax(2,1) = pars(1);
C50(2,1) = pars(2);

x = fr(start_val:end,ampa(2),inhib(2)-1);
pars = tp_fitnakarushton(x,y,pars);
plot(y,x,'r.','markersize',marker); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'r')
axis square ; tp_editplots
Rmax(2,2) = pars(1);
C50(2,2) = pars(2);

x = fr(start_val:end,ampa(2),inhib(2)+1);
pars = tp_fitnakarushton(x,y,pars);
plot(y,x,'b.','markersize',marker); hold on
plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'b')
axis square ; tp_editplots
Rmax(2,3) = pars(1);
C50(2,3) = pars(2);

xlabel('Contrast'); %ylabel('Firing rate [Hz]'); tp_editplots
axis([-0.1 1 0 15])
set(gca,'xtick',[0 0.5 1],'xticklabel',num2cell([0 0.5 1]))



print(gcf,'-depsc2',sprintf('~/spiking/plots/spiking_ei_gain_io_v%d.eps',v))

%% FITTING LOGISTIC FUNCTION
% pars(1): maximum value
% pars(2): growth/slope
% pars(3): x0

marker = 8;
figure;  set(gcf,'color','w'); hold on
subplot(1,5,1); hold on

pars_init = [1 1 1];
y = 0:0.1:1;
start_val = 1;

x = fr(start_val:end,ampa(1),inhib(1));
pars = tp_fitlogistic(x,y,pars_init);
plot(y,x,'.','color',[0.7 0.7 0.7],'markersize',marker); hold on
val = pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3)));
plot(y,pars(1) ./ (1 + exp (-pars(2)*(y - pars(3)))),'color',[0.7 0.7 0.7])
L(1,1) = pars(1);
slope(1,1) = pars(2);

x = fr(start_val:end,ampa(1),inhib(1)-1);
pars = tp_fitlogistic(x,y,pars_init);
plot(y,x,'.','color',[0.7 0.7 0.7],'markersize',marker); hold on
val = pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3)));
plot(y,pars(1) ./ (1 + exp (-pars(2)*(y - pars(3)))),'color',[0.7 0.7 0.7])
L(1,2) = pars(1);
slope(1,2) = pars(2);

x = fr(start_val:end,ampa(1),inhib(1)+1);
pars = tp_fitlogistic(x,y,pars_init);
plot(y,x,'.','color',[0.7 0.7 0.7],'markersize',marker); hold on
val = pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3)));
plot(y,pars(1) ./ (1 + exp (-pars(2)*(y - pars(3)))),'color',[0.7 0.7 0.7])
L(1,3) = pars(1);
slope(1,3) = pars(2);
set(gca,'xtick',[0 0.5 1],'xticklabel',num2cell([0 0.5 1]))
xlabel('Contrast'); ylabel('Firing rate [Hz]'); tp_editplots
axis([-0.1 1 0 15]); axis square 
% 
% % second subplot
% subplot(1,5,2); hold on
% pars_init = [10 0.1 5 50];
% 
% x = fr(start_val:end,ampa(2),inhib(2));
% pars = tp_fitnakarushton(x,y,pars_init);
% plot(y,x,'.','color',[0.7 0.7 0.7],'markersize',marker); hold on
% plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'color',[0.7 0.7 0.7])
% axis square; tp_editplots
% Rmax(2,1) = pars(1);
% C50(2,1) = pars(2);
% 
% x = fr(start_val:end,ampa(2),inhib(2)-1);
% pars = tp_fitnakarushton(x,y,pars);
% plot(y,x,'r.','markersize',marker); hold on
% plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'r')
% axis square ; tp_editplots
% Rmax(2,2) = pars(1);
% C50(2,2) = pars(2);
% 
% x = fr(start_val:end,ampa(2),inhib(2)+1);
% pars = tp_fitnakarushton(x,y,pars);
% plot(y,x,'b.','markersize',marker); hold on
% plot(y,pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4),'b')
% axis square ; tp_editplots
% Rmax(2,3) = pars(1);
% C50(2,3) = pars(2);
% 
% xlabel('Contrast'); %ylabel('Firing rate [Hz]'); tp_editplots
% axis([-0.1 1 0 15])
% set(gca,'xtick',[0 0.5 1],'xticklabel',num2cell([0 0.5 1]))


%% BARPLOTS 
% Plot Rmax and C50

figure;  set(gcf,'color','w'); hold on
for i = 1: 2
  
  subplot(4,3,i); hold on
  
  bar(1,Rmax(i,3),'edgecolor','w','facecolor',[0 0 0])
  bar(2,Rmax(i,1),'edgecolor','k','facecolor',[1 1 1])
  bar(3,Rmax(i,2),'edgecolor','w','facecolor',[.5 .5 .5])
  tp_editplots
  
  
  axis([0.5 3.5 min(Rmax(i,:))-0.7*min(Rmax(i,:)) max(Rmax(i,:))+0.2*min(Rmax(i,:))])
  ylabel('R_{max}'); axis square
  
  tp_editplots
  
end
  
for i = 1: 2
  
  subplot(4,3,i+3); hold on
  bar(1,C50(i,3),'edgecolor','w','facecolor',[0 0 0])
  bar(2,C50(i,1),'edgecolor','k','facecolor',[1 1 1])
  bar(3,C50(i,2),'edgecolor','w','facecolor',[.5 .5 .5])
  tp_editplots
  
  axis([0.5 3.5 min(C50(i,:))-0.7*min(C50(i,:)) max(C50(i,:))+0.2*min(C50(i,:))])
  ylabel('C_{50}'); axis square
  
    tp_editplots

end
  
print(gcf,'-dpdf',sprintf('~/spiking/plots/spiking_ei_gain_barplots_C50_v2.pdf'))

  
%% MODEL FREE ESTIMATION OF I/O SLOPES
% Similar to Murphy & Miller (2003) J Neurosci
clear val

val(:,1)=fr([1 end],ampa(1),inhib(1)+1);
val(:,2)=fr([1 end],ampa(1),inhib(1));
val(:,3)=fr([1 end],ampa(1),inhib(1)-1);

delta_slope_incrEI(1) = 100*(val(2,3)-val(1,3)/size(fr,1)-val(2,2)-val(1,2)/size(fr,1))/val(2,2)-val(1,2)/size(fr,1)
delta_slope_decrEI(1) = 100*(val(2,1)-val(1,1)/size(fr,1)-val(2,2)-val(1,2)/size(fr,1))/val(2,2)-val(1,2)/size(fr,1)

val(:,1)=fr([1 end],ampa(2),inhib(2)+1);
val(:,2)=fr([1 end],ampa(2),inhib(2));
val(:,3)=fr([1 end],ampa(2),inhib(2)-1);

delta_slope_incrEI(2) = 100*(val(2,3)-val(1,3)/size(fr,1)-val(2,2)-val(1,2)/size(fr,1))/val(2,2)-val(1,2)/size(fr,1)
delta_slope_decrEI(2) = 100*(val(2,1)-val(1,1)/size(fr,1)-val(2,2)-val(1,2)/size(fr,1))/val(2,2)-val(1,2)/size(fr,1)


