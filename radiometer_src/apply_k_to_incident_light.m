%load in pre-processed k info:

load /Users/kristenhunter-cevera/MVCO_light_at_depth/radiometer_src/k_lite.mat


%% Interpolate avg k values for each day of year

%not for stations 7 and 8:
%wrap year around:
x=[kaggregate(end,1)-365; kaggregate(:,1); kaggregate(1,1)+365];
y=[-kaggregate(end,2); -kaggregate(:,2); -kaggregate(1,2)];
YY=interp1(x,y,1:366);

x=[kaggregate_hilow(end,1)-365; kaggregate_hilow(:,1); kaggregate_hilow(1,1)+365];
y=[-kaggregate_hilow(end,2); -kaggregate_hilow(:,2); -kaggregate_hilow(1,2)];
YY_low=interp1(x,y,1:366);

x=[kaggregate_hilow(end,1)-365; kaggregate_hilow(:,1); kaggregate_hilow(1,1)+365];
y=[-kaggregate_hilow(end,3); -kaggregate_hilow(:,3); -kaggregate_hilow(1,3)];
YY_high=interp1(x,y,1:366);


%% Link to and comparison with syn stuff!

load('/Users/kristenhunter-cevera/MVCO_light_at_depth/syn_data_analysis/mvco_envdata_15Dec2017.mat')
load('/Users/kristenhunter-cevera/MVCO_light_at_depth/syn_data_analysis/syndata_04Jan2017.mat')

%%
figure, plot(1:366, light_avg,'.-'), hold on
E_d = light_avg.*exp(4*-YY');  %DOUBLE CHECK THAT THIS IS INDEED THE RIGHT WAY TO CALCULATE THIS!!!
plot(1:366,E_d,'.-')

E_d = light_avg.*exp(4*-YY_low');  %DOUBLE CHECK THAT THIS IS INDEED THE RIGHT WAY TO CALCULATE THIS!!!
plot(1:366,E_d,'.-')

E_d = light_avg.*exp(4*-YY_high');  %DOUBLE CHECK THAT THIS IS INDEED THE RIGHT WAY TO CALCULATE THIS!!!
plot(1:366,E_d,'.-')


xlim([1 366])
ylabel('Daily average radiation (MJ m{-2})')
xlabel('Year day')
set(gca,'fontsize',14)
%%
figure
plot(Tcorr_avg, mu_avg,'o')

%%
clf
% subplot(1,2,1,'replace')
% scatter(light_avg, mu_avg,30,ydmu,'filled')
% subplot(1,2,2,'replace')
% scatter(E_d, mu_avg,30,ydmu,'filled')

subplot(1,2,1,'replace')
scatter(light_avg, mu_avg,30,ydmu,'filled')
ylabel('Division rate (d^{-1})')
xlabel('Average radiation (MJ m{-2})')
caxis([1 366])
colormap jet
set(gca,'fontsize',14,'box','on')
title('Incident')

subplot(1,2,2,'replace')
scatter(E_d, mu_avg,30,ydmu,'filled')
xlabel('Average radiation at 4m depth (MJ m{-2})')
caxis([1 366])
hbar=colorbar; set(hbar,'Ydir','reverse'); ylabel(hbar,'Year day')
set(gca,'fontsize',14,'box','on')
title('At depth')
%%
clf
% subplot(1,2,1,'replace')
% scatter(light_avg, PE_avg,30,ydmu,'filled')
% subplot(1,2,2,'replace')
% scatter(E_d, PE_avg,30,ydmu,'filled')

subplot(1,2,1,'replace')
scatter(light_avg, PE_avg,30,Tcorr_avg,'filled')
xlabel('Average radiation (MJ m{-2})')
ylabel('PE fluorescence')
colormap jet
caxis([-2 22])
set(gca,'fontsize',14,'box','on')
xlim([0 30])

subplot(1,2,2,'replace')
scatter(E_d, PE_avg,30,Tcorr_avg,'filled')
caxis([-2 22])
xlim([0 10])
xlabel('Average radiation at 4m depth (MJ m{-2})')

set(gca,'fontsize',14,'box','on')
hbar=colorbar; ylabel(hbar,'Temperature')


%% size seemed to show relationship with light rather than division rate...
%thought about showing size with division rate...but light seems to be
%better!

figure, scatter(E_d,SSC_avg,30,ydmu,'filled')
colormap jet
figure, scatter(mu_avg,SSC_avg,30,ydmu,'filled')