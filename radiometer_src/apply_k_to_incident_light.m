%load in pre-processed k info:

load '/Users/kristenhunter-cevera/mvco_light_env/radiometer_src/k_lite.mat'

%% Interpolate avg k values for each day of year

%not for stations 7 and 8:
%wrap year around:
x=[kaggregate(end,1)-365; kaggregate(:,1); kaggregate(1,1)+365];
y=[-kaggregate(end,2); -kaggregate(:,2); -kaggregate(1,2)];
k_interp=interp1(x,y,1:366);

x=[kaggregate_hilow(end,1)-365; kaggregate_hilow(:,1); kaggregate_hilow(1,1)+365];
y=[-kaggregate_hilow(end,2); -kaggregate_hilow(:,2); -kaggregate_hilow(1,2)];
k_interp_low=interp1(x,y,1:366);

x=[kaggregate_hilow(end,1)-365; kaggregate_hilow(:,1); kaggregate_hilow(1,1)+365];
y=[-kaggregate_hilow(end,3); -kaggregate_hilow(:,3); -kaggregate_hilow(1,3)];
k_interp_high=interp1(x,y,1:366);


%% load in light

load('/Volumes/Lab_data/MVCO/FCB/Syn_and_MVCO_packaged_data/mvco_envdata_current.mat','light_avg')

%% LIGHT AT 4m depth:

%just to see
% figure, plot(1:366, light_avg,'.-'), hold on
% E_d = light_avg.*exp(4*-k_interp');  
% plot(1:366,E_d,'.-')
% 
% E_d = light_avg.*exp(4*-YY_k_interp_low');  
% plot(1:366,E_d,'.-')
% 
% E_d = light_avg.*exp(4*-YY_k_interp_high');  
% plot(1:366,E_d,'.-')
% 
% xlim([1 366])
% ylabel('Daily average radiation (MJ m{-2})')
% xlabel('Year day')
% set(gca,'fontsize',14)


%% Average light over all the depths:

%numerically:
stepsize=0.2;
depths=0:stepsize:15; 

for j=1:366
    light_depth_n(j) = light_avg(j)./15 * stepsize * trapz(exp(depths.*-k_interp(j)'));  
    light_depth_low_n(j) = light_avg(j)./15 * stepsize * trapz(exp(depths.*-k_interp_low(j)')); 
    light_depth_high_n(j) = light_avg(j)./15 * stepsize * trapz(exp(depths.*-k_interp_high(j)'));  
end

% or analytically:
for j=1:366
    light_depth(j) = (light_avg(j)./(15*k_interp(j))) * (1-exp(15*-k_interp(j)));  
    light_depth_low(j) = (light_avg(j)./(15*k_interp_low(j))) * (1-exp(15*-k_interp_low(j))); 
    light_depth_high(j) = (light_avg(j)./(15*k_interp_high(j))) * (1-exp(15*-k_interp_high(j)));  
end


%% to see :)
figure, plot(1:366, light_avg,'.-'), hold on
plot(1:366,light_depth,'.-')
plot(1:366,light_depth_low,'.-')
plot(1:366,light_depth_high,'.-')

%%
save /Users/kristenhunter-cevera/mvco_light_env/radiometer_src/k_interp.mat k_interp* light_depth light_depth_low light_depth_high

