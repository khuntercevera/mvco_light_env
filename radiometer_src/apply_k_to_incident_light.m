%load in pre-processed k info:

load /Users/kristenhunter-cevera/MVCO_light_at_depth/radiometer_src/k_lite.mat


%% 

%loaction plot:
subplot(1,2,1,'replace'), hold on
%plot station points:
plot(-70.567,41.325,'o','markersize',16,'color',[0.5 0.5 0.5]) %tower
plot(-70.555,41.335,'o','markersize',16,'color',[0.5 0.5 0.5]) %node
plot(-70.6275,41.325,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.505,41.325,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.45,41.3275,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.567,41.255,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.567,41.2,'o','markersize',16,'color',[0.5 0.5 0.5]) %station 7
plot(-70.567,41.145,'o','markersize',16,'color',[0.5 0.5 0.5]) %station 8

scatter(k_values(:,4),k_values(:,3),30,k_values(:,8),'filled')
set(gca,'box','on','fontsize',14)
xlabel('Longitude')
ylabel('Latitude')
%excellent - just the one outlier!


%% So, i think average all of the stations, except for the most outer two...and three?

%acutally, average casts that were less than 10 yeardays apart:
unq_yrdy=unique(k_avg(:,1));
kaggregate=nan(8,4);
kaggregate_hilow=nan(8,4);
count=0;

for q=1:length(unq_yrdy)
    
    if ismember(unq_yrdy(q),[66 72 73 82 255 267])
        switch unq_yrdy(q)
            case 66
                count=count+1;
                qq=find(k_avg(:,1) >= 66 & k_avg(:,1) <= 82 & k_avg(:,8)~=8 & k_avg(:,8)~=7 & k_avg(:,8)~=0 & k_avg(:,8)~=1); %exclude outer shelf casts
                kaggregate(count,1)=mean([66 72 73 82]);
                kaggregate_hilow(count,1)=mean([66 72 73 82]);
            case 267 %255 only has casts at 7 and 8...
                count=count+1;
                qq=find(k_avg(:,1) >= 255 & k_avg(:,1) <= 267 & k_avg(:,8)~=8 & k_avg(:,8)~=7 & k_avg(:,8)~=0 & k_avg(:,8)~=1); %exclude outer shelf casts
                kaggregate(count,1)=267;
                kaggregate_hilow(count,1)=267;
        end
        
        kaggregate(count,2)=nanmean(k_avg(qq,5));
        kaggregate(count,3)=nanstd(k_avg(qq,5));
        kaggregate(count,4)=length(qq); 
        
        kaggregate_hilow(count,2)=min(k_low(qq,5));
        kaggregate_hilow(count,3)=max(k_high(qq,5));
        kaggregate_hilow(count,4)=length(qq); 
        
    else
    count=count+1;    
    qq=find(k_avg(:,1)==unq_yrdy(q) & k_avg(:,8)~=8 & k_avg(:,8)~=7 & k_avg(:,8)~=0 & k_avg(:,8)~=1); %exclude outer shelf casts
    kaggregate(count,1)=unq_yrdy(q);
    kaggregate(count,2)=nanmean(k_avg(qq,5));
    kaggregate(count,3)=nanstd(k_avg(qq,5));
    kaggregate(count,4)=length(qq); 
    
    kaggregate_hilow(count,1)=unq_yrdy(q);
    kaggregate_hilow(count,2)=min(k_low(qq,5));
    kaggregate_hilow(count,3)=max(k_high(qq,5));
    kaggregate_hilow(count,4)=length(qq); 
    
    end

end
%average by month:
% yrdy_list=nan(12,1);
% for mn=1:12
%     yrdy_list(mn)=find_yearday(datenum([num2str(mn) '-1-03']));
% end
% yrdy_list=[yrdy_list; 367];
% 
% %% use k_avg as this is the average within an event:
% 
% mn_avg=nan(12,3);
% for q=1:12
%     qq=find(k_avg(:,1) >= yrdy_list(q) & k_avg(:,1) < yrdy_list(q+1) & k_avg(:,8)~=8 & k_avg(:,8)~=7 & k_avg(:,8)~=0 & k_avg(:,8)~=1); %exclude outer shelf casts
%     mn_avg(q,1)=nanmean(k_avg(qq,5));
%     mn_avg(q,2)=nanstd(k_avg(qq,5));
%     mn_avg(q,3)=length(qq);
% end

%% and plot :) 

subplot(1,2,2,'replace'), hold on


% tt=find(k_avg(:,8)~=7 & k_avg(:,8)~=8 & k_avg(:,8)~=0 & k_avg(:,8)~=1);
% scatter(k_avg(tt,1),-k_avg(tt,5),30,k_values(tt,8),'filled')

%plot the averages: with connections:
%calc slope...
% m1=kaggregate(end,2);
% m2=kaggregate(1,2);
% m=(m2-m1)/(365-347 + 16); 
% x=[1; kaggregate(:,1); 365];
% y=[-(m*(366-347)+m1); -kaggregate(:,2); -(m*(365-347)+m1)];

%or just pad:
x=[kaggregate(end,1)-365; kaggregate(:,1); kaggregate(1,1)+365];
y=[-kaggregate(end,2); -kaggregate(:,2); -kaggregate(1,2)];
plot(x, y,':','color',[0.3 0.3 0.3],'linewidth',2)

x=[kaggregate_hilow(end,1)-365; kaggregate_hilow(:,1); kaggregate_hilow(1,1)+365];
y=[-kaggregate_hilow(end,2); -kaggregate_hilow(:,2); -kaggregate_hilow(1,2)];
plot(x, y,':','color',[0.6 0.6 0.6],'linewidth',2)

x=[kaggregate_hilow(end,1)-365; kaggregate_hilow(:,1); kaggregate_hilow(1,1)+365];
y=[-kaggregate_hilow(end,3); -kaggregate_hilow(:,3); -kaggregate_hilow(1,3)];
plot(x, y,':','color',[0.6 0.6 0.6],'linewidth',2)

%plot(kaggregate(:,1), -kaggregate(:,2),'.','markersize',12,'color',[0.5 0.5 0.5])

tt=find(k_values(:,8)~=7 & k_values(:,8)~=8 & k_values(:,8)~=0 & k_values(:,8)~=1);
scatter(k_values(tt,1),-k_values(tt,5),30,k_values(tt,8),'filled')

colorbar
caxis([1 8])
set(gca,'fontsize',14,'box','on')
ylabel('Attenuation coefficient, K')
xlabel('Year Day')
xlim([1 366])
title('k for all stations')

colormap jet
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