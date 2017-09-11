%% Relationships with K at MVCO and other shelf areas

%Ideally, we'd like a relationship between K_d and chl that we can use to
%extrapolate for the rest of the year to get a rough estimate of what in
%situ light levels were like over the time series...

%The below script goes through each separate radiometer processed file and
%pulls out k-values for easier plotting

% some plots of k
%imports chl data and plots this against k-values...

%% A time series plot of K_d values!

%sourcepath=fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\processed_radiometer_files\'); % path to folders with raw data...
sourcepath=fullfile('/Volumes/Lab_data/MVCO/HyperPro_Radiometer/');
processed_path=fullfile(sourcepath,'/processed_radiometer_files/');

%How many of these do we have?
d = dir(processed_path);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'\d{1,2}\w{3,4}20\d{2}'); %find all folders with this naming scheme
datafolders = foldernames(cellfun('isempty',temp)==0);

load(fullfile(sourcepath,'good_data_folders.mat'))

%% Gather K values from all the data to plot by time and by lat & lon...

k_values=[];
k_record={};

for foldernum=good_data' %folders with viable casts in them
    
    %load in all the data!
    matsource=fullfile(processed_path,datafolders{foldernum},'/mat_outfiles/');
    
    eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
    eval(['tempdata=data_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'K_PAR_' datafolders{foldernum} '.mat'])
    eval(['K_PAR=K_PAR_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'location_' datafolders{foldernum} '.mat'])
    eval(['location=location_' datafolders{foldernum} ';'])
    
    %because some of the casts were split, need to account for this:
    for filenum=1:length(K_PAR);
        
        if K_PAR(filenum).flag==0;
            
            k_values=[k_values; datenum(datestr(datafolders{foldernum})) location(filenum).lat location(filenum).lon K_PAR(filenum).K(2) K_PAR(filenum).stats(1) K_PAR(filenum).flag];
            k_record=[k_record; {datafolders{foldernum}} {foldernum} {K_PAR(filenum).file} {location(filenum).file} location(filenum).eventnum {K_PAR(filenum).NOTES} {location(filenum).notes}];
            
        elseif K_PAR(filenum).flag==3;
            
            k_values=[k_values; datenum(datestr(datafolders{foldernum})) location(filenum).lat location(filenum).lon K_PAR(filenum).K1(2) NaN K_PAR(filenum).flag];
            k_record=[k_record; {datafolders{foldernum}} {foldernum} {K_PAR(filenum).file} {location(filenum).file} location(filenum).eventnum {K_PAR(filenum).NOTES} {location(filenum).notes}];
        end
    end
    
end

%add in yearday:
k_values=[find_yearday(k_values(:,1)) k_values];

%% organize these points based on position:

%a quick plot of where these points lie:
subplot(2,3,1,'replace'), hold on
%plot station points:
plot(-70.567,41.325,'o','markersize',16,'color',[0.5 0.5 0.5]) %tower
plot(-70.555,41.335,'o','markersize',16,'color',[0.5 0.5 0.5]) %node
plot(-70.6275,41.325,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.505,41.325,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.45,41.3275,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.567,41.255,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.567,41.2,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.567,41.145,'o','markersize',16,'color',[0.5 0.5 0.5]) %station

plot(k_values(:,4),k_values(:,3),'r.')

%% Now code each cast more or less to a station number:

approx_station_loc=[4 -70.567 41.325; %tower
    3 -70.5564 41.3366; %node
    1 -70.45 41.3275; %station 1
    2 -70.505 41.325; %station 2
    5 -70.6275 41.325;  %station 5
    6 -70.567 41.255; %station 6
    7 -70.567 41.2; %station 7
    8 -70.567 41.145]; %station 8

shore_pos=[41.3499,-70.5267];


%% set up the indexes:

for j=1:length(approx_station_loc)
    temp_ind=find(k_values(:,3) > approx_station_loc(j,3)-0.01 & k_values(:,3) < approx_station_loc(j,3)+0.01 ... %lat
        & k_values(:,4) > approx_station_loc(j,2)-0.01 & k_values(:,4) < approx_station_loc(j,2)+0.01); %lon
    k_values(temp_ind,8)=approx_station_loc(j,1);
end

%and color coordinated plot!
scatter(k_values(:,4),k_values(:,3),30,k_values(:,8),'filled')
set(gca,'box','on','fontsize',14)
xlabel('Longitude')
ylabel('latitude')
%excellent - just the one outlier!

%% now examine k, by yearday, location, etc...

subplot(2,3,2,'replace'), hold on

scatter(k_values(:,1),-k_values(:,5),30,k_values(:,8),'filled')
colorbar
set(gca,'fontsize',14,'box','on')
ylabel('Attenuation coefficient, K')
xlabel('Year Day')
xlim([1 366])
title('k for all stations')
% unqdays=unique(k_values(:,1));
% for q=1:length(unqdays)
%     line([unqdays(q) unqdays(q)],ylim,'color',[0.6 0.6 0.6])
% end

subplot(2,3,3,'replace'), hold on
%light gray vertical lines for each sampling day:
xlim([1 366]); ylim([0.15 0.5])
unqdays=unique(k_values(:,1));
for q=1:length(unqdays)
    line([unqdays(q) unqdays(q)],ylim,'color',[0.6 0.6 0.6])
end
ii=find(k_values(:,8)==3 | k_values(:,8)==4);
scatter(k_values(ii,1),-k_values(ii,5),30,k_values(ii,8),'filled')
caxis([0 8]), colorbar
set(gca,'fontsize',14,'box','on')
ylabel('Attenuation coefficient, K')
xlabel('Year Day')
xlim([1 366])
title('k for just tower and node locations')



%% RELATIONSHIP TO CHLOROPHYLLS at MVCO:

load /Users/kristenhunter-cevera/Documents/MATLAB/MVCO_Syn_analysis/mvco_chlrep.mat

%%
subplot(2,3,4,'replace'), hold on

chl_yearday=find_yearday(matdate);
%plot(chl_yearday,chlavg(:,1),'.','color',[0 0.7 0])
set(gca,'fontsize',14,'box','on')
xlim([1 366])
xlabel('Year day')
ylabel('Chlorophyll Avg')
%Also a nice seasonal pattern!

%and for plotting, find all chl that came from box around sites close to
%tower:
chl_tower=find(lat > approx_station_loc(1,3)-0.01 & lat < approx_station_loc(1,3)+0.01 & ...
    lon > approx_station_loc(1,2)-0.01 & lon < approx_station_loc(1,2)+0.01);
chl_node=find(lat > approx_station_loc(2,3)-0.01 & lat < approx_station_loc(2,3)+0.01 & ...
    lon > approx_station_loc(2,2)-0.01 & lon < approx_station_loc(2,2)+0.01);
plot(chl_yearday(chl_tower),chlavg(chl_tower,1),'.','color',[0 0.7 0])
plot(chl_yearday(chl_node),chlavg(chl_node,1),'.','color',[0 0.5 0])
title('Avg total Chl values at tower and node')

%% RELATIONSHIP? MATCH CHL TO CASTS

%use event number to match...
chl_match=nan(length(k_record),2);
for q=1:length(k_record)
    
    ii=find(cellfun('isempty',regexp(cellstr(event),k_record{q,5}))==0);
    if ~isempty(ii)
        chl_match(q,1)=k_values(q,2);
        chl_match(q,2)=nanmean(chlavg(ii,1)); %for now, just use the average of all the measurements at depth
    else
        keyboard
    end
end

  
%%
subplot(2,3,5,'replace')
hold on
plot(chl_match(:,2),-k_values(:,5),'o','markersize',4,'color',[0.5 0.5 0.5])
%plot(chl_match(:,3),-cast_record(:,2),'ro')
ylim([0 0.6])
%how about the ones just at the tower?
jj=find(k_values(:,8)==3 | k_values(:,8)==4);
plot(chl_match(jj,2),-k_values(jj,5),'.','markersize',14,'color',[0 0.5 1])

set(gca,'box','on','fontsize',14)
xlabel('Chlorophyll Avg')
ylabel('Attenuation coefficient, K')

% a line through the points? Maybe doesn't make much sense...isn't an
% exponential relationship? See Morel...
% [bcoeffs,~,~,~,stats]=regress(-k_values(:,5),[ones(size(chl_match(:,2))) chl_match(:,2)]);
% x=sort(chl_match(:,2));
% plot(x,bcoeffs(1)+bcoeffs(2)*x,'-','color',[0.5 0.5 0.5])
% text(0.5,0.47,['Significant linear fit, but with R2 of: ' num2str(1e-3*round(1000*stats(1)))])
% [bcoeffs,~,~,~,stats]=regress(-k_values(jj,5),[ones(size(chl_match(jj,2))) chl_match(jj,2)]);
% x=sort(chl_match(jj,2));
% plot(x,bcoeffs(1)+bcoeffs(2)*x,'-','color',[0 0.5 1])
% text(0.5,0.45,['Significant linear fit, but with R2 of: ' num2str(1e-3*round(1000*stats(1)))])


%% and a plot of K_par vs Chl with Morel 1988 equation and then fitted equation of that form:

%CAREFUL - XSCALE IS IN LOG SCLAE IN MOREL 1988!!!! Dah!!!!

% The curve from Morel 1988:
%K_PAR= 0.121 * C(mg/m3)^).428
jj=find(k_values(:,8)==3 | k_values(:,8)==4);
plot(sort(chl_match(jj,2)),0.121*sort(chl_match(jj,2)).^0.428,'-','linewidth',2,'color',[0.5 0.5 0.5])

% [bcoeffs,~,~,~,stats]=regress(-cast_record(ii,2),[ones(size(chl_match(ii,2))) chl_match(ii,2)]);
% x=sort(chl_match(ii,2));

x = chl_match(ii,2);
y = -cast_record(ii,2);
jj=find(~isnan(x) & ~isnan(y));

[X1,~,~,EXITFLAG] = lsqnonlin(@(theta) fit_powercurve(theta,x(jj),y(jj)),[0.2, 2.2]);
plot(sort(x),X1(1)*sort(x).^X1(2),'-','color',[0 0.5 1])

x = chl_match(:,2);
y = -cast_record(:,2);
jj=find(~isnan(x) & ~isnan(y));
[X2,~,~,EXITFLAG] = lsqnonlin(@(theta) fit_powercurve(theta,x(jj),y(jj)),[0.2, 2.2]);
plot(sort(x),X2(1)*sort(x).^X2(2),'-','color',[0 0 0])

%this also works!
% f = fit(x(jj),y(jj),'power')
% figure,
% plot(f,x,y)

set(gca,'xscale','log')
text(0.5,0.45,['Tower data: ' num2str(X1(1)) '+x^{' num2str(X1(2)) '}'],'color',[0 0.5 1])
text(0.5,0.42,['All data: ' num2str(X2(1)) '+x^{' num2str(X2(2)) '}'])
text(0.5,0.40,'Morel 1988: 0.121+x^{0.428}','color',[0.4 0.4 0.4])






%% Mean relaionships? Monthly climatologies?

% addpath /Users/kristenhunter-cevera/Documents/MATLAB/mvco_tools/
% 
% %mean weekly chl:
% [time_chl, daily_chl, chl_years] = timeseries2ydmat(matdate(chl_box), chlavg(chl_box,1));
% [chl_wk_avg, chl_wk_std, yd_wk, chl_mn_avg, chl_mn_std, yd_mn] = dy2wkmn_climatology(daily_chl, chl_years);
% 
% % and add the average from the tower/node big box for k-values:
% 
% k_avg=[]; chl_avgB=[];
% for j=1:length(cast_days)
%     jj=find(cast_record(ii,1)==cast_days(j));
%     if ~isempty(jj)
%         k_avg=[k_avg; unique(year_day(ii(jj))) nanmean(cast_record(ii(jj),2))];
%         chl_avgB(j)=nanmean(chl_match(ii(jj),2));
%     end
% end
% 
% %
% [~, is]=sort(k_avg(:,1));
% k_avg=k_avg(is,:);
% %
% subplot(2,3,5,'replace'), hold on
% [ax, h1, h2]=plotyy(yd_wk,chl_wk_avg,k_avg(:,1),-k_avg(:,2));
% set(ax(1),'xlim',[1 366],'ycolor',[0 0.7 0],'fontsize',14)
% set(ax(2),'xlim',[1 366],'ycolor',[0 0.5 1],'fontsize',14)
% ylabel(ax(1),'Chlorophyll Avg')
% ylabel(ax(2),'K')
% set(h1,'color',[0 0.7 0],'marker','.','linewidth',2,'markersize',10)
% set(h2,'color',[0 0.5 1],'marker','.','linewidth',2,'markersize',10)
% xlabel('Year Day')
% 
% title('Weekly Chl avg, within day K avg')
% 
% %%
% addpath /Users/kristenhunter-cevera/Documents/MATLAB/matlab_tools/export_fig_2016/
% set(gcf,'color','w')
% export_fig /Users/kristenhunter-cevera/Desktop/k-values.pdf
% 