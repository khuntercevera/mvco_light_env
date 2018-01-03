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
lambdas=[];
k_lambdas=[];

for foldernum=good_data' %folders with viable casts in them
    
    %load in all the data!
    matsource=fullfile(processed_path,datafolders{foldernum},'/mat_outfiles/');
    
    eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
    eval(['tempdata=data_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'K_PAR_' datafolders{foldernum} '.mat'])
    eval(['K_PAR=K_PAR_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'location_' datafolders{foldernum} '.mat'])
    eval(['location=location_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'k_lambda_' datafolders{foldernum} '.mat'])
    eval(['k_lambda=k_lambda_' datafolders{foldernum} ';'])
    
    %because some of the casts were split, need to account for this:
    for filenum=1:length(K_PAR);
        
        if K_PAR(filenum).flag==0;
            
            k_values=[k_values; datenum(datestr(datafolders{foldernum})) location(filenum).lat location(filenum).lon K_PAR(filenum).K(2) K_PAR(filenum).stats(1) K_PAR(filenum).flag];
            k_record=[k_record; {datafolders{foldernum}} {foldernum} {K_PAR(filenum).file} {location(filenum).file} location(filenum).eventnum {K_PAR(filenum).NOTES} {location(filenum).notes}];
            k_lambdas=[k_lambdas; -k_lambda(filenum).k_wv(:,3)'];  %attenuation coefficient for each wavelength
            lambdas=[lambdas; k_lambda(filenum).k_wv(:,1)']; %wavelength used
            
        elseif K_PAR(filenum).flag==3;
            
            k_values=[k_values; datenum(datestr(datafolders{foldernum})) location(filenum).lat location(filenum).lon K_PAR(filenum).K1(2) NaN K_PAR(filenum).flag];
            k_record=[k_record; {datafolders{foldernum}} {foldernum} {K_PAR(filenum).file} {location(filenum).file} location(filenum).eventnum {K_PAR(filenum).NOTES} {location(filenum).notes}];
            k_lambdas=[k_lambdas; -k_lambda(filenum).k_wv1(:,3)']; %attenuation coefficient for each wavelength
            lambdas=[lambdas; k_lambda(filenum).k_wv1(:,1)']; %wavelength used
        end
    end
    
end

if sum(sum(lambdas + repmat(-lambdas(1,:),size(lambdas,1),1)))==0
    disp('all recorded wavelengths are the same!')
    lambdas=lambdas(1,:);
else
    disp('Ruh oh - not all wavelengths are the same?')
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

% set up the indexes:

for j=1:length(approx_station_loc)
    temp_ind=find(k_values(:,3) > approx_station_loc(j,3)-0.01 & k_values(:,3) < approx_station_loc(j,3)+0.01 ... %lat
        & k_values(:,4) > approx_station_loc(j,2)-0.01 & k_values(:,4) < approx_station_loc(j,2)+0.01); %lon
    k_values(temp_ind,8)=approx_station_loc(j,1);
end

%and color coordinated plot!
scatter(k_values(:,4),k_values(:,3),30,k_values(:,8),'filled')
set(gca,'box','on','fontsize',14)
xlabel('Longitude')
ylabel('Latitude')
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

%chlavg - those titles are:
%col 1: average of all replicates for whole chl
%col 2: average of all replicates for < 10 um chl
%col 3: average of all replicates for < 80 um chl
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
plot(chl_yearday(chl_tower),chlavg(chl_tower,1),'.','color',[0 0.6 0])
plot(chl_yearday(chl_node),chlavg(chl_node,1),'.','color',[0 0.3 0])
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


%% and a plot of K_par vs Chl with Morel 1988 equation and then fitted equation of that form:

%CAREFUL - XSCALE IS IN LOG SCLAE IN MOREL 1988!!!! Dah!!!!

% The curve from Morel 1988:
%K_PAR= 0.121 * C(mg/m3)^).428
plot(0.1:0.1:10,0.121*(0.1:0.1:10).^0.428,'-','linewidth',2,'color',[0.5 0.5 0.5])
title('note log scale for x-axis')
set(gca,'xscale','log')
xlim([0 10])

% fit a power curve through the data:
%all the data:
x = chl_match(:,2);
y = -k_values(:,5);
nn=find(~isnan(x) & ~isnan(y));

[X1,~,~,EXITFLAG] = lsqnonlin(@(theta) fit_powercurve(theta,x(nn),y(nn)),[0.2, 2.2]);
plot(sort(x),X1(1)*sort(x).^X1(2),'-','color',[0 0 0])

%just the tower and node:
tn=find(k_values(:,8)==3 | k_values(:,8)==4); %tower or node values
x = chl_match(tn,2);
y = -k_values(tn,5);
nn=find(~isnan(x) & ~isnan(y));
[X2,~,~,EXITFLAG] = lsqnonlin(@(theta) fit_powercurve(theta,x(nn),y(nn)),[0.2, 2.2]);
plot(sort(x),X2(1)*sort(x).^X2(2),'-','color',[0 0.5 1])

%this also works!
% f = fit(x(jj),y(jj),'power')
% figure,
% plot(f,x,y)

text(0.2,0.45,['Tower data: ' num2str(X1(1)) 'x^{' num2str(X1(2)) '}'],'color',[0 0.5 1])
text(0.2,0.42,['All data: ' num2str(X2(1)) 'x^{' num2str(X2(2)) '}'])
text(0.2,0.39,'Morel 1988: 0.121x^{0.428}','color',[0.4 0.4 0.4])



%% Climatologies...

addpath /Users/kristenhunter-cevera/Documents/MATLAB/mvco_tools/

% mean weekly chl:
[time_chl, daily_chl, chl_years] = timeseries2ydmat(matdate, chlavg(:,1));
[chl_wk_avg, chl_wk_std, yd_wk, chl_mn_avg, chl_mn_std, yd_mn] = dy2wkmn_climatology(daily_chl, chl_years);

[~,is]=sort(k_values(:,2));
[time_k, daily_k, k_years] = timeseries2ydmat(k_values(is,2), -k_values(is,5));
[k_wk_avg, k_wk_std, k_yd_wk, k_mn_avg, k_mn_std, k_yd_mn] = dy2wkmn_climatology(daily_k, k_years);

itn=union(chl_tower,chl_node); %indexes into chl for tower and node:
[time_chltn, daily_chltn, chltn_years] = timeseries2ydmat(matdate(itn), chlavg(itn,1));
[chltn_wk_avg, chltn_wk_std, yd_wk, chltn_mn_avg, chltn_mn_std, yd_mn] = dy2wkmn_climatology(daily_chltn, chltn_years);

[~,is]=sort(k_values(tn,2));
[time_ktn, daily_ktn, ktn_years] = timeseries2ydmat(k_values(tn(is),2), -k_values(tn(is),5));
[ktn_wk_avg, ktn_wk_std, ktn_yd_wk, ktn_mn_avg, ktn_mn_std, ktn_yd_mn] = dy2wkmn_climatology(daily_ktn, ktn_years);

%%
subplot(2,3,6,'replace')
nn=find(~isnan(k_wk_avg));
[ax, h1, h2]=plotyy(yd_wk,chl_wk_avg,k_yd_wk(nn),k_wk_avg(nn));
hold(ax(2))
hold(ax(1))
set(ax(1),'xlim',[1 366],'ycolor',[0 0.7 0],'fontsize',14)
set(ax(2),'xlim',[1 366],'ycolor',[0 0.5 1],'fontsize',14,'visible','on')
ylabel(ax(1),'Chlorophyll Avg')
ylabel(ax(2),'K')
set(h1,'color',[0 0.7 0],'marker','.','linewidth',2,'markersize',10)
set(h2,'color',[0 0.5 1],'marker','o','linestyle','none','markersize',8)
xlabel('Year Day')

h3=plot(ax(1),yd_wk,chltn_wk_avg,'.-','linewidth',2,'color',[0 0.3 0]);
h4=plot(ax(2),ktn_yd_wk,ktn_wk_avg,'o','markersize',8,'color',[0 0.3 0]);

legend([h1(1); h2(1); h3(1); h4(1)],'all chl','all k','tower/node chl','tower/node k')
title('Weekly averages')

%% save that figure!

addpath /Users/kristenhunter-cevera/Documents/MATLAB/matlab_tools/export_fig_2016/
set(gcf,'color','w')
export_fig /Users/kristenhunter-cevera/MVCO_light_at_depth/PAR_k_values.pdf

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%K(LAMBDAS)...comparison with Morel...

%So, each k at each wavelength was plotted against chl value and then a
%power function was fit through it (after subtracting kw)...

% import the values for comparison:
filename = '/Users/kristenhunter-cevera/MVCO_light_at_depth/radiometer_src/Morel_Maritorena_2001_k_fits.txt';
delimiter = '\t';
formatSpec = '%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);

Morel2001=cell2mat(dataArray(1:4));
titles_Morel2001={'wavelength';'Kw';'e param';'chi param'};

%initial plot for attenuation of water:
figure
plot(Morel2001(:,1), Morel2001(:,2),'.-')

%interpolate these to match lambdas used on our radiometer:
k_water=interp1(Morel2001(:,1), Morel2001(:,2),lambdas);

%% Ok, so let's see what those the Kbio looks like after has substracted Kw plotted against chlorophyll:
%for first pass, just compare directly wiht Morel2001, but should try to calculate kw
%directly from aw+1/2(bw)...

% [~, im]=min(abs(Morel2001(:,1)-lambdas(w))); %if don't want to interpolate and just use closest wavelength....
%k_bio=k_lambdas(:,w)-Morel2001(im,2);

k_bio=k_lambdas-repmat(k_water,92,1);

%do a linear regression of log-log transformed data to match with
%Morel:
rec=nan(length(lambdas),6);
for w=1:length(lambdas)
    
    x=[ones(size(chl_match(:,2))) log(chl_match(:,2))];
    y=log(k_bio(:,w)); %wavelength by wavelength
    [b,~,~,~,stats]=regress(y,x);
    
    [~, im]=min(abs(Morel2001(:,1)-lambdas(w))); %find closest to Morel
    
    if ~any(imag(b(:))) && ~any(isinf(b(:))) && ~any(isnan(b(:)))
        
        rec(w,:)=[lambdas(w) b(1) b(2) stats(1) stats(3) im]; %Morel2001(im,1:4)
        
        subplot(1,2,1,'replace'), hold on %linear view
        plot(log(chl_match(:,2)),log(k_bio(:,w)),'.')
        line([-2 10], b(1) + b(2)*[-2 10])
        title(['\lambda ' num2str(lambdas(w))])
        xlabel('Chlorophyll mg/m^{3}')
        ylabel('k_{bio} [k(\lambda) - k_w(\lambda)]')
        set(gca,'box','on','fontsize',14)
        
        %power function view:
        subplot(1,2,2,'replace'), hold on %linear view
        plot(chl_match(:,2),k_bio(:,w),'.')
        plot(test,exp(b(1)).*test.^b(2),'.-')
        plot(test,Morel2001(im,4).*(test.^Morel2001(im,3)),'.-')
        set(gca,'xscale','log','box','on','fontsize',14)
        set(gca,'yscale','log')
        xlim([0.01 100])
        ylim([0.001 10])
        title(['Morel match: ' num2str(Morel2001(im,1))])
        xlabel('Chlorophyll mg/m^{3}')
        ylabel('k_{bio} [k(\lambda) - k_w(\lambda)]')
        %because non negative y's can't appear on log scale
        
        set(gcf,'color','w')
        
        %pause
        
    end
end

rec_titles={'lambda' 'log chi' 'e' 'R2' 'p-value' 'index to closest lambda in Morel'};

%% a figure to show the consistent differences:

ll=[420 480 550 650]; %wavelengths to see
for i=1:4
    
    [~, w]=min(abs(lambdas-ll(i)));
    [~, im]=min(abs(Morel2001(:,1)-lambdas(w)));
    
    subplot(2,2,i,'replace'), hold on 
    plot(chl_match(:,2),k_bio(:,w),'.')
    plot(test,exp(rec(w,2)).*test.^rec(w,3),'-')
    plot(test,Morel2001(im,4).*(test.^Morel2001(im,3)),'-')
    set(gca,'xscale','log','box','on','fontsize',14)
    set(gca,'yscale','log')
    xlim([0.01 100])
    ylim([0.001 10])
    title(['\lambda: ' num2str(lambdas(w)) ' nm, Morel match: ' num2str(Morel2001(im,1)) ' nm'])
    xlabel('Chlorophyll mg/m^{3}')
    ylabel('k_{bio} [k(\lambda) - k_w(\lambda)]')
    %because non negative y's can't appear on log scale
end
set(gcf,'color','w')


%% hmmm...these plots look a bit different from Morel...

%figure 4 in Morel & Maritorena:
figure
subplot(1,2,1,'replace'), hold on
plot(rec(:,1),rec(:,3),'.-') %MVCO e
hold on
plot(Morel2001(:,1),Morel2001(:,3),'.-') %Morel e
title('e parameter')

subplot(1,2,2,'replace'), hold on
plot(rec(:,1),exp(rec(:,2)),'.-') %MVCO chi
plot(Morel2001(:,1),Morel2001(:,4),'.-') %Morel chi
title('Chi parameter')

%% let's construct their k_total vs wavelength plot and compare too:

%Figure 5 in Morel & Maritorena:

figure, hold on
k_MM = Morel2001(:,4).*(0.03).^Morel2001(:,3) + Morel2001(:,2);
plot(Morel2001(:,1),k_MM,'b.-')
k_MM = Morel2001(:,4).*(0.3).^Morel2001(:,3) + Morel2001(:,2);
plot(Morel2001(:,1),k_MM,'r.-')
k_MM = Morel2001(:,4).*(1).^Morel2001(:,3) + Morel2001(:,2);
plot(Morel2001(:,1),k_MM,'g.-')
k_MM = Morel2001(:,4).*(3).^Morel2001(:,3) + Morel2001(:,2);
plot(Morel2001(:,1),k_MM,'m.-')

k_MVCO = exp(rec(:,2)).*(0.03).^rec(:,3) + k_water';
plot(rec(:,1),k_MVCO,'b.:')
k_MVCO = exp(rec(:,2)).*(0.3).^rec(:,3) + k_water';
plot(rec(:,1),k_MVCO,'r.:')
k_MVCO = exp(rec(:,2)).*(1).^rec(:,3) + k_water';
plot(rec(:,1),k_MVCO,'g.:')
k_MVCO = exp(rec(:,2)).*(3).^rec(:,3) + k_water';
plot(rec(:,1),k_MVCO,'m.:')

xlim([400 700])



%% Interpolate avg k values for each day of year

%for just tower and node now:
%and the spline through those:
nn=find(~isnan(ktn_wk_avg));
%YY = spline(ktn_yd_wk(nn),ktn_wk_avg(nn),1:366);

%wrap year around:
x=ktn_yd_wk(nn); x=[x(end)-366; x(2:end); x(2)+366];
y=ktn_wk_avg(nn); y=[y(end); y(2:end); y(2)];
%for the moment, exclude low first value:
YY=interp1(x,y,1:366);
figure, hold on

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

plot(ktn_yd_wk(nn),ktn_wk_avg(nn),'o','markerfacecolor',[0.5 0.5 0.5])
plot(1:366,YY,'-','color',[0.5 0.5 0.5])








%% Link to and comparison with syn stuff!

load('/Users/kristenhunter-cevera/MVCO_light_at_depth/syn_data_analysis/mvco_envdata_15Dec2017.mat')
load('/Users/kristenhunter-cevera/MVCO_light_at_depth/syn_data_analysis/syndata_04Jan2017.mat')

%%
figure, plot(1:366, light_avg,'.-'), hold on
E_d = light_avg.*exp(4*-YY');
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
