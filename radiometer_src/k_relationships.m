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
sourcepath=fullfile('/Volumes/Lab_data/MVCO/HyperPro_Radiometer/processed_radiometer_files/');

%How many of these do we have?
d = dir(sourcepath);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'\d{1,2}\w{3,4}20\d{2}'); %find all folders with this naming scheme
datafolders = foldernames(cellfun('isempty',temp)==0);

%%
cast_record=[];

for foldernum=[2:3 5:10 12 14:16 18:20];

    %find dat data!
    matsource=fullfile(sourcepath,datafolders{foldernum},'/mat_outfiles/');
    
    eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
    eval(['tempdata=data_' datafolders{foldernum} ';'])    
    
    eval(['load ' matsource 'K_PAR_' datafolders{foldernum} '.mat'])
    eval(['K_PAR=K_PAR_' datafolders{foldernum} ';'])

    eval(['load ' matsource 'location_' datafolders{foldernum} '.mat'])
    eval(['location=location_' datafolders{foldernum} ';'])
    
    %timestamp={tempdata.timestamp}'; timestamp=cell2mat(timestamp);
    %INCORRECT DATES LISTED HERE!!!
    
    timestamp=repmat(datenum(datestr(datafolders{foldernum})),1,length(tempdata))';
    if foldernum ~= 18
        lat={location.lat}'; 
        lon={location.lon}';
        station={location.station}';  
    else
        lat=cell(length(timestamp),1);
        lon=lat;
        station=lat;
    end
     K={K_PAR.K}'; 
    
    %replace empty cells with nan
%     qq=find(cellfun('isempty',K)==1);
%     if ~isempty(qq)
%         K(qq)={NaN}; 
%     end
    
    K((cellfun('isempty',K)==1))={NaN};
    lat((cellfun('isempty',lat)==1))={NaN}; %replace empty cells with NaN
    lon((cellfun('isempty',lon)==1))={NaN};
    station((cellfun('isempty',station)==1))={NaN};
        
    K=cell2mat(K);
    lat=cell2mat(lat);
    lon=cell2mat(lon);
    station=cell2mat(station);
    
    cast_record=[cast_record; timestamp K lat lon station];

end

[~, is]=sort(cast_record(:,1));
cast_record=cast_record(is,:);


%% a bit more complicated graphs based on position:

tower=[41.325 -70.5667];
node=[41.3366 -70.5564];
shore=[41.3499,-70.5267];

[year_day]=find_yearday(cast_record(:,1));

%let's grid all the points based on position:
box1=[-70.58 -70.55 41.135 41.155]; 
box2=[-70.58 -70.55 41.19 41.21];
box3=[-70.59 -70.52 41.23 41.26];
box_tower=[-70.58 -70.53 41.315 41.33];
box_node=[-70.58 -70.53 41.33 41.345];
box6=[-70.7 -70.6 41.315 41.33];
box7=[-70.52 -70.47 41.315 41.33];
box8=[-70.46 -70.42 41.315 41.33];
box9=[-69.85 -69.75 41.315 41.33];
box10=[-69.4 -69.3 41.315 41.33];

boxlist=whos('box*');

%% set up the indexes:
for j=1:length(boxlist)
    eval(['temp=' boxlist(j).name ';'])
    temp_ind=find(cast_record(:,3) > temp(3) & cast_record(:,3) < temp(4) & cast_record(:,4) > temp(1) & cast_record(:,4) < temp(2));
    varname=regexp(boxlist(j).name,'(?<=box)\w*','match'); varname=char(varname);
    eval(['b' varname '=temp_ind;'])
end

%%
subplot(2,3,1,'replace'), hold on

% scatter(cast_record(:,4),cast_record(:,3),40,year_day,'filled')

plot(cast_record(b1,4),cast_record(b1,3),'k.','markersize',12)
plot(cast_record(b2,4),cast_record(b2,3),'k.','markersize',12)
plot(cast_record(b3,4),cast_record(b3,3),'k.','markersize',12)
plot(cast_record(b9,4),cast_record(b9,3),'r.','markersize',12)
plot(cast_record(b10,4),cast_record(b10,3),'r.','markersize',12)

hold on

% box close to tower:
box=[-70.7 -70.4 41.3 41.35]; %left right bottom top bounds
ii=find(cast_record(:,3) > box(3) & cast_record(:,3) < box(4) & cast_record(:,4) > box(1) & cast_record(:,4) < box(2));
patch([-70.4 -70.4 -70.7 -70.7],[41.3 41.35 41.35 41.3],'k','facecolor','none')
plot(cast_record(ii,4),cast_record(ii,3),'.','markersize',12,'color',[0 0.5 1])
%legend('Casts','Tower','Node','Shore','location','southeast')
set(gca,'box','on','fontsize',14)
xlabel('Longitude')
ylabel('Latitude')
title('Locations of available radiomater casts')
%% draw boxes for each indivual station:

%draw a box for each:
for j=1:length(boxlist)
    eval(['temp=' boxlist(j).name ';'])
    x=[temp(1) temp(2) temp(2) temp(1)];
    y=[temp(3) temp(3) temp(4) temp(4)];
    patch(x,y,'k','facecolor','none')
end

%%
subplot(2,3,2,'replace'), hold on

o1=plot(year_day,-cast_record(:,2),'o','color',[0.5 0.5 0.5],'markersize',4);

cc=jet(10);
% plot(year_day(b1),-cast_record(b1,2),'.','markersize',14,'color',cc(1,:))
% plot(year_day(b2),-cast_record(b2,2),'.','markersize',14,'color',cc(2,:))
% plot(year_day(b3),-cast_record(b3,2),'.','markersize',14,'color',cc(3,:))
% plot(year_day(b_tower),-cast_record(b_tower,2),'s','markersize',10,'color',cc(4,:))
% plot(year_day(b_node),-cast_record(b_node,2),'p','markersize',10,'color',cc(5,:))
% plot(year_day(b6),-cast_record(b6,2),'.','markersize',14,'color',cc(6,:))
% plot(year_day(b7),-cast_record(b7,2),'.','markersize',14,'color',cc(7,:))
% plot(year_day(b8),-cast_record(b8,2),'.','markersize',14,'color',cc(8,:))
% plot(year_day(b9),-cast_record(b9,2),'.','markersize',14,'color',cc(9,:))
% plot(year_day(b10),-cast_record(b10,2),'.','markersize',14,'color',cc(10,:))


o3=plot(year_day(b1),-cast_record(b1,2),'.','markersize',14,'color','k');
plot(year_day(b2),-cast_record(b2,2),'.','markersize',14,'color','k')
plot(year_day(b3),-cast_record(b3,2),'.','markersize',14,'color','k')
o2=plot(year_day(ii),-cast_record(ii,2),'.','markersize',14,'color',[0 0.5 1]);
o4=plot(year_day(b9),-cast_record(b9,2),'.','markersize',14,'color','r');
plot(year_day(b10),-cast_record(b10,2),'.','markersize',14,'color','r')

set(gca,'fontsize',16,'box','on')
ylabel('Attenuation coefficient, K')
xlabel('Year Day')
xlim([1 366])
legend([o1; o2;o3;o4],'no lat/lon match yet','casts inside box, close to tower','casts south','casts east','location','NorthOutside')

%% Let's see how these relate to chlorophylls at MVCO:

load /Users/kristenhunter-cevera/Documents/MATLAB/MVCO_Syn_analysis/mvco_chlrep.mat

%%
subplot(2,3,4,'replace'), hold on

chl_yearday=find_yearday(matdate);
%plot(chl_yearday,chlavg(:,1),'.','color',[0 0.7 0])
set(gca,'fontsize',16,'box','on')
xlim([1 366])
xlabel('Year day')
ylabel('Chlorophyll Avg')
%Also a nice seasonal pattern!

%and for plotting, find all chl that came from box around sites close to
%tower:
chl_box=find(lat > box(3) & lat < box(4) & lon > box(1) & lon < box(2));
plot(chl_yearday(chl_box),chlavg(chl_box,1),'.','color',[0 0.7 0])
title('All chl values close to tower')

%% okay, so we would like to match chl samples to casts based on
% lat/lon/date...here we go!

chl_ind=cell(length(cast_record),1);
chl_match=nan(length(cast_record),4);
cast_days=unique(cast_record(:,1));
cc=0;

for j=1:length(cast_days)
    
    day=cast_days(j);
    qq=find(cast_record(:,1)==day); %matches to casts
    chl=find(floor(matdate)==day); %matches to chl
    
    if ~isempty(chl)
        %for each cast find, the closest chl in distance!
        for i=1:length(qq)
            cc=cc+1;
            
            lat_test=abs(cast_record(qq(i),3)-lat(chl));
            la = find(lat_test==min(lat_test)); %index of closest latitude
            
            lon_test=abs(cast_record(qq(i),4)-lon(chl));
            lo = find(lon_test==min(lon_test));  %index of closest lognitude
            
            chl2use=chl(intersect(lo,la)); %hopefully have some matching indexes!
            
            dd=find(depth(chl2use) <= 10  & depth(chl2use) > 2);                       
            
            chl_match(cc,1)=day;
            chl_match(cc,2)=nanmean(chlavg(chl2use,1)); %average of all measurements
            chl_match(cc,3)=nanmean(chlavg(chl2use(dd),1)); %those in between 2 and 10 m
            chl_match(cc,4)=nanmean(depth(chl2use(dd)));
            
            chl_ind(cc)={chl2use};           
        end
    
    else %[mm im]=min(abs(matdate-day));
        disp('cannot find this day?')
        keyboard
    end
        
end


%%
subplot(2,3,6,'replace')
hold on
plot(chl_match(:,2),-cast_record(:,2),'o','markersize',4,'color',[0.5 0.5 0.5])
%plot(chl_match(:,3),-cast_record(:,2),'ro')

%how about the ones just at the tower?
plot(chl_match(ii,2),-cast_record(ii,2),'.','markersize',14,'color',[0 0.5 1])
plot(chl_match(b1,2),-cast_record(b1,2),'.','markersize',14,'color','k')
plot(chl_match(b2,2),-cast_record(b2,2),'.','markersize',14,'color','k')
plot(chl_match(b3,2),-cast_record(b3,2),'.','markersize',14,'color','k')
plot(chl_match(b9,2),-cast_record(b9,2),'.','markersize',14,'color','r')
plot(chl_match(b10,2),-cast_record(b10,2),'.','markersize',14,'color','r')

set(gca,'box','on','fontsize',14)
xlabel('Chlorophyll Avg')
ylabel('Attenuation coefficient, K')

% a line through those tower ones?
% [bcoeffs,~,~,~,stats]=regress(-cast_record(ii,2),[ones(size(chl_match(ii,2))) chl_match(ii,2)]);
% x=sort(chl_match(ii,2));
% plot(x,bcoeffs(1)+bcoeffs(2)*x,'-','color',[0 0.5 1])
% text(0.5,0.47,['Significant linear fit, but with R2 of: ' num2str(1e-3*round(1000*stats(1)))])

%% Mean relaionships? Monthly climatologies?

addpath /Users/kristenhunter-cevera/Documents/MATLAB/mvco_tools/

%mean weekly chl:
[time_chl, daily_chl, chl_years] = timeseries2ydmat(matdate(chl_box), chlavg(chl_box,1));
[chl_wk_avg, chl_wk_std, yd_wk, chl_mn_avg, chl_mn_std, yd_mn] = dy2wkmn_climatology(daily_chl, chl_years);

% and add the average from the tower/node big box for k-values:

k_avg=[]; chl_avgB=[];
for j=1:length(cast_days)    
    jj=find(cast_record(ii,1)==cast_days(j));
    if ~isempty(jj)
        k_avg=[k_avg; unique(year_day(ii(jj))) nanmean(cast_record(ii(jj),2))];
        chl_avgB(j)=nanmean(chl_match(ii(jj),2));
    end
end

%
[~, is]=sort(k_avg(:,1));
k_avg=k_avg(is,:);
%
subplot(2,3,5,'replace'), hold on
[ax, h1, h2]=plotyy(yd_wk,chl_wk_avg,k_avg(:,1),-k_avg(:,2));
set(ax(1),'xlim',[1 366],'ycolor',[0 0.7 0],'fontsize',14)
set(ax(2),'xlim',[1 366],'ycolor',[0 0.5 1],'fontsize',14)
ylabel(ax(1),'Chlorophyll Avg')
ylabel(ax(2),'K')
set(h1,'color',[0 0.7 0],'marker','.','linewidth',2,'markersize',10)
set(h2,'color',[0 0.5 1],'marker','.','linewidth',2,'markersize',10)
xlabel('Year Day')

title('Weekly Chl avg, within day K avg')

%%
addpath /Users/kristenhunter-cevera/Documents/MATLAB/matlab_tools/export_fig_2016/
set(gcf,'color','w')
export_fig /Users/kristenhunter-cevera/Desktop/k-values.pdf

%% and a plot of K_par vs Chl with Morel 1988 equation and then fitted equation of that form:

%CAREFUL - XSCALE IS IN LOG SCLAE IN MOREL 1988!!!! Dah!!!!
clf
hold on
plot(chl_match(:,2),-cast_record(:,2),'o','markersize',4,'color',[0.5 0.5 0.5])
%plot(chl_match(:,3),-cast_record(:,2),'ro')

%how about the ones just at the tower?
plot(chl_match(ii,2),-cast_record(ii,2),'.','markersize',14,'color',[0 0.5 1])
plot(chl_match(b1,2),-cast_record(b1,2),'.','markersize',14,'color','k')
plot(chl_match(b2,2),-cast_record(b2,2),'.','markersize',14,'color','k')
plot(chl_match(b3,2),-cast_record(b3,2),'.','markersize',14,'color','k')
plot(chl_match(b9,2),-cast_record(b9,2),'.','markersize',14,'color','r')
plot(chl_match(b10,2),-cast_record(b10,2),'.','markersize',14,'color','r')

set(gca,'box','on','fontsize',14)
xlabel('Chlorophyll Avg')
ylabel('Attenuation coefficient, K')

% and the line from Morel 1988:
%K_PAR= 0.121 * C(mg/m3)^).428

plot(sort(chl_match(ii,2)),0.121*sort(chl_match(ii,2)).^0.428,'-','linewidth',2,'color',[0.5 0.5 0.5])


% [bcoeffs,~,~,~,stats]=regress(-cast_record(ii,2),[ones(size(chl_match(ii,2))) chl_match(ii,2)]);
% x=sort(chl_match(ii,2));

%%
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






%% wavelength plots

%find measured wavelengths:
%regexp('ED(348.76)','ED\((?<lambda>\d{3}\.\d{2})\)','names')
temp=regexp(edl_hdr(:,1),'\d{3}\.\d{2}','match');
ind=find(cellfun('isempty',temp)==0);
temp=temp(ind);
temp=[temp{:}]';
lambdas=cellfun(@(x) str2num(x),temp);

%extract just light measurments:
wv=cell2mat(edl(:,5:end-10));

