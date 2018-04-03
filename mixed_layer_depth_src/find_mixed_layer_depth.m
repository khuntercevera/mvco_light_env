%Now that we have the downcast and a rough cut on quality control,
%we can start to identify where the mixed layer depth might be...

%use of two metrics:
%   Brunt Vaisala frequency (average, median, or max)
%   Depth at whcich temperature deviates from 0.1 0r 0.X degrees from surface

%for calculating Brunt-Vaisala:
addpath /Users/kristenhunter-cevera/MVCO_light_at_depth/seawater_ver3_2/
addpath /Users/kristenhunter-cevera/Documents/MATLAB/mvco_tools/

plotflag=1;
warning off

%first, find all trips around MVCO:

%% find only trips to tower or node:

box_tower=[-70.58 -70.53 41.315 41.33];
box_node=[-70.58 -70.53 41.33 41.345];

mvco_ind=find(downcast_lat > 41.3 & downcast_lat < 41.35 & downcast_lon < -70.53 & downcast_lon > -70.60 & cellfun('isempty',{CTD_QC(:).flag})==1);
%cellfun('isempty',regexp({CTD_QC(:).cast_name}','(deck)|(test)'))==1)
%this should be about ~104 casts...

% and if curious, see where all the casts come from...
figure
plot(downcast_lon(mvco_ind),downcast_lat(mvco_ind),'.')
patch([-70.53 -70.53 -70.60 -70.60],[41.3 41.35 41.35 41.3],'k','facecolor','none')

%%
for q=1:length(mvco_ind);
    
    col_hdr=CTD_QC(mvco_ind(q)).data_hdr;
    temp_data=CTD_QC(mvco_ind(q)).data;
    
    pdens_deltas=[0.005 0.01 0.05 0.1];
    temp_deltas=[0.2 0.4 0.6 0.8 1];
    
    temp_ref=[];
    pdens_ref=[];
    
    rec_pdens=[pdens_deltas' nan(size(pdens_deltas))'];
    rec_temp=[temp_deltas' nan(size(temp_deltas))'];
    
    %could do this via indexing....
    j=1; %counter for depth bins
    h=1; %counter for values in pdens threshold
    g=1; %counter for values in temp threshold
    while j < length(downcast_bins)-1
        j=j+1;
        if abs(pdens_ref - downcast_pdens(j,q)) > pdens_deltas(h)
            rec_pdens(h,2)=downcast_bins(j,q);
            h=h+1;
        end
        if abs(temp_ref - downcast_temp(j,q)) > temp_deltas(g)
            rec_temp(h,2)=downcast_bins(j,q);
            g=g+1;
        end
    end
    
    %calculate N2 from binned data:
    N2=sw_bfrq(downcast_sal(:,q),downcast_temp(:,q),downcast_press(:,q),downcast_lat(q));
    
    clf %see what this metric is highlighting - over all looks pretty good!
    subplot(2,3,1,'replace'),  hold on
    %plot(pdens,depth,'k.-')
    plot(downcast_pdens(:,q),downcast_bins(:,q),'.-','color',[0 0.5 1])
    plot(pdensU,depthU,'.','color',[1 0.5 0])
    plot(binned_dataD(:,6),binned_dataD(:,1),'.-')
    line(xlim,[binned_dataD(im,1) binned_dataD(im,1)],'color','r')
    line(xlim,[4 4],'color',[0.5 0.5 0.5])
    set(gca,'ydir','reverse','fontsize',14)
    xlabel('Potential Density (kg/m^3)') %ylabel(col_hdr{6})
    title([CTD(mvco_ind(q)).cast_name ':' datestr(file_time) ' : ' num2str(q) ' out of ' num2str(length(mvco_ind))],'interpreter','none')
    %legend('Obs \rho','Potential \rho','location','NorthWest')
    
    subplot(2,3,2,'replace'), hold on
    plot(cast_time,depth,'k.-')
    plot(cast_timeD,depthD,'.','color',[0 0.5 1])
    plot(cast_timeU,depthU,'.','color',[1 0.5 0])
    set(gca,'ydir','reverse','fontsize',14)
    xlabel('Time') %ylabel(col_hdr{6})
    title('CTD position with time')
    
    subplot(2,3,3,'replace'), hold on
    plot(N2D,binned_dataD(1:end-1,1),'.-','color',[0 0.5 1])
    plot(N2U,binned_dataU(1:end-1,1),'r.-','color',[1 0.5 0])
    line([1e-4 1e-4], get(gca,'ylim'),'color','r') %anything higher than this and you are probably stratified
    line([avgN2 avgN2], get(gca,'ylim'),'color','c') %average N2
    
    plot(N2(im),binned_dataD(im,1),'bp','markersize',12) %max only considering below 3.75m
    set(gca,'ydir','reverse','fontsize',14)
    xlabel('Brunt-Vaisala Freq')
    xlim([-1e-3 1e-2])
    line(xlim,[4 4],'color',[0.5 0.5 0.5])
    title('N2 with depth')
    
    subplot(2,3,4,'replace'),  hold on
    %plot(temperature,depth,'k.-')
    plot(temperatureD,depthD,'.','color',[0 0.5 1])
    plot(binned_dataD(:,4),binned_dataD(:,1),'.-')
    plot(temperatureU,depthU,'.','color',[1 0.5 0])
    plot(binned_dataU(:,4),binned_dataU(:,1),'r.-')
    set(gca,'ydir','reverse','fontsize',14)
    xlabel('Temperature (\circC)') %ylabel(col_hdr{6})
    line(xlim,[4 4],'color',[0.5 0.5 0.5])
    title('Temperature with Depth')
    
    subplot(2,3,5,'replace'),  hold on
    plot(sal,depth,'k.-')
    plot(salD,depthD,'.','color',[0 0.5 1])
    plot(binned_dataD(:,3),binned_dataD(:,1),'.-')
    plot(salU,depthU,'.','color',[1 0.5 0])
    plot(binned_dataU(:,3),binned_dataU(:,1),'.-')
    set(gca,'ydir','reverse','fontsize',14)
    xlabel('Salinity') %ylabel(col_hdr{6})
    xlim([median(salD)-1 median(salD)+1])
    line(xlim,[4 4],'color',[0.5 0.5 0.5])
    title('Salinity with Depth')
    
    subplot(2,3,6,'replace'),  hold on
    plot(cast_time,sal,'k.-')
    plot(cast_timeD,salD,'.','color',[0 0.5 1])
    plot(cast_timeU,salU,'.','color',[1 0.5 0])
    plot(cast_time,depth,'.-')
    ylabel('Salinity') %ylabel(col_hdr{6})
    
    
    
    %Is the water stratified? Or wouldn't be well mixed?
    %Looking at Young-Oh's slide, it looks like the min would be around
    %.25 * 10^-4 N2 (S^{-2})...
    
    disp(['Max N2: ' num2str(N2(im))])
    disp(['Avg N2: ' num2str(avgN2)])
    disp(['Median N2: ' num2str(medN2)])
    
    keyboard
    
    %record info:
    ctd_mix(q).cast_name = CTD(mvco_ind(q)).cast_name;
    ctd_mix(q).empty_flag = 0;
    ctd_mix(q).file_time =file_time;
    
    %record key processing steps:
    ctd_mix(q).descent_ind=dsc;
    ctd_mix(q).upward_ind=usc;
    ctd_mix(q).binned_downdata=binned_dataD;
    ctd_mix(q).binned_updata=binned_dataU;
    ctd_mix(q).bfrq_down=N2D;
    ctd_mix(q).bfrq_up=N2U;
    ctd_mix(q).bin_titles=binned_data_titles;
    
    %record some metrics and decisions:
    ctd_mix(q).user_call= user_call;
    ctd_mix(q).max_N2 = N2(im);
    ctd_mix(q).depth_maxN2= binned_data(im,1);
    ctd_mix(q).avg_N2 = avgN2;
    ctd_mix(q).med_N2= medN2;
    ctd_mix(q).data_at_4m=[binned_data(ii4,1) N2(ii4) binned_data(ii4,4) binned_data(ii4,6)];
    ctd_mix(q).data_at_12m=[binned_data(ii12,1) N2(ii12) binned_data(ii12,4) binned_data(ii12,6)];
    ctd_mix(q).cast_used = upcast_flag;
    
    ctd_mix(q).notes=notes;
    
    if strcmp(user_call,'mixed')
        ctd_mix(q).stratification_driver='';
    else
        ctd_mix(q).stratification_driver=reason;
    end
    
    %         [Q] = input('Does the automatic finding of N2 look reasonable? If so, enter "auto", if not, enter "manual" \n');
    %         if any(diff(cast_time) < 0), disp('Something wrong with time sync?'), end
    
    
end %for loop


%% Okay and now some plots to try to make sense of this mess...

%Compare difference of temp to diff of density, color coded for stratified
%or not....

%first remove empty rows where something went awry:
jj=find(cell2mat({ctd_mix(:).empty_flag}')==0); %use only good data casts

mix=ctd_mix(jj);

%exclude points that are not close enough to 12m?

%ii=find(abs(cell2mat(mld(jj,10))-12) < 1);

%% label points based on mixed or stratified:

user_call={mix(:).user_call}';
mm0=find(strcmp(user_call,'mixed'));
ss0=find(strcmp(user_call,'stratified'));

%by max:
maxN2=cell2mat({mix(:).max_N2}');
ss1=find(maxN2 >= 2e-5);
mm1=find(maxN2 < 2e-5);

% by average:
avgN2=cell2mat({mix(:).avg_N2}');
ss2=find(avgN2 >= 2e-5);
mm2=find(avgN2 < 2e-5);

%by median:
medN2=cell2mat({mix(:).med_N2}');
ss3=find(medN2 >= 2e-5);
mm3=find(medN2 < 2e-5);

%choose an a method for indexing:
mm=mm0;
ss=ss0;

time=cell2mat({mix(:).file_time}');
yrdy=find_yearday(time);
depth_maxN2 = cell2mat({mix(:).depth_maxN2}');
data_at_4m=cell2mat({mix(:).data_at_4m}');
data_at_12m=cell2mat({mix(:).data_at_12m}');


%% diff in temp vs. diff in dens

clf
% subplot(2,4,1,'replace')
% plot(cell2mat(mldBuse(ss,2)),cell2mat(mldBuse(ss,3)),'r.','markersize',12), hold on
% plot(cell2mat(mldBuse(mm,2)),cell2mat(mldBuse(mm,3)),'b.','markersize',12)
% ylabel('Max N2 (s^{-2})')
% line(xlim,[3e-4 3e-4],'color',[0.5 0.5 0.5])
% datetick

%Max N2 by year day:
subplot(2,4,1,'replace')
plot(yrdy,maxN2,'.','markersize',12), hold on
% plot(yrdy(ss),maxN2(ss),'r.','markersize',12), hold on
% plot(yrdy(mm),maxN2(mm),'b.','markersize',12)
ylabel('Max N2 (s^{-2})')
xlim([0 365])
xlabel('Yearday')
line(xlim,[2e-5 2e-5],'color',[0.5 0.5 0.5])

%Average N2 by year day:
subplot(2,4,2,'replace')
plot(yrdy,avgN2,'.','markersize',12), hold on
% plot(yrdy(ss),avgN2(ss),'r.','markersize',12), hold on
% plot(yrdy(mm),avgN2(mm),'b.','markersize',12)
ylabel('Average N2 (s^{-2})')
xlim([0 365])
xlabel('Yearday')
line(xlim,[2e-5 2e-5],'color',[0.5 0.5 0.5])
line(xlim,[5e-5 5e-5],'color',[0.5 0.5 0.5])

%Depth of max N2:
subplot(2,4,3,'replace')
plot(maxN2,depth_maxN2,'.','markersize',12), hold on
% plot(maxN2(ss),depth_maxN2(ss),'r.','markersize',12), hold on
% plot(maxN2(mm),depth_maxN2(mm),'b.','markersize',12)
xlabel('Max N2 (s^{-2})')
ylabel('Depth of max N2 (m)')
set(gca,'Ydir','reverse','Ygrid','on')
line([2e-5 2e-5],ylim,'color',[0.5 0.5 0.5])

%% Relationship between average N2 and max N2:
subplot(2,4,4,'replace')
plot(avgN2,maxN2,'.','markersize',12), hold on
xlabel('Average N2 (s^{-2})')
ylabel('Max N2 (s^{-2})')

%% Density difference and max N2
subplot(2,4,5,'replace')
% plot(data_at_4m(ss,4)-data_at_12m(ss,4),maxN2(ss),'r.','markersize',12), hold on
% plot(data_at_4m(mm,4)-data_at_12m(mm,4),maxN2(mm),'b.','markersize',12)
plot(data_at_4m(:,4)-data_at_12m(:,4),maxN2(:),'.','markersize',12)
xlabel('Density difference between 4m and 12m (\circC)')
ylabel('Max N2 (s^{-2})')

subplot(2,4,6,'replace')
%plot(data_at_4m(ss,4)-data_at_12m(ss,4),avgN2(ss),'r.','markersize',12), hold on
%plot(data_at_4m(mm,4)-data_at_12m(mm,4),avgN2(mm),'b.','markersize',12)
plot(data_at_4m(:,4)-data_at_12m(:,4),avgN2(:),'.','markersize',12), hold on
xlabel('Density difference between 4m and 12m (\circC)')
ylabel('Average N2 (s^{-2})')
line(xlim,[2e-5 2e-5],'color',[0.5 0.5 0.5])
line(xlim,[5e-5 5e-5],'color',[0.5 0.5 0.5])

%% Temperature difference and max N2:

% subplot(2,4,7,'replace')
% % plot(cell2mat(mldBuse(ss,8))-cell2mat(mldBuse(ss,12)),cell2mat(mldBuse(ss,3)),'r.','markersize',12), hold on
% % plot(cell2mat(mldBuse(mm,8))-cell2mat(mldBuse(mm,12)),cell2mat(mldBuse(mm,3)),'b.','markersize',12)
% plot(cell2mat(mldBuse(:,8))-cell2mat(mldBuse(:,12)),cell2mat(mldBuse(:,3)),'r.','markersize',12)
%
% xlabel('Temperature difference between 4m and 12m (\circC)')
% ylabel('Max N2 (s^{-2})')

%legend('N2 >= 1e-4', 'N2 < 1e-4')
%%
% subplot(2,4,2,'replace')
% plot(cell2mat(mldBuse(ss,2)),cell2mat(mldBuse(ss,9))-cell2mat(mldBuse(ss,13)),'r.','markersize',12), hold on
% plot(cell2mat(mldBuse(mm,2))-cell2mat(mldBuse(mm,12)),cell2mat(mldBuse(mm,9))-cell2mat(mldBuse(mm,13)),'b.','markersize',12)
% ylabel('Density difference between 4m and 12m (\circC)')
% datetick
% set(gca,'xgrid','on')

%% temperature difference between 4-12m and dens difference between 4-12m
subplot(2,4,7,'replace')
% plot(data_at_4m(ss,4)-data_at_12m(ss,4),data_at_4m(ss,3)-data_at_12m(ss,3),'r.','markersize',12), hold on
% plot(data_at_4m(mm,4)-data_at_12m(mm,4),data_at_4m(mm,3)-data_at_12m(mm,3),'b.','markersize',12)
plot(data_at_4m(:,4)-data_at_12m(:,4),data_at_4m(:,3)-data_at_12m(:,3),'.','markersize',12)
ylabel('Temperature difference between 4m and 12m (\circC)')
xlabel('Density difference between 4m and 12m (\circC)')


subplot(2,4,8,'replace')
% plot(avgN2(ss),data_at_4m(ss,3)-data_at_12m(ss,3),'r.','markersize',12), hold on
% plot(avgN2(mm),data_at_4m(mm,3)-data_at_12m(mm,3),'b.','markersize',12)
plot(avgN2(:),data_at_4m(:,3)-data_at_12m(:,3),'.','markersize',12)
ylabel('Temperature difference between 4m and 12m (\circC)')
xlabel('Average N2')

%
% subplot(2,4,8,'replace')
% plot(cell2mat(mldBuse(ss,9))-cell2mat(mldBuse(ss,13)),cell2mat(mldBuse(ss,3)),'r.','markersize',12), hold on
% plot(cell2mat(mldBuse(mm,9))-cell2mat(mldBuse(mm,13)),cell2mat(mldBuse(mm,3)),'b.','markersize',12)
% xlabel('Density difference between 4m and 12m (\circC)')
% ylabel('Max N2 (s^{-2})')

% subplot(2,4,6,'replace')
% plot(cell2mat(mldBuse(ss,8))-cell2mat(mldBuse(ss,12)),cell2mat(mldBuse(ss,3)),'r.','markersize',12), hold on
% plot(cell2mat(mldBuse(mm,8))-cell2mat(mldBuse(mm,12)),cell2mat(mldBuse(mm,3)),'b.','markersize',12)
% xlabel('Temperature difference between 4m and 12m (\circC)')
% ylabel('Max N2 (s^{-2})')


%%
figure
%subplot(2,4,7,'replace')
scatter(cell2mat(mldBuse(:,8))-cell2mat(mldBuse(:,12)),cell2mat(mldBuse(:,9))-cell2mat(mldBuse(:,13)),40,find_yearday(cell2mat(mldBuse(:,2))),'filled')
xlabel('Temperature difference between 4m and 12m (\circC)')
ylabel('Density difference between 4m and 12m (\circC)')

colormap(jet)
colorbar
%in general this is a tight relationship, and gives confidence to use temp
%diff as a proxy for density diff

%%

%but the real question is how does density (or temp diff) for that matter
%correlate to actual stratification?

%Brunt Vaisala freq does show a nice correlation with density difference
X=cell2mat(mldBuse(:,9))-cell2mat(mldBuse(:,13));
Y=cell2mat(mldBuse(:,3));
nn=find(~isnan(X) & ~isnan(Y)); X=X(nn); Y=Y(nn);
[B,BINT,~,~,STATS] = regress(Y,[ones(size(X)) X]);

subplot(2,4,7)
plot(sort(X),B(1)+B(2)*sort(X),'k-')
text(0,0.01,['R2: ' num2str(1e-3*round(1000*STATS(1)))])

subplot(2,4,8)
plot(sort(X),B(1)+B(2)*sort(X),'k-')


%% Brunt Vaisala freq does show a nice correlation with density difference,
%but even more beautiful relationship as log!
X=cell2mat(mldBuse(:,13))-cell2mat(mldBuse(:,9)); %density difference
Y=cell2mat(mldBuse(:,3));
nn=find(~isnan(X) & ~isnan(Y)); X=X(nn); Y=Y(nn);
[xmin,resnorm,residual] = lsqnonlin(@(theta) lsq_lightcurve(theta,X,Y),[-2 0.01],[-10 -100],[10 1000]);


%%
Y_hat=xmin(1)*(1-exp(xmin(2)*sort(X)));
%Y_hat=xmin(1)*(1-exp(-xmin(2)*sort(X)));
%Y_hat=xmin(1)*(1-exp((xmin(2)./sort(X))));
plot(sort(X),Y_hat,'r-')

