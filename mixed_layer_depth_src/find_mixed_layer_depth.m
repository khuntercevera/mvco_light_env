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

downcast_N2=nan(size(downcast_bins));
downcast_N2=downcast_N2(1:end-1,:);

mvco_pdens=nan(4,5,length(mvco_ind));
mvco_delta=nan(3,4,length(mvco_ind));
mvco_temp=nan(4,5,length(mvco_ind));
%%
for q=47:length(mvco_ind);
    
    col_hdr=CTD_QC(mvco_ind(q)).data_hdr;
    temp_data=CTD_QC(mvco_ind(q)).data;
    
    %look at when density or temp crosses delta values from a reference: --------------------------------------------
    pdens_deltas=[0.05 0.1 0.2 0.5];
    temp_deltas=[0.2 0.5 0.8 1];
    
    %reference as close to surface as possibe:
    tt=find(~isnan(downcast_temp(:,mvco_ind(q))));
    temp_ref1=downcast_temp(tt(1),mvco_ind(q)); %first non-nan measurement
    qq=find(~isnan(downcast_pdens(:,mvco_ind(q))));
    pdens_ref1=downcast_pdens(qq(1),mvco_ind(q));
    
    %reference at 4 m:
    rr=find(downcast_bins(:,mvco_ind(q))==4);
    temp_ref2=downcast_temp(rr,mvco_ind(q)); %first non-nan measurement
    qq1=find(~isnan(downcast_pdens(:,mvco_ind(q))));
    pdens_ref2=downcast_pdens(rr,mvco_ind(q));
       
    rec_pdens=[pdens_deltas' repmat(tt(1),size(pdens_deltas))' nan(size(pdens_deltas))' repmat(rr,size(pdens_deltas))' nan(size(pdens_deltas))'];
    rec_temp=[temp_deltas' repmat(qq(1),size(temp_deltas))' nan(size(temp_deltas))' repmat(rr,size(temp_deltas))' nan(size(temp_deltas))'];
    
    for w=1:2
        
        eval(['pdens_ref=pdens_ref' num2str(w) ';'])
        eval(['temp_ref=temp_ref' num2str(w) ';'])
        
        j=rec_pdens(1,2*w); %counter for depth bins
        h=1; %counter for values in pdens threshold
        while j < length(downcast_bins)-1 && h <= length(pdens_deltas)
            
            if abs(pdens_ref - downcast_pdens(j,mvco_ind(q))) > pdens_deltas(1,h)
                rec_pdens(h,2*w+1)=downcast_bins(j,mvco_ind(q));
                h=h+1;
            else
              j=j+1;  
            end
            
        end
        
        j=rec_temp(1,2*w); %counter for depth bins
        g=1; %counter for values in temp threshold
        while j < length(downcast_bins)-1 && g <= length(temp_deltas)
            
            if abs(temp_ref - downcast_temp(j,mvco_ind(q))) > temp_deltas(g)
                rec_temp(g,2*w+1)=downcast_bins(j,mvco_ind(q));
                g=g+1;
            else
                j=j+1;
            end
        end
    end
    
        %density changes based on temp_deltas; just for comparison:
    for k=1:length(temp_deltas)
        pdens_deltas(2,k)=sw_pden(downcast_sal(qq(1),mvco_ind(q)),downcast_temp(qq(1),mvco_ind(q))-temp_deltas(k),downcast_press(qq(1),mvco_ind(q)),0) - ...
            sw_pden(downcast_sal(qq(1),mvco_ind(q)),downcast_temp(qq(1),mvco_ind(q)),downcast_press(qq(1),mvco_ind(q)),0);%-downcast_pdens(qq(1),mvco_ind(q));
        pdens_deltas(3,k)=sw_pden(downcast_sal(rr,mvco_ind(q)),downcast_temp(rr,mvco_ind(q))-temp_deltas(k),downcast_press(rr,mvco_ind(q)),0) - ...
            sw_pden(downcast_sal(rr,mvco_ind(q)),downcast_temp(rr,mvco_ind(q)),downcast_press(rr,mvco_ind(q)),0);% -downcast_pdens(rr,mvco_ind(q));
    end
    
    %calculate N2 from binned data: --------------------------------------------------
    N2=sw_bfrq(downcast_sal(:,mvco_ind(q)),downcast_temp(:,mvco_ind(q)),downcast_press(:,mvco_ind(q)),downcast_lat(mvco_ind(q)));
    downcast_N2(:,mvco_ind(q))=N2;
    avgN2=nanmean(N2);
    medN2=nanmedian(N2);
    [~,im]=max(N2);
    
    %     disp(['Max N2: ' num2str(N2(im))])
    %     disp(['Avg N2: ' num2str(avgN2)])
    %     disp(['Median N2: ' num2str(medN2)])
    
    
    
    %plots!
    %
    clf %see what metrics are highlighting - over all looks pretty good!
    subplot(1,4,1,'replace'),  hold on
    %plot(pdens,depth,'k.-')
    plot(downcast_pdens(:,mvco_ind(q)),downcast_bins(:,mvco_ind(q)),'.-','color','k')
    xlim([min(downcast_pdens(:,mvco_ind(q)))-0.5  max(downcast_pdens(:,mvco_ind(q)))+0.5])
    %xlim([1024 1026])
    line(xlim,[rec_pdens(1,3) rec_pdens(1,3)],'color',[0.8 0.8 0.8])
    line(xlim,[rec_pdens(2,3) rec_pdens(2,3)],'color',[0.6 0.6 0.6])
    line(xlim,[rec_pdens(3,3) rec_pdens(3,3)],'color',[0.4 0.4 0.4])
    line(xlim,[rec_pdens(4,3) rec_pdens(4,3)],'color',[0.2 0.2 0.2])
    
    text(min(downcast_pdens(:,mvco_ind(q)))-0.5, rec_pdens(1,3)-0.1,num2str(pdens_deltas(1,1)))
    text(max(downcast_pdens(:,mvco_ind(q)))+0.5, rec_pdens(2,3)-0.1,num2str(pdens_deltas(1,2)))
    text(min(downcast_pdens(:,mvco_ind(q)))-0.5, rec_pdens(3,3)-0.1,num2str(pdens_deltas(1,3)))
    text(max(downcast_pdens(:,mvco_ind(q)))+0.5, rec_pdens(4,3)-0.1,num2str(pdens_deltas(1,4)))
    
    line(xlim,[rec_pdens(1,5) rec_pdens(1,5)],'linestyle',':','color',[0 0.8 0])
    line(xlim,[rec_pdens(2,5) rec_pdens(2,5)],'linestyle',':','color',[0 0.6 0])
    line(xlim,[rec_pdens(3,5) rec_pdens(3,5)],'linestyle',':','color',[0 0.4 0])
    line(xlim,[rec_pdens(4,5) rec_pdens(4,5)],'linestyle',':','color',[0 0.2 0])
    
    text(min(downcast_pdens(:,mvco_ind(q)))-0.5, rec_pdens(1,5)-0.1,num2str(pdens_deltas(1,1)),'color',[0 0.5 0])
    text(max(downcast_pdens(:,mvco_ind(q)))+0.5, rec_pdens(2,5)-0.1,num2str(pdens_deltas(1,2)),'color',[0 0.5 0])
    text(min(downcast_pdens(:,mvco_ind(q)))-0.5, rec_pdens(3,5)-0.1,num2str(pdens_deltas(1,3)),'color',[0 0.5 0])
    text(max(downcast_pdens(:,mvco_ind(q)))+0.5, rec_pdens(4,5)-0.1,num2str(pdens_deltas(1,4)),'color',[0 0.5 0])
    
    set(gca,'ydir','reverse','fontsize',14)
    xlabel('Potential Density (kg/m^3)') %ylabel(col_hdr{6})
    title([CTD(mvco_ind(q)).cast_name ':' datestr(CTD_QC(mvco_ind(q)).upload_time) ' : ' num2str(q) ' out of ' num2str(length(mvco_ind))],'interpreter','none')
    %legend('Obs \rho','Potential \rho','location','NorthWest')
    
    
    subplot(1,4,2,'replace'), hold on
    plot(downcast_temp(:,mvco_ind(q)),downcast_bins(:,mvco_ind(q)),'.-','color',[1 0.5 0])
    xlim([min(downcast_temp(:,mvco_ind(q)))-0.5  max(downcast_temp(:,mvco_ind(q)))+0.5])
    
    line(xlim,[rec_temp(1,3) rec_temp(1,3)],'color',[0.8 0.8 0.8])
    line(xlim,[rec_temp(2,3) rec_temp(2,3)],'color',[0.6 0.6 0.6])
    line(xlim,[rec_temp(3,3) rec_temp(3,3)],'color',[0.4 0.4 0.4])
    line(xlim,[rec_temp(4,3) rec_temp(4,3)],'color',[0.2 0.2 0.2])
        
    text(min(downcast_temp(:,mvco_ind(q)))-0.4, rec_temp(1,3)-0.1,[num2str(temp_deltas(1,1)) ' - ' num2str(0.001*round(1000*pdens_deltas(2,1)))])
    text(max(downcast_temp(:,mvco_ind(q)))+0.3, rec_temp(2,3)-0.1,[num2str(temp_deltas(1,2)) ' - ' num2str(0.001*round(1000*pdens_deltas(2,2)))])
    text(min(downcast_temp(:,mvco_ind(q)))-0.4, rec_temp(3,3)-0.1,[num2str(temp_deltas(1,3)) ' - ' num2str(0.001*round(1000*pdens_deltas(2,3)))])
    text(max(downcast_temp(:,mvco_ind(q)))+0.3, rec_temp(4,3)-0.1,[num2str(temp_deltas(1,4)) ' - ' num2str(0.001*round(1000*pdens_deltas(2,4)))])
        
    line(xlim,[rec_temp(1,5) rec_temp(1,5)],'linestyle',':','color',[0 0.8 0])
    line(xlim,[rec_temp(2,5) rec_temp(2,5)],'linestyle',':','color',[0 0.6 0])
    line(xlim,[rec_temp(3,5) rec_temp(3,5)],'linestyle',':','color',[0 0.4 0])
    line(xlim,[rec_temp(4,5) rec_temp(4,5)],'linestyle',':','color',[0 0.2 0])
    
    text(max(downcast_temp(:,mvco_ind(q)))+0.3, rec_temp(1,5)-0.1,[num2str(temp_deltas(1,1)) ' - ' num2str(0.001*round(1000*pdens_deltas(3,1)))],'color',[0 0.5 0])
    text(min(downcast_temp(:,mvco_ind(q)))-0.4, rec_temp(2,5)-0.1,[num2str(temp_deltas(1,2)) ' - ' num2str(0.001*round(1000*pdens_deltas(3,2)))],'color',[0 0.5 0])
    text(max(downcast_temp(:,mvco_ind(q)))+0.3, rec_temp(3,5)-0.1,[num2str(temp_deltas(1,3)) ' - ' num2str(0.001*round(1000*pdens_deltas(3,3)))],'color',[0 0.5 0])
    text(min(downcast_temp(:,mvco_ind(q)))-0.4, rec_temp(4,5)-0.1,[num2str(temp_deltas(1,4)) ' - ' num2str(0.001*round(1000*pdens_deltas(3,4)))],'color',[0 0.5 0])
      
    set(gca,'ydir','reverse','fontsize',14)
    xlabel('Temperature')
    title('Temperature with depth')
    
    subplot(1,4,3,'replace'), hold on
    plot(downcast_sal(:,mvco_ind(q)),downcast_bins(:,mvco_ind(q)),'.-','color',[0 0.5 1])
    xlim([min(downcast_sal(:,mvco_ind(q)))-0.5  max(downcast_sal(:,mvco_ind(q)))+0.5])
    set(gca,'ydir','reverse','fontsize',14)
    xlabel('Salinity')
    title('Salinity with depth')
    
    subplot(1,4,4,'replace'), hold on
    plot(downcast_N2(:,mvco_ind(q)),downcast_bins(1:end-1,mvco_ind(q)),'.-','color',[0 0.7 0])
    set(gca,'ydir','reverse','fontsize',14)
    xlabel('N2')
    title('Brunt-Vaisala with depth')
    
    line([1e-4 1e-4], get(gca,'ylim'),'linestyle',':','color','r') %anything higher than this and you are probably stratified
    line([avgN2 avgN2], get(gca,'ylim'),'color',[0.5 0.5 0.5]) %average N2
    line([medN2 medN2], get(gca,'ylim'),'color',[0.2 0.2 0.2]) %average N2
    plot(N2(im),downcast_bins(im,mvco_ind(q)),'pk','markersize',12) %max only considering below 3.75m
    
    
    %Is the water stratified? Or wouldn't be well mixed?
    %Looking at Young-Oh's slide, it looks like the min would be around
    %.25 * 10^-4 N2 (S^{-2})...
    
    
    %
    keyboard
    
    mvco_pdens(:,:,q)=rec_pdens;
    mvco_temp(:,:,q)=rec_temp;
    mvco_delta(:,:,q)=pdens_deltas;
end %for loop


%% and now to plot!

%okay, it seems we have a few types of casts:
%   surface mini layer, where below 3/4 m, it would be well mixed
%   a mixed layer somewhere mid depth
%   just a slow, 14m almost of stratification
%   nice and fully mixed :)

%it does seem that a cutoff of 0.2 for a delta in density captures a 'mixed
%layer nicely':

%identify each type of water column:
time=cell2mat({CTD_QC(mvco_ind).upload_time}');
yrdy=find_yearday(time);
surf=squeeze(mvco_pdens(3,3,:)); %from surface ref
four=squeeze(mvco_pdens(3,5,:)); %from 4m ref

%% well mixed:
mm1=find(isnan(surf) & isnan(four));
mm2=find(surf >= 12 | four >= 12); %ends up being captured just by four

%only surface piece stratified:
ss=find(surf < 4 & isnan(four));

%layers somewhere in the middle:
ms=find((surf >= 4 & surf < 12) | four < 12);


%alright, have all the points!

clf, hold on
plot(yrdy(mm1),15*ones(size(mm1)),'b.','markersize',12)
plot(yrdy(mm2),four(mm2),'b.','markersize',12)
plot(yrdy(ss),surf(ss),'.','markersize',12,'color',[1 0.3 0])
plot(yrdy(ms),four(ms),'.','markersize',12,'color',[0.5 0.5 0.5])


%% okay, and now the density differences at 4m and 12m:

i4=find(downcast_bins(:,1)==4);
i12=find(downcast_bins(:,1)==12);

figure
subplot(1,3,1,'replace')
hist(downcast_pdens(i12,mvco_ind([mm1;mm2]))-downcast_pdens(i4,mvco_ind([mm1;mm2])));
title('well mixed')
subplot(1,3,2,'replace')
hist(downcast_pdens(i12,mvco_ind(ms))-downcast_pdens(i4,mvco_ind(ms)));
title('mixed layer at depth')
subplot(1,3,3,'replace')
hist(downcast_pdens(i12,mvco_ind(ss))-downcast_pdens(i4,mvco_ind(ss)));
title('surface stratification')

%%

figure
subplot(1,3,1,'replace')
hist(downcast_temp(i12,mvco_ind([mm1;mm2]))-downcast_temp(i4,mvco_ind([mm1;mm2])));
title('well mixed')
subplot(1,3,2,'replace')
hist(downcast_temp(i12,mvco_ind(ms))-downcast_temp(i4,mvco_ind(ms)));
title('mixed layer at depth')
subplot(1,3,3,'replace')
hist(downcast_temp(i12,mvco_ind(ss))-downcast_temp(i4,mvco_ind(ss)));
title('surface stratification')
%%

figure, subplot(1,2,1,'replace')
colormap jet
pdens_diff=[downcast_pdens(i12,mvco_ind)-downcast_pdens(i4,mvco_ind)]';
temp_diff=[downcast_temp(i12,mvco_ind)-downcast_temp(i4,mvco_ind)]';

[b,~,~,~,stats]=regress(temp_diff,[ones(length(pdens_diff),1) pdens_diff]);
scatter(pdens_diff,temp_diff,30,yrdy,'filled')

line([0 1],[b(1) b(1)+b(2)])

subplot(1,2,2,'replace'), hold on
plot(0.2,squeeze(mvco_delta(3,1,:)),'b.')
plot(0.5,squeeze(mvco_delta(3,2,:)),'k.')
plot(0.8,squeeze(mvco_delta(3,3,:)),'r.')
plot(1,squeeze(mvco_delta(3,4,:)),'g.')
title('calculated denisty changes for a given temperature change')

plot(0.2,nanmean(squeeze(mvco_delta(3,1,:))),'bp')
plot(0.5,nanmean(squeeze(mvco_delta(3,2,:))),'kp')
plot(0.8,nanmean(squeeze(mvco_delta(3,3,:))),'rp')
plot(1,nanmean(squeeze(mvco_delta(3,4,:))),'gp')

%[b2,~,~,~,stats]=regress(temp_deltas',[ones(4,1) nanmean((mvco_delta(3,:,:)),3)'])

%% climatology difference over time:

load /Volumes/Lab_data/MVCO/EnvironmentalData/Tday_beam.mat %from the MVCO tower beam
eval('beam_years=yearlist;')
eval('beam_date=mdate_mat;')
load /Volumes/Lab_data/MVCO/EnvironmentalData/Tall_day.mat %from undersea node
eval('node_date=mdate;')
eval('node_years=year;')
eval('Tday_node=Tday;')
%%
beam_avg=nanmean(Tday_beam,2);
node_avg=nanmean(Tday_node(:,4:end),2);
avg_diff=nanmean((Tday_beam-Tday_node(:,4:end)),2);
%%
clf, plot(1:366,avg_diff,'.'), hold on
plot(1:366,beam_avg-node_avg,'.')