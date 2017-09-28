%Okay, before can use CTD casts, need to QC data and remove troubling start
%points and residual fresh water in salinity sensor at start:

%for calculating potential density (pressure affect removed):
addpath /Users/kristenhunter-cevera/MVCO_light_at_depth/seawater_ver3_2/
addpath /Users/kristenhunter-cevera/Documents/MATLAB/mvco_tools/
mldB={};
plotflag=1;
warning off

mld_titlesB={'cast name'; 'file_time'; 'max N2'; 'depth of max N2'; 'mean N2'; ...
                'closest depth to 4m'; 'N2 closest to 4m'; 'temperature closest to 4m'; 'pdens closest to 4m'; ... %depth, N2, temp and dens at 4m
                'closest depth to 12m'; 'N2 closest to 12m'; 'temperature closest to 12m'; 'pdens closest to 12m'; 'upcast_flag';'manual edit'};  %depth, N2, temp and dens at 12m
          

%mld_titles={'cast name'; 'file_time'; 'max N2'; 'depth of max N2'; 'stratified?'; ...
%                 'closest depth to 4m'; 'N2 closest to 4m'; 'temperature closest to 4m'; 'pdens closest to 4m'; ... %depth, N2, temp and dens at 4m
%                 'closest depth to 12m'; 'N2 closest to 12m'; 'temperature closest to 12m'; 'pdens closest to 12m'; 'upcast_flag';'manual edit'};  %depth, N2, temp and dens at 12m
%           
%%
for q=45:length(mvco_ind);
    
    col_hdr=CTD(mvco_ind(q)).data_hdr;
    temp_data=CTD(mvco_ind(q)).data;
    
    if any(~cellfun('isempty',temp_data))
        
        %find and unload desired variables from cell array:
        dp=find(cellfun('isempty',regexp( col_hdr,'Depth'))==0);
        s=find(cellfun('isempty',regexp( col_hdr,'Salinity'))==0);
        t=find(cellfun('isempty',regexp( col_hdr,'Temperature'))==0);
        p=find(cellfun('isempty',regexp( col_hdr,'Pressure'))==0);
        de=find(cellfun('isempty',regexp( col_hdr,'[^(Potential] Density'))==0);
        pd=find(cellfun('isempty',regexp( col_hdr,'Potential'))==0);
        ct=find(cellfun('isempty',regexp(col_hdr,'Time, Elapsed'))==0);
        
        depth=temp_data{dp};
        press=temp_data{p};
        temperature=temp_data{t};
        sal=temp_data{s};
        dens=temp_data{de};
        pdens=temp_data{pd};
        cast_time=temp_data{ct};
        
        file_time=CTD(mvco_ind(q)).upload_time;
        
        %Now to QC data:
        %Remove any salinities with less than 25 ppt
        %typically there is a time delay at the beginning before CTD starts
        %dropping, go with data that's after at least a 30s delay...
        %occasionally, there are time points out of order...resort
        
        %find descent:
        smooth_depth=smooth(depth,200); %just a bit easier for direct logic tests
        dc=find(diff(smooth_depth) > 0.0025); %falling
        
        %find longest continuous stretch-this should be the descent:
        nonconsec=find(diff(dc)~=1); %find the values that are not contiguous
        [mm, ind]=max(diff(nonconsec));
        dsc_ind1=dc(nonconsec(ind)+1);
        dsc_ind2=dc(nonconsec(ind+1));
        dsc=dsc_ind1:dsc_ind2; %indexes for main descent
        usc=dsc_ind2+1:length(depth); %and the upcast!
        
        %calc Brunt-Vaisala frequency to help determine MLD, but first bin
        if ~isempty(dc) && max(depth) > 4 && length(dsc) > 1
            
            %for easier handling:
            depthD=depth(dsc);
            salD=sal(dsc);
            temperatureD=temperature(dsc);
            pressD=press(dsc);
            pdensD=pdens(dsc);
            cast_timeD=cast_time(dsc);
            
            depthU=depth(usc);
            salU=sal(usc);
            temperatureU=temperature(usc);
            pressU=press(usc);
            pdensU=pdens(usc);
            cast_timeU=cast_time(usc);
            
            
            %bin variables into ~0.25 m bins:
            max_depth=max([depthD;depthU]);
            min_depth=min([depthD;depthU]);
            
            depth_bins=floor(min_depth):0.25:ceil(max_depth);
            binned_dataD=nan(length(depth_bins),6);
            binned_dataU=nan(length(depth_bins),6);
            for i=1:length(depth_bins)-1
                ii=find(depthD >= depth_bins(i) & depthD < depth_bins(i+1));
                binned_dataD(i,1)=depth_bins(i);
                binned_dataD(i,2)=nanmean(depthD(ii));
                binned_dataD(i,3)=nanmean(salD(ii));
                binned_dataD(i,4)=nanmean(temperatureD(ii));
                binned_dataD(i,5)=nanmean(pressD(ii));
                binned_dataD(i,6)=nanmean(pdensD(ii));
                
                jj=find(depthU >= depth_bins(i) & depthU < depth_bins(i+1));
                binned_dataU(i,1)=depth_bins(i);
                binned_dataU(i,2)=nanmean(depthU(jj));
                binned_dataU(i,3)=nanmean(salU(jj));
                binned_dataU(i,4)=nanmean(temperatureU(jj));
                binned_dataU(i,5)=nanmean(pressU(jj));
                binned_dataU(i,6)=nanmean(pdensU(jj));
            end
            
            %calculate N2 from binned data:
            N2D=sw_bfrq(binned_dataD(:,3),binned_dataD(:,4),binned_dataD(:,5));
            N2U=sw_bfrq(binned_dataU(:,3),binned_dataU(:,4),binned_dataU(:,5));
            %sw_bfrq(sal (psu), temperature (deg C), pressure (db))
            %we need to see what the max N2 is, and what N2 is at 4m...
            
            %If salinity data doesn't really match on the upcast, then exclude this in the downcast:
            %maybe only consider N2 values that are below 3m to avoid false positives:
            
            delta_salU=max(salU(find(depthU < 4)))-min(salU(find(depthU < 4)));
            delta_salD=max(salD(find(depthD < 4)))-min(salD(find(depthD < 4)));
            
            if abs(delta_salD./delta_salU) > 2
                disp('Salinity on the down cast questionable...going with upcast for now...')
                upcast_flag=1;
                binned_data=binned_dataU;
                N2=N2U;
            else
                disp('Salinity difference is reasonable, going with downcast for now...')
                binned_data=binned_dataD;
                N2=N2D;
                upcast_flag=0;
            end
            
            i3=find(binned_data(1:end-1,1) > 3);
            [mm, is]=sort(N2(i3),'descend');
            nn=find(~isnan(mm)); %avoid nan's
            mm=mm(nn); is=is(nn); %to make indexing cleaner
            im=i3(is(1)); %index corresponds to max N2 value below 3 m...
            
            %mean N2 over data:
            avgN2=nanmean(N2(i3));
            
            %and find the indicies for 4m and 12m (or closet depth):
            [d4,ii4]=min(abs(binned_data(:,1)-4));
            [d12,ii12]=min(abs(binned_data(:,1)-12));
            
        end
        
        % Hmm, need to add column of indexes, instabilities (N2 < 0) and median?
        
        if plotflag==1
            clf %see what this metric is highlighting - over all looks pretty good!
            subplot(2,3,1,'replace'),  hold on
            %plot(pdens,depth,'k.-')
            plot(pdensD,depthD,'.','color',[0 0.5 1])
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
            
            %pause
        end %plotflag
               
        %Is the water stratified? Or wouldn't be well mixed?
        %Looking at Young-Oh's slide, it looks like the min would be around
        %.25 * 10^-4 N2 (S^{-2})...
        
        %mld2=[mld2; {CTD(mvco_ind(q)).cast_name} {file_time} {below3_N2(im)} {below3_bins(im})];
            disp(['Max N2: ' num2str(N2(im))])
            disp(['Avg N2: ' num2str(avgN2)])
            
        [Q] = input('Does the automatic finding of N2 look reasonable? If so, enter "auto", if not, enter "manual" \n');
        
        if strcmp(Q,'auto')
        mldB=[mldB; {CTD(mvco_ind(q)).cast_name} {file_time} {N2(im)} {binned_data(im,1)} {avgN2} ...
                {binned_data(ii4,1)} {N2(ii4)} {binned_data(ii4,4)} {binned_data(ii4,6)} ... %depth, N2, temp and dens at 4m
                {binned_data(ii12,1)} {N2(ii12)} {binned_data(ii12,4)} {binned_data(ii12,6)} {upcast_flag} {'auto'}];  %depth, N2, temp and dens at 12m

        elseif strcmp(Q,'manual')
            disp('Please adjust and enter in data for mld cell array')
            keyboard
%          mldB=[mldB; {CTD(mvco_ind(q)).cast_name} {file_time} {N2(im)} {binned_data(im,1)} {avgN2} ...
%                 {binned_data(ii4,1)} {N2(ii4)} {binned_data(ii4,4)} {binned_data(ii4,6)} ... %depth, N2, temp and dens at 4m
%                 {binned_data(ii12,1)} {N2(ii12)} {binned_data(ii12,4)} {binned_data(ii12,6)} {upcast_flag} {'auto'}];  %depth, N2, temp and dens at 12m
%      
        end
            
%         if N2(im) <= 3e-4 %well mixed:
%             mldB=[mldB; {CTD(mvco_ind(q)).cast_name} {file_time} {N2(im)} {binned_data(im,1)} {'mixed'} ...
%                 {binned_data(ii4,1)} {N2(ii4)} {binned_data(ii4,4)} {binned_data(ii4,6)} ... %depth, N2, temp and dens at 4m
%                 {binned_data(ii12,1)} {N2(ii12)} {binned_data(ii12,4)} {binned_data(ii12,6)} {upcast_flag} {'auto'}];  %depth, N2, temp and dens at 12m
%             disp('Mixed!')
%             disp(num2str(N2(im)))
%         else
%             mldB=[mldB; {CTD(mvco_ind(q)).cast_name} {file_time} {N2(im)} {binned_data(im,1)} {'stratified'} ...
%                 {binned_data(ii4,1)} {N2(ii4)} {binned_data(ii4,4)} {binned_data(ii4,6)} ... %depth, N2, temp and dens at 4m
%                 {binned_data(ii12,1)} {N2(ii12)} {binned_data(ii12,4)} {binned_data(ii12,6)} {upcast_flag} {'auto'}];  %depth, N2, temp and dens at 12m
%             disp('Stratified!')
%             disp(num2str(N2(im)))
%         end
          
        
        if any(diff(cast_time) < 0), disp('Something wrong with time sync?'), end
        clear temp_data col_hdr
        
    end %if not empty
    
end %for loop

%mld_hdr={'cast name';'time of file';'max N2 (<3.75m)'; 'Depth for ax N2'; 'call on water column';'reason';'4m temp';'12m temp';'4m pdens';'12m pdens';'If not 12m, lowest depth'};



%% Okay and now some plots to try to make sense of this mess...

%Compare difference of temp to diff of density, color coded for stratified
%or not....

%first remove empty rows where something went awry:
mld=mldB;
jj=find(cellfun('isempty',mld(:,3))==0);
%ii=find(cellfun('isempty',mld(:,7))==0 & cellfun('isempty',mld(:,8))==0 & cellfun('isempty',mld(:,9))==0 & cellfun('isempty',mld(:,10))==0);
ii=find(abs(cell2mat(mld(jj,10))-12) < 1);

%%
mldBuse=mld(jj(ii),:);
% ss=find(cellfun('isempty',regexp(mldBuse(:,5),'stratified'))==0);
% mm=find(cellfun('isempty',regexp(mldBuse(:,5),'stratified'))==1);

ss=find(cell2mat(mldBuse(:,3)) >= 1e-5);
mm=find(cell2mat(mldBuse(:,3)) < 1e-5);

%% by average:
ss=find(cell2mat(mldBuse(:,5)) >= 1e-5);
mm=find(cell2mat(mldBuse(:,5)) < 1e-5);
%% diff in temp vs. diff in dens

clf
% subplot(2,4,1,'replace')
% plot(cell2mat(mldBuse(ss,2)),cell2mat(mldBuse(ss,3)),'r.','markersize',12), hold on
% plot(cell2mat(mldBuse(mm,2)),cell2mat(mldBuse(mm,3)),'b.','markersize',12)
% ylabel('Max N2 (s^{-2})')
% line(xlim,[3e-4 3e-4],'color',[0.5 0.5 0.5])
% datetick

subplot(2,4,1,'replace')
plot(find_yearday(cell2mat(mldBuse(ss,2))),cell2mat(mldBuse(ss,3)),'r.','markersize',12), hold on
plot(find_yearday(cell2mat(mldBuse(mm,2))),cell2mat(mldBuse(mm,3)),'b.','markersize',12)
ylabel('Max N2 (s^{-2})')
xlim([0 365])
xlabel('Yearday')
line(xlim,[1e-4 1e-4],'color',[0.5 0.5 0.5])

subplot(2,4,2,'replace')
plot(cell2mat(mldBuse(ss,2)),cell2mat(mldBuse(ss,9))-cell2mat(mldBuse(ss,13)),'r.','markersize',12), hold on
plot(cell2mat(mldBuse(mm,2))-cell2mat(mldBuse(mm,12)),cell2mat(mldBuse(mm,9))-cell2mat(mldBuse(mm,13)),'b.','markersize',12)
ylabel('Density difference between 4m and 12m (\circC)')
datetick
set(gca,'xgrid','on')

subplot(2,4,7,'replace')
plot(cell2mat(mldBuse(ss,9))-cell2mat(mldBuse(ss,13)),cell2mat(mldBuse(ss,3)),'r.','markersize',12), hold on
plot(cell2mat(mldBuse(mm,9))-cell2mat(mldBuse(mm,13)),cell2mat(mldBuse(mm,3)),'b.','markersize',12)
xlabel('Density difference between 4m and 12m (\circC)')
ylabel('Max N2 (s^{-2})')

subplot(2,4,8,'replace')
plot(cell2mat(mldBuse(ss,9))-cell2mat(mldBuse(ss,13)),cell2mat(mldBuse(ss,3)),'r.','markersize',12), hold on
plot(cell2mat(mldBuse(mm,9))-cell2mat(mldBuse(mm,13)),cell2mat(mldBuse(mm,3)),'b.','markersize',12)
xlabel('Density difference between 4m and 12m (\circC)')
ylabel('Max N2 (s^{-2})')


subplot(2,4,3,'replace')
plot(cell2mat(mldBuse(ss,3)),cell2mat(mldBuse(ss,4)),'r.','markersize',12), hold on
plot(cell2mat(mldBuse(mm,3)),cell2mat(mldBuse(mm,4)),'b.','markersize',12)
xlabel('Max N2 (s^{-2})')
ylabel('Depth of max N2 (m)')
set(gca,'Ydir','reverse')

subplot(2,4,6,'replace')
plot(cell2mat(mldBuse(ss,8))-cell2mat(mldBuse(ss,12)),cell2mat(mldBuse(ss,3)),'r.','markersize',12), hold on
plot(cell2mat(mldBuse(mm,8))-cell2mat(mldBuse(mm,12)),cell2mat(mldBuse(mm,3)),'b.','markersize',12)
xlabel('Temperature difference between 4m and 12m (\circC)')
ylabel('Max N2 (s^{-2})')

subplot(2,4,5,'replace')
plot(cell2mat(mldBuse(ss,8))-cell2mat(mldBuse(ss,12)),cell2mat(mldBuse(ss,9))-cell2mat(mldBuse(ss,13)),'r.','markersize',12), hold on
plot(cell2mat(mldBuse(mm,8))-cell2mat(mldBuse(mm,12)),cell2mat(mldBuse(mm,9))-cell2mat(mldBuse(mm,13)),'b.','markersize',12)
xlabel('Temperature difference between 4m and 12m (\circC)')
ylabel('Density difference between 4m and 12m (\circC)')

legend('N2 >= 1e-4', 'N2 < 1e-4')
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

