%Okay, before can use CTD casts, need to QC data and remove troubling start
%points and residual fresh water in salinity sensor at start:

%for calculating potential density (pressure affect removed):
addpath /Users/kristenhunter-cevera/MVCO_light_at_depth/seawater_ver3_2/
mld2={};

%%
for q=1:length(mvco_ind);
    
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
            max_depth=max(depthD);
            min_depth=min(depthD);
            
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
            %we need to see what the max N2 is, and what N2 is at 4m...
            
            %If salinity data doesn't really match on the upcast, then exclude this in the downcast:
            %maybe only consider N2 values that are below 3m to avoid false positives:
                
            delta_salU=max(salU(find(depthU < 4)))-min(salU(find(depthU < 4)));
            delta_salD=max(salD(find(depthD < 4)))-min(salD(find(depthD < 4)));
            
            if abs(delta_salD./delta_salU) > 2
                disp('Salinity on the down cast questionable...going with upcast for now...')
                upcast_flag=1;
                binne_data=binned_dataU;
                N2=N2U;
            else
                binne_data=binned_dataD;
                N2=N2D;
                upcast_flag=0;
            end
            
            i3=find(binned_dataD(1:end-1,1) > 3);
            [mm, is]=sort(N2(i3),'descend');
            nn=find(~isnan(mm)); %avoid nan's
            mm=mm(nn); is=is(nn); %to make indexing cleaner
            im=i3(nn(is(1))); %index corresponds to max N2 value below 3 m...
            
            %and find the indicies for 4m and 12m (or closet depth):
            [d4,ii4]=min(abs(binned_dataD(:,1)-4));
            [d12,ii12]=min(abs(binned_dataD(:,1)-12));
            
        end
        
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
        title([CTD(mvco_ind(q)).cast_name ':' datestr(file_time)],'interpreter','none')
        %legend('Obs \rho','Potential \rho','location','NorthWest')
        
        subplot(2,3,2,'replace'), hold on
        plot(cast_time,depth,'k.-')
        plot(cast_timeD,depthD,'.','color',[0 0.5 1])
        plot(cast_timeU,depthU,'.','color',[1 0.5 0])
        set(gca,'ydir','reverse','fontsize',14)
        xlabel('Time') %ylabel(col_hdr{6})
        title('CTD position with time')
        
        subplot(2,3,3,'replace'), hold on
        plot(N2,binned_dataD(1:end-1,1),'k.-')
        plot(N2U,binned_dataU(1:end-1,1),'r.-')
        line([2.5e-4 2.5e-4], get(gca,'ylim'),'color','r') %anything higher than this and you are probably stratified
        plot(N2(im),binned_dataD(im,1),'rp') %max only considering below 3.75m
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
        
        %keyboard
        
        %Is the water stratified? Or wouldn't be well mixed?
        %Looking at Young-Oh's slide, it looks like the min would be around
        %.25 * 10^-4 N2 (S^{-2})...
        
        %mld2=[mld2; {CTD(mvco_ind(q)).cast_name} {file_time} {below3_N2(im)} {below3_bins(im})];
        
        if N2(im) <= 2.5e-4 %well mixed:
            mld2=[mld2; {CTD(mvco_ind(q)).cast_name} {file_time} {N2(im)} {binned_data(im,1)} {'mixed'} ...
                {binned_data(ii4,1)} {N2(ii4)} {binned_data(ii4,4)} {binned_data(ii4,6)} ... %depth, N2, temp and dens at 4m
                {binned_data(ii12,1)} {N2(ii12)} {binned_data(ii12,4)} {binned_data(ii12,6)} {upcast_flag}];  %depth, N2, temp and dens at 12m
        else
            mld2=[mld2; {CTD(mvco_ind(q)).cast_name} {file_time} {N2(im)} {binned_data(im,1)} {'stratified'} ...
                {binned_data(ii4,1)} {N2(ii4)} {binned_data(ii4,4)} {binned_data(ii4,6)} ... %depth, N2, temp and dens at 4m
                {binned_data(ii12,1)} {N2(ii12)} {binned_data(ii12,4)} {binned_data(ii12,6)} {upcast_flag}];  %depth, N2, temp and dens at 12m
        end
        
        
    end
    
    pause
    
    if any(diff(cast_time) < 0), disp('Something wrong with time sync?'), end
    clear temp_data col_hdr
    
end

%mld_hdr={'cast name';'time of file';'max N2 (<3.75m)'; 'Depth for ax N2'; 'call on water column';'reason';'4m temp';'12m temp';'4m pdens';'12m pdens';'If not 12m, lowest depth'};



%% Okay and now some plots to try to make sense of this mess...

%Compare difference of temp to diff of density, color coded for stratified
%or not....

%first remove empty rows where something went awry:
ii=find(cellfun('isempty',mld(:,7))==0 & cellfun('isempty',mld(:,8))==0 & cellfun('isempty',mld(:,9))==0 & cellfun('isempty',mld(:,10))==0);
%%
mld2use=mld(ii,:);
ss=find(cellfun('isempty',regexp(mld2use(:,5),'stratified'))==0);
mm=find(cellfun('isempty',regexp(mld2use(:,5),'stratified'))==1);
%%
plot(cell2mat(mld2use(ss,7))-cell2mat(mld2use(ss,8)),cell2mat(mld2use(ss,9))-cell2mat(mld2use(ss,10)),'r.'), hold on
plot(cell2mat(mld2use(mm,7))-cell2mat(mld2use(mm,8)),cell2mat(mld2use(mm,9))-cell2mat(mld2use(mm,10)),'b.')

%%

scatter(cell2mat(mld2use(ss,7))-cell2mat(mld2use(ss,8)),cell2mat(mld2use(ss,9))-cell2mat(mld2use(ss,10)),50,cell2mat(mld2use(ss,3)),'filled'), hold on
scatter(cell2mat(mld2use(mm,7))-cell2mat(mld2use(mm,8)),cell2mat(mld2use(mm,9))-cell2mat(mld2use(mm,10)),50,cell2mat(mld2use(mm,3)),'filled')
plot(cell2mat(mld2use(ss,7))-cell2mat(mld2use(ss,8)),cell2mat(mld2use(ss,9))-cell2mat(mld2use(ss,10)),'ro'), hold on
plot(cell2mat(mld2use(mm,7))-cell2mat(mld2use(mm,8)),cell2mat(mld2use(mm,9))-cell2mat(mld2use(mm,10)),'bo')
colorbar

%%
clf
plot(cell2mat(mld2use(mm,9))-cell2mat(mld2use(mm,10)),cell2mat(mld2use(mm,3)),'bo')
hold on
plot(cell2mat(mld2use(ss,9))-cell2mat(mld2use(ss,10)),cell2mat(mld2use(ss,3)),'ro')
xlabel('Density difference between 4m and 12m')
ylabel('Max N2 somewhere along water column')


%% OLDER MATERIAL BEFORE JULY 2017:
%% fancy plot:
clf, hold on
for q=1:length(mld)
    plot(mld{q,2},mld{q,3},'b.','markersize',14)
    if strcmp(mld(q,end),'up to')
        line([mld{q,2} mld{q,2}],[mld{q,3} 16],'color',[1 0.3 0.3]) %up to
    else
        line([mld{q,2} mld{q,2}],[0 mld{q,3}],'color',[0.3 0.3 0.3])    %down to
    end
end
datetick('x','mmm/yy')
set(gca,'Ydir','reverse')

%% we also need to extract temperature records at 4m and 12m...

for q=1:length(mld)
    
    ind=mld{q,1}; %index into data
    
    col_hdr=hdr{ind};
    temp_data=data{ind};
    
    
    %unload variables from cell array:
    depth=temp_data{1};
    press=temp_data{2};
    temperature=temp_data{3};
    sal=temp_data{5};
    dens=temp_data{6};
    pdens=temp_data{end};
    time=file_time{j,3};
    
    %find time
    ii=find(cellfun('isempty',regexp(col_hdr,'Time, Elapsed'))==0);
    cast_time=temp_data{ii};
    
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
    
    %calc Brunt-Vaisala frequency to help determine MLD, but first bin
    if ~isempty(dc)
        
        %for easier handling:
        depthD=depth(dsc);
        salD=sal(dsc);
        temperatureD=temperature(dsc);
        pressD=press(dsc);
        pdensD=pdens(dsc);
        cast_timeD=cast_time(dsc);
        
        ss=find(depthD>=4);
        if ~isempty(ss)
            mld{q,6}=temperatureD(ss(1));
            mld{q,7}=pdensD(ss(1));
        else
            mld{q,6}=NaN;
            mld{q,7}=NaN;
        end
        
        qq=find(depthD>=12);
        if ~isempty(qq)
            mld{q,8}=temperatureD(qq(1));
            mld{q,9}=pdensD(qq(1));
        else
            mld{q,8}=NaN;
            mld{q,9}=NaN;
        end
    end
end

%%
qq=find(cellfun(@(x) x< 14,mld(:,3))~=0); %something with a stratification happening!
clf, hold on
plot(cell2mat(mld(:,7))-cell2mat(mld(:,9)),cell2mat(mld(:,6))-cell2mat(mld(:,8)),'.','markersize',12)
plot(cell2mat(mld(qq,7))-cell2mat(mld(qq,9)),cell2mat(mld(qq,6))-cell2mat(mld(qq,8)),'.','markersize',12)

%%
clf
scatter(cell2mat(mld(:,7))-cell2mat(mld(:,9)),cell2mat(mld(:,6))-cell2mat(mld(:,8)),40,cell2mat(mld(:,3)),'filled')
