%Okay, before can use CTD casts, need to QC data and remove troubling start
%points and residual fresh water in salinity sensor at start:


%for calculating potential density (pressure affect removed):
addpath /Users/kristenhunter-cevera/Dropbox/MVCO_mixed_layer_depth/seawater_ver3_2/
mld={};
%%
for q=1:length(tower_ind);
    
    j=tower_ind(q);
    col_hdr=hdr{j};
    temp_data=data{j};
    
    if any(~cellfun('isempty',temp_data))
        
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
            
            %bin variables into ~0.25 m bins:
            max_depth=max(depthD);
            min_depth=min(depthD);
            
            depth_bins=floor(min_depth):0.25:ceil(max_depth);
            binned_data=nan(length(depth_bins),6);
            for i=1:length(depth_bins)-1
                ii=find(depthD >= depth_bins(i) & depthD < depth_bins(i+1));
                binned_data(i,1)=depth_bins(i);
                binned_data(i,2)=nanmean(depthD(ii));
                binned_data(i,3)=nanmean(salD(ii));
                binned_data(i,4)=nanmean(temperatureD(ii));
                binned_data(i,5)=nanmean(pressD(ii));
                binned_data(i,6)=nanmean(pdensD(ii));
            end
            
            N2=sw_bfrq(binned_data(:,3),binned_data(:,4),binned_data(:,5));
            [mm, is]=sort(N2,'descend');
            nn=~isnan(mm);
            mm=mm(nn); is=is(nn);
            im=is(1);
            
            
            %            %N2=sw_bfrq(smooth(sal(dsc),100),smooth(temperature(dsc),100),smooth(press(dsc),100));
            %
            %         ii=find(depthD > 3 & depth2 < max_depth-1);
            %         [mm, im]=max(N2(ii(1:end-1)));
            
            clf %see what this metric is highlighting - over all looks pretty good!
            subplot(2,3,1,'replace'),  hold on
            %plot(pdens,depth,'k.-')
            plot(pdensD,depthD,'.','color',[0 0.5 1])
            plot(binned_data(:,6),binned_data(:,1),'.-')
            line(xlim,[binned_data(im,1) binned_data(im,1)],'color','r')
            set(gca,'ydir','reverse','fontsize',14)
            xlabel('Density (kg/m^3)') %ylabel(col_hdr{6})
            title([file_time{j,2} ' - ' datestr(file_time{j,3})])
            %legend('Obs \rho','Potential \rho','location','NorthWest')
            
            subplot(2,3,2,'replace'), hold on
            %plot(cast_time,depth,'k.-')
            plot(cast_timeD,depthD,'.','color',[0 0.5 1])
            set(gca,'ydir','reverse','fontsize',14)
            xlabel('Time') %ylabel(col_hdr{6})
            
            subplot(2,3,3,'replace'), hold on
            plot(N2,binned_data(1:end-1,1),'k.-')
            line([1e-4 1e-4], get(gca,'ylim'),'color','r')
            plot(N2(im),binned_data(im,1),'rp')
            set(gca,'ydir','reverse','fontsize',14)
            xlabel('Brunt-Vaisala Freq')
            xlim([-1e-3 1e-2])
            
            subplot(2,3,4,'replace'),  hold on
            %plot(temperature,depth,'k.-')
            plot(temperatureD,depthD,'.','color',[0 0.5 1])
            plot(binned_data(:,4),binned_data(:,1),'.-')
            set(gca,'ydir','reverse','fontsize',14)
            xlabel('Temperature (\circC)') %ylabel(col_hdr{6})
            
            subplot(2,3,5,'replace'),  hold on
            plot(sal,depth,'k.-')
            plot(salD,depthD,'.','color',[0 0.5 1])
            plot(binned_data(:,3),binned_data(:,1),'.-')
            set(gca,'ydir','reverse','fontsize',14)
            xlabel('Salinity') %ylabel(col_hdr{6})
            xlim([median(salD)-1 median(salD)+1])
            
            subplot(2,3,6,'replace'),  hold on
            plot(cast_time,sal,'k.-')
            plot(cast_time,depth,'.-')
            plot(cast_timeD,salD,'.','color',[0 0.5 1])
            ylabel('Salinity') %ylabel(col_hdr{6})
            
            %okay, now if everything looks good, record the depth where we
            %think some stratification is happening:
            keyboard
            %mld=[mld; {j time binned_data(im) N2(im)} {'down to'}];
            
        end
        
        
        if any(diff(cast_time) < 0), disp('Something wrong with time sync?'), end
        
        
    end
    
end

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
        
        jj=find(depthD>=4);
        if ~isempty(jj)
            mld{q,6}=temperatureD(jj(1));
            mld{q,7}=pdensD(jj(1));
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
