%Okay, before can use CTD casts, should go through the data to:

%flag or remove any bad data 
(%remove profiles with no data (less than 10 data points in vertical) or
%just bad cast (temp < -10))
%and
%extract down cast

%Now, Steve Lentz had already gone through manually identified top and
%bottom of downcast, but we can double check to see where his markings
%ended up

%Total CTDs (raw):
load '/Volumes/Lab_data/MVCO/processed_CTD_casts/CTD_30Mar2018.mat';

%Steve's data for comparison:
load('/Users/kristenhunter-cevera/MVCO_light_at_depth/mixed_layer_depth_src/processing_from_SteveLentz/mvco_ship_ctd_qc.mat')

%Question of binning data in a mat file? Might be useful...
z_grid=0:0.2:38; %to map data points onto...

%%

CTD_QC=CTD;
downcast_temp=nan(length(z_grid),length(CTD)); 
downcast_sal=downcast_temp; 
downcast_press=downcast_temp; 
downcast_pdens=downcast_temp;
downcast_bins=downcast_temp;
downcast_files={};
downcast_lat=[];
downcast_lon=[];

%%
for q=[154 155 160] %1:length(CTD); [1 10 20 51 53 107 108 150 151 152 154 155 160]
    
    col_hdr=CTD(q).data_hdr;
    temp_data=CTD(q).data;
    
    if any(~cellfun('isempty',temp_data)) %any data?
        
        %find and unload desired variables from cell array:
        %column indexes
        dp=find(cellfun('isempty',regexp(col_hdr,'Depth'))==0);
        s=find(cellfun('isempty',regexp(col_hdr,'Salinity'))==0);
        t=find(cellfun('isempty',regexp(col_hdr,'Temperature'))==0);
        p=find(cellfun('isempty',regexp(col_hdr,'Pressure'))==0);
        de=find(cellfun('isempty',regexp(col_hdr,'[^(Potential] Density'))==0);
        pd=find(cellfun('isempty',regexp(col_hdr,'Potential'))==0);
        ct=find(cellfun('isempty',regexp(col_hdr,'Time, Elapsed'))==0);
        
        depth=temp_data{dp}; %take off the first 50 data points?
        press=temp_data{p};
        temperature=temp_data{t};
        sal=temp_data{s};
        dens=temp_data{de};
        pdens=temp_data{pd};
        time=temp_data{ct};
        
        file_time=CTD(q).upload_time;
        
        %some basic QC checks:
        if max(depth) < 4
            CTD_QC(q).flag='not a cast';
            
        else %continue with QC!
            
            %find descent automatically:
            smooth_depth=smooth(depth,200); %just a bit easier for direct logic tests
            dc=find(diff(smooth_depth(1:end)) > 0.0025); %falling
                        
            %find longest continuous stretch-this should be the descent:
            nonconsec=find(diff(dc)~=1); %find the values that are not contiguous in the index
            [mm, ind]=max(diff(nonconsec)); %[mm, ind]=max(diff(dc(nonconsec)));
            dsc_ind1=dc(nonconsec(ind)+1);
            dsc_ind2=dc(nonconsec(ind+1));
            dsc=dsc_ind1:dsc_ind2; %indexes for main descent
            usc=dsc_ind2+1:length(depth); %and the upcast!
            
            %any bad salinity points?

            %a few products: binned values and interpolated values onto a grid:
            
            %bin variables into 0.2 m bins:
            max_depth=max(depth(dsc));
            min_depth=min(depth(dsc));
            
            depth_bins=floor(min_depth):0.2:ceil(max_depth);          
            for i=1:length(depth_bins)-1
                ii=find(depth(dsc) >= depth_bins(i) & depth(dsc) < depth_bins(i+1));                
                %binned_depth(i)=nanmean(depth(dsc(ii)));
                binned_sal(i)=nanmean(sal(dsc(ii)));
                binned_temp(i)=nanmean(temperature(dsc(ii)));
                binned_press(i)=nanmean(press(dsc(ii)));
                binned_pdens(i)=nanmean(pdens(dsc(ii)));
            end
            
            %move these into a grid matrix:
            zz=find(z_grid==depth_bins(1));
            downcast_temp(zz:(zz+length(depth_bins)-2),q)=binned_temp; 
            downcast_sal(zz:(zz+length(depth_bins)-2),q)=binned_sal; 
            downcast_press(zz:(zz+length(depth_bins)-2),q)=binned_press; 
            downcast_pdens(zz:(zz+length(depth_bins)-2),q)=binned_pdens;
            downcast_bins(zz:(zz+length(depth_bins)-1),q)=depth_bins;
            
            downcast_lat(q)=CTD(q).lat;
            downcast_lon(q)=CTD(q).lon;
            downcast_files{q}=CTD(q).cast_name;
            
            CTD_QC(q).downcast_ind = dsc;
            CTD_QC(q).upcast_ind = usc; %accounts for datapoints removed at start
            
            %If want to manually go for the points:
%             figure, plot(depth,'.'), set(gca,'YDir','reverse') %to get indices
%             [x, ~]=ginput(2); 
%             dsc_ind1=floor(x(1)); dsc_ind2=ceil(x(2));
    
            %I don't know if I want a grid...curious interpolation
%             %map to grid and to structure file:
%             dscB=dsc;
%             for ki=1:50; %the for loop is to repeat this process until
%             you have unqiue data points (and 50 seems to be enough to do
%             it...)
%                dscB=dscB(diff(depth(dscB))>0); %needs to have unqiue values for interpolation:
%             end %i think
%             
%             press_grid(:,q)=interp1(depth(dscB),press(dscB),z_grid);
%             temperature_grid(:,q)=interp1(depth(dscB),temperature(dscB),z_grid);
%             sal_grid(:,q)=interp1(depth(dscB),sal(dscB),z_grid);
%             dens_grid(:,q)=interp1(depth(dscB),dens(dscB),z_grid);
%             pdens_grid(:,q)=interp1(depth(dscB),pdens(dscB),z_grid);
%             time_grid(q)=file_time;
%             filename=[filename; {CTD{q}.cast_name}];
%             
            %and plot:
            %compare to Steve's manual points:
            jj=find(mdayul_ctd==file_time);
            if ~isempty(jj)
                ind_top=jt_ctd(jj); ind_bottom=jb_ctd(jj);
                steve_depth=z_ctd; steve_temp=wt_ctd(jj,:);
            else
                ind_top=1; ind_bottom=length(depth);
                steve_depth=NaN; steve_temp=NaN;
            end
            
            clf %see what metrics are doing - over all looks pretty good!
            
            subplot(2,2,1,'replace'), hold on
            plot(time,depth,'k.-')
            plot(time(dsc),depth(dsc),'.','color',[0 0.5 1])
            plot(time(usc),depth(usc),'.','color',[1 0.5 0])
            line([time(ind_top) time(ind_top)],ylim)
            line([time(ind_bottom) time(ind_bottom)],ylim)
            set(gca,'ydir','reverse','fontsize',14)
            xlabel('Time') %ylabel(col_hdr{6})
            title('CTD position with time')
            
            subplot(2,2,2,'replace'),  hold on
            %plot(pdens,depth,'k.-')
            plot(pdens(dsc),depth(dsc),'.','color',[0 0.5 1])
            plot(pdens(usc),depth(usc),'.','color',[1 0.5 0])
            set(gca,'ydir','reverse','fontsize',14)
            xlabel('Potential Density (kg/m^3)') %ylabel(col_hdr{6})
            title([CTD(q).cast_name ':' datestr(file_time) ' : ' num2str(q) ' out of ' num2str(length(CTD))],'interpreter','none')
            %legend('Obs \rho','Potential \rho','location','NorthWest')
            
            subplot(2,2,3,'replace'),  hold on
            %plot(temperature,depth,'k.-')
            plot(temperature(dsc),depth(dsc),'.','color',[0 0.5 1])
            plot(temperature(usc),depth(usc),'.','color',[1 0.5 0])
            plot(binned_temp,depth_bins(1:end-1),'ro-')
            plot(steve_temp,steve_depth,'k.-')
            set(gca,'ydir','reverse','fontsize',14)
            xlabel('Temperature (\circC)') %ylabel(col_hdr{6})
            line(xlim,[4 4],'color',[0.5 0.5 0.5])
            title('Temperature with Depth')
            
            subplot(2,2,4,'replace'),  hold on
            %plot(sal,depth,'k.-')
            plot(sal(dsc),depth(dsc),'.','color',[0 0.5 1])
            plot(sal(usc),depth(usc),'.','color',[1 0.5 0])
            set(gca,'ydir','reverse','fontsize',14)
            xlabel('Salinity') %ylabel(col_hdr{6})
            line([mean(sal(dsc))-2*std(sal(dsc)) mean(sal(dsc))-2*std(sal(dsc))],ylim)
            line([mean(sal(dsc))+2*std(sal(dsc)) mean(sal(dsc))+2*std(sal(dsc))],ylim)
            line(xlim,[4 4],'color',[0.5 0.5 0.5])
            title('Salinity with Depth')
            
            keyboard  
            
        end
    else %is empty file
        CTD_QC(q).flag='not a cast';
    end
    
    clearvars -except CTD CTD_QC downcast_* q z_grid *_ctd %remove temporary variables from each iteration
    
end

%% a bit of clean up:

downcast_bins=downcast_bins(1:end-1,:); %remove extra bin

jj=find(downcast_bins==0); %replace zeros with nan's:
downcast_bins(jj)=NaN;

downcast_pdens(jj)=NaN;
downcast_press(jj)=NaN;
downcast_temp(jj)=NaN;
downcast_sal(jj)=NaN;

downcast_lat(downcast_lat==0)=NaN;
downcast_lon(downcast_lon==0)=NaN;

%% double check them datas!

J=[];
for j=1:length(CTD_QC)
   if ~isempty(CTD_QC(j).flag)
       J=[J;j];
   end  
end

%good data:
ii=setxor(J,1:length(CTD_QC));

%% double check on bins:
for j=1:size(downcast_bins,1)
    if round(sum(diff(downcast_bins(j,downcast_bins(j,:)~=0)))) ~= 0 
        keyboard
    end
end

%should not keyboard!

%% if need to remove bad data points:

%only came across two cases:
% qq=find(depth < -1);
% qq=find(temperature > 80);
% 
% ii=setxor(qq,1:length(depth));
% for j=1:size(temp_data,2)
%     temp_data{j}=temp_data{j}(ii,:); %exclude bad data point...
% end
% 
% CTD_QC(q).data=temp_data;

