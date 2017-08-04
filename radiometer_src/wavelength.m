% Script to look at how the spectral quality of light changes with time:

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
wv_record=[];
plotflag=0;
extraplotflag=0;

for foldernum=[2:3 5:10 12 14:16 18:20]; %are these the ones with lat/lon?
    
    %find dat data!
    matsource=fullfile(sourcepath,datafolders{foldernum},'/mat_outfiles/');
    
    eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
    eval(['tempdata=data_' datafolders{foldernum} ';'])
    
    %     eval(['load ' matsource 'K_PAR_' datafolders{foldernum} '.mat'])
    %     eval(['K_PAR=K_PAR_' datafolders{foldernum} ';'])
    %
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
        lat=num2cell(nan(length(timestamp),1));%cell(length(timestamp),1);
        lon=lat;
        station=lat;
    end
    
    %now that we have the data in, let's look at those wavelenghts!
    
    %find the depth for 4m, get the wave distribution, then keep and find
    %max...
    for k=1:length(tempdata)
        
        if tempdata(k).emptyflag==0 %a legit cast, 1 is empty file, 2 is data, but not a cast
            wv=tempdata(k).wavelen_downwell;
            edl=tempdata(k).adj_edl;
            depth=tempdata(k).depth;
            edl_ind=tempdata(k).edl_ind;
            time=tempdata(k).mprtime;
            
            wv_sol=tempdata(k).wavelen_solarst;
            esl=tempdata(k).adj_esl;
            
            %Look only at the down cast:
            dz=smooth2(diff(depth),15); %change in depth, smoothed to reduce noise
            dc=find(dz > 0.025); %dc for downcast
            
            jj=find(depth(dc) > 2 & depth(dc) < max(depth)-1); %exclude top and bottom measurements
            
            %matching indexes for wavelengths:
            [~,impr,ipar]=intersect(dc(jj),edl_ind);
            impr=dc(jj(impr)); %indices for depth that belong to light measurements and downcast
            %ipar directly indexes into edl
            
            %record solar std spectrum:
            [~, im]=max(esl,[],2);
            max_wv=[0 mode(wv_sol(im))];
            
            i4=find(depth(impr,1)>3.5 & depth(impr,1)<4.5);
            i8=find(depth(impr,1)>7.5 & depth(impr,1)<8.5);
            i12=find(depth(impr,1)>11.5 & depth(impr,1)<12.5);
            
            [~, ii]=max(edl(ipar,:),[],2); %find max at each depth record
            max_wv=[max_wv; depth(impr) wv(ii)]; %record
            
            if isempty(i4), r4=[NaN NaN]; else r4=[nanmean(max_wv(i4+1,1)) mode(max_wv(i4+1,2))]; end
            if isempty(i8), r8=[NaN NaN]; else r8=[nanmean(max_wv(i8+1,1)) mode(max_wv(i8+1,2))]; end
            if isempty(i12), r12=[NaN NaN]; else r12=[nanmean(max_wv(i12+1,1)) mode(max_wv(i12+1,2))]; end
            
            wv_record=[wv_record; tempdata(k).timestamp lat{k} lon{k} max_wv(1,:) r4 r8 r12];
            
            if plotflag
                
                subplot(1,4,1,'replace')
                %a look at the cast:
                plot(time,depth,'-'), hold on
                plot(time(impr),depth(impr),'.','markersize',12)
                set(gca,'YDir','reverse')
                datetick
                
                subplot(1,4,2,'replace'), hold on
                plot(wv_sol,esl,'-','color',[0.5 0.5 0.5])
                title(tempdata(k).file,'Interpreter','none')
                
                for j=1:length(ipar)
                    scatter(wv,edl(ipar(j),:),40,depth(impr(j))*ones(size(wv)),'filled')
                end
                
                subplot(1,4,3,'replace'),hold on
                if ~isempty(i4), plot(wv,edl(ipar(i4),:),'-','color',[0 0.5 1]); end
                if ~isempty(i8),plot(wv,edl(ipar(i8),:),'-','color',[0 0 1]); end
                if ~isempty(i12),plot(wv,edl(ipar(i12),:),'-','color',[0 0 0.7]); end
                
                subplot(1,4,4,'replace'), hold on
                plot(max_wv(:,1),max_wv(:,2),'.','markersize',12)
                if ~isempty(i4), plot(wv_record(end,6),wv_record(end,7),'p','color',[0 0.5 1]); end
                if ~isempty(i8),plot(wv_record(end,8),wv_record(end,9),'p','color',[0 0 1]); end
                if ~isempty(i12),plot(wv_record(end,10),wv_record(end,11),'p','color',[0 0 0.7]); end
                
                xlabel('Depth (m)')
                ylabel('Wavelength of maximum energy (nm)')
                
                keyboard
            end
            
        end  
    end
end

%What I'd also like to do is record a set of casts, and then plot these
%against each other, color coded by time of year to identify any shifts??
%% find the tower casts:

box=[-70.7 -70.4 41.3 41.35];
ii=find(wv_record(:,3) > box(1) & wv_record(:,3) < box(2) & wv_record(:,2) > box(3) & wv_record(:,2) < box(4));

%%
figure
plot(find_yearday(wv_record(ii,1)),wv_record(ii,7),'.','markersize',12), hold on
plot(find_yearday(wv_record(ii,1)),wv_record(ii,9),'.','markersize',12)
plot(find_yearday(wv_record(ii,1)),wv_record(ii,11),'.','markersize',12)
